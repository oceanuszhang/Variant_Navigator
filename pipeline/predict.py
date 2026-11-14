"""ABL1/imatinib measurement lookup, legacy held-out estimates, and runtime inference.

Measurements take precedence. Saved ESM test estimates are exposed separately.
Runtime RandomForest inference does not inherit the ESM validation metrics.
"""

from __future__ import annotations

from typing import Any
from pathlib import Path
import json
import hashlib

from pipeline import dms, features

# Fitted once per process; fitting takes a couple of seconds.
_MODEL: Any = None
_SCALER: Any = None
_MEASURED: dict[str, Any] | None = None

REFERENCE = Path(__file__).resolve().parents[1] / "data/cancer-exp/ABL1_NM_005157.6_DMS_reference_aa.txt"
VALIDATION = Path(__file__).resolve().parents[1] / "output/model/validation_predictions.json"


def held_out(mutant: str) -> dict[str, Any] | None:
    """Return a saved test estimate, separately from assay lookup and inference."""
    if not VALIDATION.exists():
        return None
    raw = VALIDATION.read_bytes()
    artifact = json.loads(raw)
    row = next((r for r in artifact["rows"] if r["mutant"] == mutant.upper()), None)
    if row is None:
        return None
    return {**row, "kind": "held_out", "model": artifact["model"],
            "head": artifact["head"], "seed": artifact["seed"],
            "metrics": artifact["metrics"][artifact["head"]],
            "artifact_sha256": hashlib.sha256(raw).hexdigest(),
            "split": "position-held-out", "provenance_limit":
            "Legacy artifact does not include training membership or sequence hash"}


class PredictError(RuntimeError):
    """Model unavailable, or the variant cannot be featurized."""


def _measured_index() -> dict[str, Any]:
    global _MEASURED
    if _MEASURED is None:
        _MEASURED = {v.mutant.upper(): v for v in dms.load_netgrowth()}
    return _MEASURED


def measured(mutant: str) -> dict[str, Any] | None:
    """The real measurement, if this variant is in the screen."""
    v = _measured_index().get(mutant.strip().upper())
    if v is None:
        return None
    return {
        "mutant": v.mutant,
        "measured_netgrowth": v.response,
        "replicate_spread": v.replicate_spread,
        "source": "measured (BCR-ABL DMS, 1200 nM imatinib)",
    }


def _fit() -> tuple[Any, Any]:
    global _MODEL, _SCALER
    if _MODEL is not None:
        return _MODEL, _SCALER
    try:
        from sklearn.ensemble import RandomForestRegressor
        from sklearn.preprocessing import StandardScaler
    except ImportError as exc:
        raise PredictError(f"scikit-learn unavailable: {exc}") from exc
    vs = dms.load_netgrowth()
    X, y, _groups, _names = features.featurize(vs, blocks="all")
    sc = StandardScaler().fit(X)
    m = RandomForestRegressor(
        n_estimators=300, min_samples_leaf=2, n_jobs=-1, random_state=0
    ).fit(sc.transform(X), y)
    _MODEL, _SCALER = m, sc
    return m, sc


def predict_variant(
    position: int, ref_aa: str, alt_aa: str, *, compound: str = "imatinib",
    allow_inference: bool = True
) -> dict[str, Any]:
    """Predicted net growth under 1200 nM imatinib for one substitution.

    Higher = more drug escape. Returns the measurement instead when the variant
    is already in the screen.
    """
    ref_aa, alt_aa = ref_aa.upper(), alt_aa.upper()
    seq = REFERENCE.read_text().strip()
    if compound.lower() != "imatinib":
        raise PredictError("Response model supports imatinib only")
    if not 242 <= position <= 512:
        raise PredictError("Response model supports ABL1 positions 242-512 only")
    if len(ref_aa) != 1 or len(alt_aa) != 1 or ref_aa not in features.AA or alt_aa not in features.AA or ref_aa == alt_aa:
        raise PredictError("Expected a non-synonymous standard amino-acid substitution")
    if seq[position - 1] != ref_aa:
        raise PredictError(f"ABL1 position {position} is {seq[position - 1]}, not {ref_aa}")
    mutant = f"{ref_aa}{position}{alt_aa}"
    hit = measured(mutant)
    if hit:
        return {**hit, "predicted": False, "kind": "measurement", "compound": "imatinib", "dose_nm": 1200}

    if not allow_inference:
        raise PredictError("No measured assay value; inference disabled in offline mode")
    model, scaler = _fit()
    row = features.substitution_features(ref_aa, alt_aa) + features.position_features(
        position, compound=compound, ref=ref_aa
    )
    yhat = float(model.predict(scaler.transform([row]))[0])
    vs = dms.load_netgrowth()
    responses = sorted(v.response for v in vs)
    pct = sum(1 for r in responses if r < yhat) / len(responses) * 100
    return {
        "mutant": mutant,
        "predicted_netgrowth": round(yhat, 5),
        "percentile_vs_screen": round(pct, 1),
        "predicted": True,
        "source": f"model (RandomForest on {len(vs)} measured variants)",
        "kind": "inference",
        "compound": "imatinib", "dose_nm": 1200,
        "features": features.FEATURE_NAMES,
        "validation_status": "unverified for this runtime fit",
        "caveat": "Runtime handcrafted-feature RandomForest; ESM validation metrics do not apply. Research only.",
    }


def available() -> bool:
    """Assay lookup is available without optional model dependencies."""
    return REFERENCE.exists() and dms.available()


def _main() -> int:
    for pos, ref, alt in ((315, "T", "I"), (396, "H", "R"), (486, "F", "S"), (999, "A", "V")):
        try:
            r = predict_variant(pos, ref, alt)
            tag = "PREDICTED" if r["predicted"] else "MEASURED "
            val = r.get("predicted_netgrowth", r.get("measured_netgrowth"))
            print(f"  {tag} {r['mutant']:8} netgrowth={val:+.5f}")
        except Exception as exc:  # noqa: BLE001
            print(f"  {ref}{pos}{alt}: {exc}")
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
