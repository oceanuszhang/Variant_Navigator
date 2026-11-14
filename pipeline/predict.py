"""Trained net-growth prediction for a single variant — the demo path.

Deliberately small. It fits the model from `scripts/train_baseline.py` on the
full BCR-ABL screen and predicts one variant, so `main.py` can answer "what
happens to this variant under imatinib?" with a number rather than only a
geometric prior.

Honesty carried with every prediction, because the number is modest:

    r = 0.245 under position-grouped CV, which is 25% of the attainable
    ceiling (0.978), with MAE 2.8x the replicate noise floor.

So a prediction is reported with its measured accuracy attached and is useful
for ranking, not for a clinical threshold. `measured()` short-circuits the whole
thing when the variant is already in the screen: 4,925 variants have a real
measurement, and a measurement beats a prediction.
"""

from __future__ import annotations

from typing import Any

from pipeline import dms, features

# Fitted once per process; fitting takes a couple of seconds.
_MODEL: Any = None
_SCALER: Any = None
_MEASURED: dict[str, Any] | None = None

# From scripts/train_baseline.py, 5-fold position-grouped CV.
CV_PEARSON_R = 0.2454
CV_MAE = 0.01955
CEILING_R = 0.9779


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
    position: int, ref_aa: str, alt_aa: str, *, compound: str = "imatinib"
) -> dict[str, Any]:
    """Predicted net growth under 1200 nM imatinib for one substitution.

    Higher = more drug escape. Returns the measurement instead when the variant
    is already in the screen.
    """
    mutant = f"{ref_aa.upper()}{position}{alt_aa.upper()}"
    hit = measured(mutant)
    if hit:
        return {**hit, "predicted": False, "cv_pearson_r": CV_PEARSON_R}

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
        "cv_pearson_r": CV_PEARSON_R,
        "cv_mae": CV_MAE,
        "ceiling_r": CEILING_R,
        "caveat": (
            f"r={CV_PEARSON_R} under position-grouped CV ({CV_PEARSON_R/CEILING_R*100:.0f}% "
            "of attainable). Useful for ranking, not as a clinical threshold."
        ),
    }


def available() -> bool:
    try:
        import sklearn  # noqa: F401
    except ImportError:
        return False
    return dms.available()


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
