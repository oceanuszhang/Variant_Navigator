"""Load trusted local ESM regression artifacts and validate their feature contract.

Joblib files can execute Python code: use only artifacts created by this project
or obtained from a trusted source. Hash checks detect corruption, not trust.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any


FEATURE_RECIPE = "mutant_minus_reference_at_variant_residue"


def load_run(run_dir: Path, *, allow_partial: bool = False) -> tuple[dict[str, Any], Any]:
    """Validate and load a trusted local artifact before spending GPU time."""
    import joblib
    import sklearn

    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    if manifest["schema_version"] != 1:
        raise ValueError("Unsupported model manifest schema")
    if manifest["partial"] and not allow_partial:
        raise ValueError("Partial smoke artifacts cannot be used as production models")
    if manifest["versions"]["sklearn"] != sklearn.__version__:
        raise ValueError("scikit-learn version differs from persisted model")
    if hashlib.sha256(manifest["reference"].encode()).hexdigest() != manifest["sequence_sha256"]:
        raise ValueError("Reference sequence checksum mismatch")
    if set(manifest["train_positions"]) & set(manifest["test_positions"]):
        raise ValueError("Training and test positions overlap")
    artifact = manifest["artifacts"][manifest["head"]]
    path = (run_dir / artifact["file"]).resolve()
    if path.parent != run_dir.resolve():
        raise ValueError("Model path escapes run directory")
    if hashlib.sha256(path.read_bytes()).hexdigest() != artifact["sha256"]:
        raise ValueError("Persisted model checksum mismatch")
    return manifest, joblib.load(path)


def _predict(manifest: dict[str, Any], fitted: Any, vector: list[float]) -> dict[str, Any]:
    import numpy as np
    array = np.asarray(vector, dtype=float)
    if array.shape != (manifest["feature_dim"],) or not np.isfinite(array).all():
        raise ValueError("Invalid embedding dimensions or nonfinite values")
    return {"kind": "inference", "predicted_netgrowth": float(fitted.predict(array[None, :])[0]),
            "model": manifest["model"], "head": manifest["head"], "partial": manifest["partial"],
            "validation_metrics": manifest["metrics"][manifest["head"]],
            "validation_scope": "saved position-held-out split; not a calibrated uncertainty estimate"}


def predict_embedding(run_dir: Path, vector: list[float], *, model: str,
                      sequence_sha256: str, feature_recipe: str,
                      allow_partial: bool = False) -> dict[str, Any]:
    """Predict from a feature vector with matching model and sequence provenance."""
    manifest, fitted = load_run(run_dir, allow_partial=allow_partial)
    for key, value in (("model", model), ("sequence_sha256", sequence_sha256),
                       ("feature_recipe", feature_recipe)):
        if manifest[key] != value:
            raise ValueError(f"Feature contract mismatch: {key}")
    return _predict(manifest, fitted, vector)


def predict_variant(run_dir: Path, variant: str, *, allow_partial: bool = False) -> dict[str, Any]:
    """Embed one ABL1 substitution on Modal and apply the exact saved head.

    This explicit model estimate is separate from the preferred measured lookup.
    Invalid variants/artifacts are rejected before importing or contacting Modal.
    """
    from pipeline.mutant import apply_variant
    manifest, fitted = load_run(run_dir, allow_partial=allow_partial)
    if manifest["feature_recipe"] != FEATURE_RECIPE:
        raise ValueError("Unsupported ESM feature recipe")
    if manifest["model"] not in {"esm2_t12_35M_UR50D", "esm2_t33_650M_UR50D"}:
        raise ValueError("Unsupported ESM model")
    lo, hi = manifest["crop"]
    reference = manifest["reference"]
    if lo != 1 or hi != len(reference):
        raise ValueError("Single-variant inference currently requires a reference starting at residue 1")
    mutant = apply_variant(reference, variant)
    if len(mutant.substitutions) != 1 or mutant.deletions:
        raise ValueError("ESM response inference requires one substitution")
    sub = mutant.substitutions[0]
    if not 242 <= sub.position <= 512 or sub.wt == sub.mt:
        raise ValueError("ESM response model supports non-synonymous ABL1 substitutions at 242-512")
    from scripts.modal_esm import app, embed_mutants
    try:
        with app.run():
            vectors = embed_mutants.remote(reference, [(sub.position, sub.mt)], lo, manifest["model"], 1)
    except Exception as exc:
        raise RuntimeError(f"Modal ESM inference failed: {exc}") from exc
    if len(vectors) != 1:
        raise ValueError("Unexpected number of ESM feature vectors")
    result = _predict(manifest, fitted, vectors[0])
    membership = ("train" if sub.position in manifest["train_positions"] else
                  "test" if sub.position in manifest["test_positions"] else "unseen_position")
    return {**result, "mutant": str(sub), "split_membership": membership,
            "compound": "imatinib", "dose_nm": 1200,
            "source": "persisted ESM head with Modal embedding",
            "sequence_sha256": manifest["sequence_sha256"],
            "model_sha256": manifest["artifacts"][manifest["head"]]["sha256"],
            "caveat": "Research model estimate; training-position results are not held-out validation."}
