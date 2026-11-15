"""Trusted local access to saved ESM embeddings and classification outputs."""

from __future__ import annotations

import csv
import hashlib
import json
from functools import lru_cache
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "data/cancer-exp/ABL1_NM_005157.6_DMS_reference_aa.txt"
OUTPUT = ROOT / "output/cancer-exp"
DEPLOYMENT = OUTPUT / "deployment"
EVALUATION = OUTPUT / "evaluation"
VARIANT_RE = re.compile(r"([ACDEFGHIKLMNPQRSTVWY])(\d{3})([ACDEFGHIKLMNPQRSTVWY])")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def reference_sequence() -> str:
    return "".join(REFERENCE.read_text().split()).upper()


def normalize_sequence(value: str) -> str:
    lines = value.strip().splitlines()
    if lines and lines[0].startswith(">"):
        lines = lines[1:]
    sequence = "".join("".join(lines).split()).upper()
    if not sequence or not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", sequence):
        raise ValueError(
            "Paste one amino-acid sequence using the 20 standard residue letters; "
            "one FASTA header is allowed."
        )
    return sequence


def validate_input(variant: str, supplied_sequence: str) -> tuple[str, int, str, str]:
    variant = variant.strip().upper()
    match = VARIANT_RE.fullmatch(variant)
    if not match:
        raise ValueError("Enter one substitution such as E450R.")
    ref, position_text, alt = match.groups()
    position = int(position_text)
    sequence = reference_sequence()
    pasted = normalize_sequence(supplied_sequence)
    if pasted != sequence:
        raise ValueError(
            f"This demo supports the exact {len(sequence)}-aa ABL1 study reference only."
        )
    if not 1 <= position <= len(sequence):
        raise ValueError(f"Position must be between 1 and {len(sequence)}.")
    if sequence[position - 1] != ref:
        raise ValueError(
            f"ABL1 study reference position {position} is "
            f"{sequence[position - 1]}, not {ref}."
        )
    if ref == alt:
        raise ValueError("The alternate amino acid must differ from the reference.")
    return variant, position, ref, alt


@lru_cache(maxsize=1)
def _artifacts() -> dict:
    manifest_path = DEPLOYMENT / "manifest.json"
    predictions_path = DEPLOYMENT / "all_predictions.csv"
    validation_path = EVALUATION / "validation_predictions.csv"
    metrics_path = EVALUATION / "metrics.json"
    for path in (manifest_path, predictions_path, validation_path, metrics_path):
        if not path.exists():
            raise FileNotFoundError(
                f"Missing {path}. Run: python main.py --train-only"
            )
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("reference_sha256") != _sha256(REFERENCE):
        raise RuntimeError("Reference sequence does not match the deployment manifest.")
    if manifest.get("embeddings_sha256") != _sha256(OUTPUT / "embeddings.npz"):
        raise RuntimeError("Embedding artifact does not match the deployment manifest.")
    with predictions_path.open(newline="") as handle:
        deployment_rows = {row["mutant"]: row for row in csv.DictReader(handle)}
    with validation_path.open(newline="") as handle:
        validation_rows = {row["mutant"]: row for row in csv.DictReader(handle)}
    return {
        "manifest": manifest,
        "deployment": deployment_rows,
        "validation": validation_rows,
        "metrics": json.loads(metrics_path.read_text()),
    }


def clear_cache() -> None:
    _artifacts.cache_clear()


def predict_saved(variant: str, supplied_sequence: str, model: str = "xgboost") -> dict:
    variant, position, ref, alt = validate_input(variant, supplied_sequence)
    if model not in {"xgboost", "ridge"}:
        raise ValueError("Model must be xgboost or ridge.")
    artifacts = _artifacts()
    row = artifacts["deployment"].get(variant)
    if row is None:
        raise ValueError(
            f"{variant} is not one of the 4,925 experimentally measured substitutions."
        )
    validation = artifacts["validation"].get(variant)
    result = {
        "mutant": variant,
        "position": position,
        "reference_aa": ref,
        "alternate_aa": alt,
        "model": model,
        "prediction": row[f"{model}_prediction"],
        "measured_class": row["measured_class"],
        "measured_netgrowth": float(row["measured_netgrowth"]),
        "split": row["split"],
        "deployment_scope": (
            "Model refit on all 4,925 rows for the demo; this prediction is not "
            "an independent accuracy estimate."
        ),
        "held_out": None,
        "sequence_context": {
            "before": reference_sequence()[max(0, position - 5):position - 1],
            "reference": ref,
            "alternate": alt,
            "after": reference_sequence()[position:position + 4],
        },
    }
    if model == "xgboost":
        result["probabilities"] = {
            name: float(row[f"xgboost_probability_{name}"])
            for name in ("decrease", "unchanged", "increase")
        }
    if validation is not None:
        held = {
            "prediction": validation[f"{model}_prediction"],
            "correct": validation[f"{model}_correct"].lower() == "true",
            "scope": "Site-grouped held-out evaluation prediction",
        }
        if model == "xgboost":
            held["probabilities"] = {
                name: float(validation[f"xgboost_probability_{name}"])
                for name in ("decrease", "unchanged", "increase")
            }
        result["held_out"] = held
    return result


def validation_summary() -> dict:
    artifacts = _artifacts()
    metrics = artifacts["metrics"]
    rows = list(artifacts["validation"].values())
    return {
        "metrics": metrics,
        "rows": rows,
        "artifact_sha256": _sha256(EVALUATION / "validation_predictions.csv"),
        "manifest": artifacts["manifest"],
    }
