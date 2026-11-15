"""Train held-out evaluation heads and full-data deployment heads."""

from __future__ import annotations

import csv
import hashlib
import importlib.metadata
import json
from pathlib import Path

import joblib
import numpy as np
from sklearn.linear_model import RidgeClassifier
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    confusion_matrix,
    f1_score,
    recall_score,
)
from sklearn.model_selection import GroupShuffleSplit
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.utils.class_weight import compute_sample_weight
from xgboost import XGBClassifier

RANDOM_SEED = 42
NEUTRAL_THRESHOLD = 0.01
CLASS_NAMES = np.asarray(["decrease", "unchanged", "increase"])


def response_classes(values: np.ndarray) -> np.ndarray:
    """Convert corrected net growth to decrease/unchanged/increase labels."""
    return np.where(
        values < -NEUTRAL_THRESHOLD,
        0,
        np.where(values > NEUTRAL_THRESHOLD, 2, 1),
    ).astype(np.int32)


def _class_counts(labels: np.ndarray) -> dict[str, int]:
    counts = np.bincount(labels, minlength=len(CLASS_NAMES))
    return {name: int(count) for name, count in zip(CLASS_NAMES, counts, strict=True)}


def _classification_metrics(measured: np.ndarray, predictions: np.ndarray) -> dict[str, object]:
    opposite = ((measured == 0) & (predictions == 2)) | (
        (measured == 2) & (predictions == 0)
    )
    recalls = recall_score(
        measured, predictions, labels=np.arange(len(CLASS_NAMES)), average=None,
        zero_division=0,
    )
    return {
        "validation_rows": int(len(measured)),
        "correct_count": int((measured == predictions).sum()),
        "accuracy": float(accuracy_score(measured, predictions)),
        "balanced_accuracy": float(balanced_accuracy_score(measured, predictions)),
        "macro_f1": float(f1_score(measured, predictions, average="macro")),
        "class_recall": {
            name: float(value)
            for name, value in zip(CLASS_NAMES, recalls, strict=True)
        },
        "opposite_direction_count": int(opposite.sum()),
        "opposite_direction_percent": float(100.0 * opposite.mean()),
        "confusion_matrix": confusion_matrix(
            measured, predictions, labels=np.arange(len(CLASS_NAMES))
        ).tolist(),
    }


def _ridge() -> Pipeline:
    return Pipeline([
        ("scale", StandardScaler()),
        ("classifier", RidgeClassifier(alpha=1.0, class_weight="balanced")),
    ])


def _xgboost() -> XGBClassifier:
    return XGBClassifier(
        objective="multi:softprob",
        num_class=len(CLASS_NAMES),
        eval_metric="mlogloss",
        n_estimators=300,
        learning_rate=0.05,
        max_depth=4,
        min_child_weight=2,
        subsample=0.8,
        colsample_bytree=0.8,
        reg_lambda=1.0,
        random_state=RANDOM_SEED,
        n_jobs=-1,
    )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n")


def train_and_evaluate(
    embeddings_path: str | Path,
    output_dir: str | Path,
    reference_path: str | Path | None = None,
) -> dict:
    """Evaluate by held-out sites, then fit separate deployment heads on all rows."""
    embeddings_path = Path(embeddings_path)
    output_dir = Path(output_dir)
    evaluation_dir = output_dir / "evaluation"
    deployment_dir = output_dir / "deployment"
    evaluation_dir.mkdir(parents=True, exist_ok=True)
    deployment_dir.mkdir(parents=True, exist_ok=True)

    with np.load(embeddings_path, allow_pickle=False) as data:
        features = data["embeddings"].astype(np.float32)
        target = data["target"].astype(np.float32)
        positions = data["protein_start"].astype(np.int32)
        replicate_1 = data["replicate_1"].astype(np.float32)
        replicate_2 = data["replicate_2"].astype(np.float32)
        mutant = data["mutant"].astype(str)
        ref_aa = data["ref_aa"].astype(str)
        alt_aa = data["alt_aa"].astype(str)

    expected_rows = len(target)
    if features.ndim != 2 or features.shape != (expected_rows, 480):
        raise ValueError(f"expected ({expected_rows}, 480) features, got {features.shape}")
    if not np.isfinite(features).all() or not np.isfinite(target).all():
        raise ValueError("embeddings or targets contain non-finite values")
    if len(set(mutant.tolist())) != expected_rows:
        raise ValueError("mutant identifiers must be unique")

    labels = response_classes(target)
    splitter = GroupShuffleSplit(n_splits=1, test_size=0.20, random_state=RANDOM_SEED)
    train_index, validation_index = next(splitter.split(features, labels, groups=positions))
    train_positions = set(positions[train_index].tolist())
    validation_positions = set(positions[validation_index].tolist())
    if train_positions & validation_positions:
        raise RuntimeError("position leakage between splits")

    x_train, y_train = features[train_index], labels[train_index]
    x_validation, y_validation = features[validation_index], labels[validation_index]

    ridge_eval = _ridge()
    xgb_eval = _xgboost()
    ridge_eval.fit(x_train, y_train)
    xgb_eval.fit(
        x_train, y_train,
        sample_weight=compute_sample_weight(class_weight="balanced", y=y_train),
    )

    majority_class = int(np.bincount(y_train).argmax())
    majority_predictions = np.full(len(validation_index), majority_class, dtype=np.int32)
    ridge_predictions = ridge_eval.predict(x_validation).astype(np.int32)
    xgb_predictions = xgb_eval.predict(x_validation).astype(np.int32)
    xgb_probabilities = xgb_eval.predict_proba(x_validation)

    metrics = {
        "task": "three-class direction prediction",
        "class_names": CLASS_NAMES.tolist(),
        "class_definition": {
            "decrease": f"netgrowth_mean_corrected < -{NEUTRAL_THRESHOLD}",
            "unchanged": (
                f"-{NEUTRAL_THRESHOLD} <= netgrowth_mean_corrected "
                f"<= {NEUTRAL_THRESHOLD}"
            ),
            "increase": f"netgrowth_mean_corrected > {NEUTRAL_THRESHOLD}",
        },
        "split": {
            "random_seed": RANDOM_SEED,
            "strategy": "80/20 GroupShuffleSplit grouped by protein_start",
            "training_rows": int(len(train_index)),
            "validation_rows": int(len(validation_index)),
            "training_positions": len(train_positions),
            "validation_positions": len(validation_positions),
            "training_position_values": sorted(train_positions),
            "validation_position_values": sorted(validation_positions),
            "training_class_counts": _class_counts(y_train),
            "validation_class_counts": _class_counts(y_validation),
        },
        "primary_metric": "balanced_accuracy",
        "models": {
            "majority_baseline": _classification_metrics(y_validation, majority_predictions),
            "ridge": _classification_metrics(y_validation, ridge_predictions),
            "xgboost": _classification_metrics(y_validation, xgb_predictions),
        },
    }
    _write_json(evaluation_dir / "metrics.json", metrics)
    _write_json(output_dir / "metrics.json", metrics)

    fieldnames = [
        "mutant", "protein_start", "ref_aa", "alt_aa", "measured_netgrowth",
        "replicate_1", "replicate_2", "measured_class", "majority_prediction",
        "majority_correct", "ridge_prediction", "ridge_correct",
        "xgboost_prediction", "xgboost_correct", "xgboost_probability_decrease",
        "xgboost_probability_unchanged", "xgboost_probability_increase",
    ]
    rows: list[dict[str, object]] = []
    for local_row, source_row in enumerate(validation_index):
        rows.append({
            "mutant": mutant[source_row],
            "protein_start": int(positions[source_row]),
            "ref_aa": ref_aa[source_row],
            "alt_aa": alt_aa[source_row],
            "measured_netgrowth": float(target[source_row]),
            "replicate_1": float(replicate_1[source_row]),
            "replicate_2": float(replicate_2[source_row]),
            "measured_class": CLASS_NAMES[y_validation[local_row]],
            "majority_prediction": CLASS_NAMES[majority_predictions[local_row]],
            "majority_correct": bool(majority_predictions[local_row] == y_validation[local_row]),
            "ridge_prediction": CLASS_NAMES[ridge_predictions[local_row]],
            "ridge_correct": bool(ridge_predictions[local_row] == y_validation[local_row]),
            "xgboost_prediction": CLASS_NAMES[xgb_predictions[local_row]],
            "xgboost_correct": bool(xgb_predictions[local_row] == y_validation[local_row]),
            "xgboost_probability_decrease": float(xgb_probabilities[local_row, 0]),
            "xgboost_probability_unchanged": float(xgb_probabilities[local_row, 1]),
            "xgboost_probability_increase": float(xgb_probabilities[local_row, 2]),
        })
    for path in (evaluation_dir / "validation_predictions.csv", output_dir / "validation_predictions.csv"):
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

    joblib.dump(ridge_eval, evaluation_dir / "ridge_model.joblib")
    xgb_eval.save_model(evaluation_dir / "xgboost_model.json")

    ridge_full = _ridge()
    xgb_full = _xgboost()
    ridge_full.fit(features, labels)
    xgb_full.fit(
        features, labels,
        sample_weight=compute_sample_weight(class_weight="balanced", y=labels),
    )
    ridge_full_path = deployment_dir / "ridge_full.joblib"
    xgb_full_path = deployment_dir / "xgboost_full.json"
    joblib.dump(ridge_full, ridge_full_path)
    xgb_full.save_model(xgb_full_path)

    ridge_full_predictions = ridge_full.predict(features).astype(np.int32)
    xgb_full_predictions = xgb_full.predict(features).astype(np.int32)
    xgb_full_probabilities = xgb_full.predict_proba(features)
    validation_mask = np.isin(positions, np.asarray(sorted(validation_positions)))

    deployment_fields = [
        "mutant", "protein_start", "ref_aa", "alt_aa", "measured_netgrowth",
        "measured_class", "split", "ridge_prediction", "ridge_correct",
        "xgboost_prediction", "xgboost_correct", "xgboost_probability_decrease",
        "xgboost_probability_unchanged", "xgboost_probability_increase",
    ]
    with (deployment_dir / "all_predictions.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=deployment_fields)
        writer.writeheader()
        for row in range(expected_rows):
            writer.writerow({
                "mutant": mutant[row],
                "protein_start": int(positions[row]),
                "ref_aa": ref_aa[row],
                "alt_aa": alt_aa[row],
                "measured_netgrowth": float(target[row]),
                "measured_class": CLASS_NAMES[labels[row]],
                "split": "validation" if validation_mask[row] else "training",
                "ridge_prediction": CLASS_NAMES[ridge_full_predictions[row]],
                "ridge_correct": bool(ridge_full_predictions[row] == labels[row]),
                "xgboost_prediction": CLASS_NAMES[xgb_full_predictions[row]],
                "xgboost_correct": bool(xgb_full_predictions[row] == labels[row]),
                "xgboost_probability_decrease": float(xgb_full_probabilities[row, 0]),
                "xgboost_probability_unchanged": float(xgb_full_probabilities[row, 1]),
                "xgboost_probability_increase": float(xgb_full_probabilities[row, 2]),
            })

    reference = Path(reference_path) if reference_path is not None else None
    manifest = {
        "task": metrics["task"],
        "row_count": expected_rows,
        "embedding_dimensions": int(features.shape[1]),
        "neutral_threshold": NEUTRAL_THRESHOLD,
        "class_names": CLASS_NAMES.tolist(),
        "deployment_training_scope": "all 4,925 experimental rows",
        "evaluation_scope": "site-grouped held-out split; deployment predictions are not evaluation",
        "embeddings_sha256": _sha256(embeddings_path),
        "reference_sha256": _sha256(reference) if reference and reference.exists() else None,
        "ridge_model_sha256": _sha256(ridge_full_path),
        "xgboost_model_sha256": _sha256(xgb_full_path),
        "versions": {
            package: importlib.metadata.version(package)
            for package in ("numpy", "scikit-learn", "xgboost", "joblib")
        },
    }
    _write_json(deployment_dir / "manifest.json", manifest)
    return metrics
