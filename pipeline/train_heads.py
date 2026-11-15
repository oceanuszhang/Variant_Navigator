"""Train and evaluate three-class treatment-response heads."""

from __future__ import annotations

import csv
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


def _classification_metrics(
    measured: np.ndarray,
    predictions: np.ndarray,
) -> dict[str, object]:
    opposite = ((measured == 0) & (predictions == 2)) | (
        (measured == 2) & (predictions == 0)
    )
    return {
        "validation_rows": int(len(measured)),
        "correct_count": int((measured == predictions).sum()),
        "accuracy": float(accuracy_score(measured, predictions)),
        "balanced_accuracy": float(balanced_accuracy_score(measured, predictions)),
        "macro_f1": float(f1_score(measured, predictions, average="macro")),
        "opposite_direction_count": int(opposite.sum()),
        "opposite_direction_percent": float(100.0 * opposite.mean()),
        "confusion_matrix": confusion_matrix(
            measured,
            predictions,
            labels=np.arange(len(CLASS_NAMES)),
        ).tolist(),
    }


def train_and_evaluate(embeddings_path: str | Path, output_dir: str | Path) -> dict:
    embeddings_path = Path(embeddings_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

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

    labels = response_classes(target)
    splitter = GroupShuffleSplit(n_splits=1, test_size=0.20, random_state=RANDOM_SEED)
    train_index, validation_index = next(
        splitter.split(features, labels, groups=positions)
    )

    train_positions = set(positions[train_index].tolist())
    validation_positions = set(positions[validation_index].tolist())
    overlap = train_positions & validation_positions
    if overlap:
        raise RuntimeError(f"position leakage between splits: {sorted(overlap)}")

    x_train = features[train_index]
    y_train = labels[train_index]
    x_validation = features[validation_index]
    y_validation = labels[validation_index]

    ridge = Pipeline(
        [
            ("scale", StandardScaler()),
            (
                "classifier",
                RidgeClassifier(alpha=1.0, class_weight="balanced"),
            ),
        ]
    )
    xgboost = XGBClassifier(
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

    ridge.fit(x_train, y_train)
    xgboost.fit(
        x_train,
        y_train,
        sample_weight=compute_sample_weight(class_weight="balanced", y=y_train),
    )

    majority_class = int(np.bincount(y_train).argmax())
    majority_predictions = np.full(len(validation_index), majority_class, dtype=np.int32)
    ridge_predictions = ridge.predict(x_validation).astype(np.int32)
    xgboost_predictions = xgboost.predict(x_validation).astype(np.int32)
    xgboost_probabilities = xgboost.predict_proba(x_validation)

    metrics = {
        "task": "three-class direction prediction",
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
            "training_class_counts": _class_counts(y_train),
            "validation_class_counts": _class_counts(y_validation),
        },
        "primary_metric": "exact three-class accuracy",
        "models": {
            "majority_baseline": _classification_metrics(
                y_validation,
                majority_predictions,
            ),
            "ridge": _classification_metrics(y_validation, ridge_predictions),
            "xgboost": _classification_metrics(y_validation, xgboost_predictions),
        },
    }

    (output_dir / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")

    predictions_path = output_dir / "validation_predictions.csv"
    with predictions_path.open("w", newline="") as handle:
        fieldnames = [
            "mutant",
            "protein_start",
            "ref_aa",
            "alt_aa",
            "measured_netgrowth",
            "replicate_1",
            "replicate_2",
            "measured_class",
            "majority_prediction",
            "majority_correct",
            "ridge_prediction",
            "ridge_correct",
            "xgboost_prediction",
            "xgboost_correct",
            "xgboost_probability_decrease",
            "xgboost_probability_unchanged",
            "xgboost_probability_increase",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for local_row, source_row in enumerate(validation_index):
            writer.writerow(
                {
                    "mutant": mutant[source_row],
                    "protein_start": int(positions[source_row]),
                    "ref_aa": ref_aa[source_row],
                    "alt_aa": alt_aa[source_row],
                    "measured_netgrowth": float(target[source_row]),
                    "replicate_1": float(replicate_1[source_row]),
                    "replicate_2": float(replicate_2[source_row]),
                    "measured_class": CLASS_NAMES[y_validation[local_row]],
                    "majority_prediction": CLASS_NAMES[majority_predictions[local_row]],
                    "majority_correct": bool(
                        majority_predictions[local_row] == y_validation[local_row]
                    ),
                    "ridge_prediction": CLASS_NAMES[ridge_predictions[local_row]],
                    "ridge_correct": bool(
                        ridge_predictions[local_row] == y_validation[local_row]
                    ),
                    "xgboost_prediction": CLASS_NAMES[xgboost_predictions[local_row]],
                    "xgboost_correct": bool(
                        xgboost_predictions[local_row] == y_validation[local_row]
                    ),
                    "xgboost_probability_decrease": float(
                        xgboost_probabilities[local_row, 0]
                    ),
                    "xgboost_probability_unchanged": float(
                        xgboost_probabilities[local_row, 1]
                    ),
                    "xgboost_probability_increase": float(
                        xgboost_probabilities[local_row, 2]
                    ),
                }
            )

    joblib.dump(ridge, output_dir / "ridge_model.joblib")
    xgboost.save_model(output_dir / "xgboost_model.json")
    return metrics
