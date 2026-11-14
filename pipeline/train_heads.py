"""Train and evaluate Ridge and XGBoost heads on frozen ESM-2 features."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import joblib
import numpy as np
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error
from sklearn.model_selection import GroupShuffleSplit
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from xgboost import XGBRegressor

RANDOM_SEED = 42


def _model_metrics(
    measured: np.ndarray,
    predictions: np.ndarray,
    replicate_1: np.ndarray,
    replicate_2: np.ndarray,
) -> tuple[dict[str, float | int], np.ndarray, np.ndarray]:
    lower = np.minimum(replicate_1, replicate_2)
    upper = np.maximum(replicate_1, replicate_2)
    in_range = (predictions >= lower) & (predictions <= upper)
    absolute_error = np.abs(predictions - measured)
    metrics: dict[str, float | int] = {
        "validation_rows": int(len(measured)),
        "in_range_count": int(in_range.sum()),
        "in_range_percent": float(100.0 * in_range.mean()),
        "mean_absolute_error": float(mean_absolute_error(measured, predictions)),
    }
    return metrics, in_range, absolute_error


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

    splitter = GroupShuffleSplit(n_splits=1, test_size=0.20, random_state=RANDOM_SEED)
    train_index, validation_index = next(splitter.split(features, target, groups=positions))

    train_positions = set(positions[train_index].tolist())
    validation_positions = set(positions[validation_index].tolist())
    overlap = train_positions & validation_positions
    if overlap:
        raise RuntimeError(f"position leakage between splits: {sorted(overlap)}")

    x_train = features[train_index]
    y_train = target[train_index]
    x_validation = features[validation_index]
    y_validation = target[validation_index]

    ridge = Pipeline(
        [
            ("scale", StandardScaler()),
            ("regressor", Ridge(alpha=1.0)),
        ]
    )
    xgboost = XGBRegressor(
        objective="reg:squarederror",
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
    xgboost.fit(x_train, y_train)

    ridge_predictions = ridge.predict(x_validation)
    xgboost_predictions = xgboost.predict(x_validation)

    ridge_metrics, ridge_in_range, ridge_error = _model_metrics(
        y_validation,
        ridge_predictions,
        replicate_1[validation_index],
        replicate_2[validation_index],
    )
    xgboost_metrics, xgboost_in_range, xgboost_error = _model_metrics(
        y_validation,
        xgboost_predictions,
        replicate_1[validation_index],
        replicate_2[validation_index],
    )

    metrics = {
        "split": {
            "random_seed": RANDOM_SEED,
            "strategy": "80/20 GroupShuffleSplit grouped by protein_start",
            "training_rows": int(len(train_index)),
            "validation_rows": int(len(validation_index)),
            "training_positions": len(train_positions),
            "validation_positions": len(validation_positions),
        },
        "range_definition": (
            "prediction is between netgr_corr_rep1.high and "
            "netgr_corr_rep2.high, inclusive"
        ),
        "models": {
            "ridge": ridge_metrics,
            "xgboost": xgboost_metrics,
        },
    }

    metrics_path = output_dir / "metrics.json"
    metrics_path.write_text(json.dumps(metrics, indent=2) + "\n")

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
            "replicate_lower",
            "replicate_upper",
            "ridge_prediction",
            "ridge_absolute_error",
            "ridge_in_range",
            "xgboost_prediction",
            "xgboost_absolute_error",
            "xgboost_in_range",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for local_row, source_row in enumerate(validation_index):
            low = min(replicate_1[source_row], replicate_2[source_row])
            high = max(replicate_1[source_row], replicate_2[source_row])
            writer.writerow(
                {
                    "mutant": mutant[source_row],
                    "protein_start": int(positions[source_row]),
                    "ref_aa": ref_aa[source_row],
                    "alt_aa": alt_aa[source_row],
                    "measured_netgrowth": float(target[source_row]),
                    "replicate_1": float(replicate_1[source_row]),
                    "replicate_2": float(replicate_2[source_row]),
                    "replicate_lower": float(low),
                    "replicate_upper": float(high),
                    "ridge_prediction": float(ridge_predictions[local_row]),
                    "ridge_absolute_error": float(ridge_error[local_row]),
                    "ridge_in_range": bool(ridge_in_range[local_row]),
                    "xgboost_prediction": float(xgboost_predictions[local_row]),
                    "xgboost_absolute_error": float(xgboost_error[local_row]),
                    "xgboost_in_range": bool(xgboost_in_range[local_row]),
                }
            )

    joblib.dump(ridge, output_dir / "ridge_model.joblib")
    xgboost.save_model(output_dir / "xgboost_model.json")
    return metrics
