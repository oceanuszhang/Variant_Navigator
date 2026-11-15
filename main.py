#!/usr/bin/env python3
"""Run the minimal BCR-ABL imatinib growth-direction experiment."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

from pipeline.train_heads import train_and_evaluate

ROOT = Path(__file__).resolve().parent
DEFAULT_SEQUENCE = ROOT / "data/cancer-exp/ABL1_NM_005157.6_DMS_reference_aa.txt"
DEFAULT_DATA = ROOT / "data/cancer-exp/high_imatinib_netgrowth_data.csv"
DEFAULT_OUTPUT = ROOT / "output/cancer-exp"


def generate_embeddings(
    sequence_path: Path,
    data_path: Path,
    embeddings_path: Path,
    batch_size: int,
) -> None:
    command = [
        sys.executable,
        "-m",
        "modal",
        "run",
        str(ROOT / "pipeline/modal_esm.py"),
        "--sequence-path",
        str(sequence_path),
        "--csv-path",
        str(data_path),
        "--output-path",
        str(embeddings_path),
        "--batch-size",
        str(batch_size),
    ]
    subprocess.run(command, cwd=ROOT, check=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sequence", type=Path, default=DEFAULT_SEQUENCE)
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument(
        "--recompute-embeddings",
        action="store_true",
        help="rerun ESM-2 even when the cached embedding file exists",
    )
    parser.add_argument(
        "--train-only",
        action="store_true",
        help="train from an existing embeddings.npz without calling Modal",
    )
    args = parser.parse_args()

    sequence_path = args.sequence.resolve()
    data_path = args.data.resolve()
    output_dir = args.output_dir.resolve()
    embeddings_path = output_dir / "embeddings.npz"
    output_dir.mkdir(parents=True, exist_ok=True)

    if not sequence_path.exists():
        parser.error(f"reference sequence not found: {sequence_path}")
    if not data_path.exists():
        parser.error(f"experiment CSV not found: {data_path}")
    if args.train_only and not embeddings_path.exists():
        parser.error(f"--train-only requires {embeddings_path}")

    if not args.train_only and (args.recompute_embeddings or not embeddings_path.exists()):
        print("[1/2] Generating frozen ESM-2 embeddings on Modal", flush=True)
        generate_embeddings(sequence_path, data_path, embeddings_path, args.batch_size)
    else:
        print(f"[1/2] Reusing {embeddings_path}", flush=True)

    print("[2/2] Training Ridge and XGBoost direction classifiers", flush=True)
    metrics = train_and_evaluate(embeddings_path, output_dir, sequence_path)

    for model_name in ("majority_baseline", "ridge", "xgboost"):
        result = metrics["models"][model_name]
        print(
            f"{model_name}: "
            f"accuracy={100.0 * result['accuracy']:.2f}%; "
            f"balanced_accuracy={100.0 * result['balanced_accuracy']:.2f}%; "
            f"macro_F1={result['macro_f1']:.3f}"
        )
    print(f"Results written to {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
