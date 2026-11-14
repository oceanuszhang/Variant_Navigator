#!/usr/bin/env python3
"""Train and evaluate on ESM-2 mutant-marginal embeddings.

    .venv/bin/python scripts/train_embeddings.py

Two heads (ridge, gradient boosting) on the 480-dim embedding differences from
`embed_variants.py`, with an 80/20 split GROUPED BY RESIDUE POSITION so that
substitutions at the same position cannot land on both sides. Fixed seed.

Metrics:
  MAE, Spearman, Pearson        against the measured netgrowth_mean_corrected
  in-range %                    prediction falls between the variant's two
                                corrected replicates

READ THE IN-RANGE NUMBER CAREFULLY. The replicate interval is narrow - the mean
gap between replicates is 0.0069 while the response spans ~0.12 - so landing
inside it is a demanding test, not a generous one. A low percentage does not
mean the model is broken; it means the model is not yet as precise as repeating
the experiment. Reported alongside the width of the interval so the two are
never confused.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

ESM_DIR = Path("data/esm")
OUT_DIR = Path("output/model")


def spearman(a, b):
    import numpy as np
    def rank(v):
        order = np.argsort(v, kind="mergesort")
        r = np.empty(len(v), dtype=float)
        r[order] = np.arange(1, len(v) + 1)
        return r
    ra, rb = rank(np.asarray(a)), rank(np.asarray(b))
    ra -= ra.mean(); rb -= rb.mean()
    denom = (np.sqrt((ra**2).sum()) * np.sqrt((rb**2).sum()))
    return float((ra * rb).sum() / denom) if denom else 0.0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--npz", default=None)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--test-frac", type=float, default=0.2)
    args = ap.parse_args()

    import numpy as np
    from sklearn.ensemble import HistGradientBoostingRegressor
    from sklearn.linear_model import RidgeCV

    path = Path(args.npz) if args.npz else next(ESM_DIR.glob("*_mutdiff.npz"))
    d = np.load(path, allow_pickle=True)
    X, y, groups = d["X"], d["y"], d["groups"]
    rep1, rep2, names = d["rep1"], d["rep2"], d["names"]
    model_name = str(d["model"][0]) if "model" in d else "?"
    print(f"{path.name}: X={X.shape}  model={model_name}")

    # Position-grouped 80/20 split, deterministic.
    positions = np.array(sorted(set(groups.tolist())))
    rng = np.random.default_rng(args.seed)
    perm = rng.permutation(len(positions))
    n_test = max(1, int(round(len(positions) * args.test_frac)))
    test_pos = set(positions[perm[:n_test]].tolist())
    te = np.array([g in test_pos for g in groups])
    tr = ~te
    print(f"split: train={tr.sum()} ({len(set(groups[tr].tolist()))} positions)  "
          f"test={te.sum()} ({len(test_pos)} positions)  seed={args.seed}")
    assert not (set(groups[tr].tolist()) & test_pos), "position leaked across the split"
    print("leakage check: 0 shared positions")

    heads = {
        "ridge": RidgeCV(alphas=[0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0]),
        "hist_gbm": HistGradientBoostingRegressor(
            max_iter=400, learning_rate=0.06, max_depth=6,
            min_samples_leaf=20, l2_regularization=1.0, random_state=args.seed,
        ),
    }

    # Replicate interval, for the in-range metric and for context.
    has_rep = ~(np.isnan(rep1) | np.isnan(rep2))
    width = np.abs(rep1 - rep2)
    print(f"\nreplicate interval: mean width {np.nanmean(width):.5f}, "
          f"median {np.nanmedian(width):.5f}, response sd {y.std():.5f}")
    print("  (a prediction must land inside this to count as in-range)")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    results = {}
    best = None
    print(f"\n{'head':10} {'MAE':>9} {'Spearman':>9} {'Pearson':>9} {'in-range':>9}")
    print("-" * 52)
    for name, head in heads.items():
        head.fit(X[tr], y[tr])
        pred = head.predict(X[te])
        mae = float(np.mean(np.abs(pred - y[te])))
        sp = spearman(pred, y[te])
        pe = float(np.corrcoef(pred, y[te])[0, 1])
        m = has_rep[te]
        lo = np.minimum(rep1[te][m], rep2[te][m])
        hi = np.maximum(rep1[te][m], rep2[te][m])
        inr = float(np.mean((pred[m] >= lo) & (pred[m] <= hi))) * 100
        results[name] = {"mae": mae, "spearman": sp, "pearson": pe,
                         "in_range_pct": inr, "n_test": int(te.sum())}
        print(f"{name:10} {mae:>9.5f} {sp:>9.4f} {pe:>9.4f} {inr:>8.1f}%")
        if best is None or sp > results[best]["spearman"]:
            best = name

    # Trivial baseline: predict the training mean.
    base = float(np.mean(y[tr]))
    bmae = float(np.mean(np.abs(base - y[te])))
    print(f"{'mean':10} {bmae:>9.5f} {0.0:>9.4f} {0.0:>9.4f} {'-':>9}")
    print("-" * 52)
    print(f"best head: {best} (Spearman {results[best]['spearman']:.4f}, "
          f"MAE {results[best]['mae']:.5f} vs {bmae:.5f} for predict-the-mean)")

    # Per-variant table for the demo, from the best head.
    head = heads[best]
    pred = head.predict(X[te])
    rows = []
    for k, i in enumerate(np.where(te)[0]):
        lo_i = min(rep1[i], rep2[i]) if has_rep[i] else None
        hi_i = max(rep1[i], rep2[i]) if has_rep[i] else None
        rows.append({
            "mutant": str(names[i]),
            "position": int(groups[i]),
            "measured": float(y[i]),
            "predicted": float(pred[k]),
            "error": float(pred[k] - y[i]),
            "abs_error": abs(float(pred[k] - y[i])),
            "rep_low": None if lo_i is None else float(lo_i),
            "rep_high": None if hi_i is None else float(hi_i),
            "in_range": bool(lo_i is not None and lo_i <= pred[k] <= hi_i),
        })
    rows.sort(key=lambda r: r["abs_error"])
    (OUT_DIR / "validation_predictions.json").write_text(json.dumps({
        "model": model_name, "head": best, "seed": args.seed,
        "metrics": results, "mean_baseline_mae": bmae, "rows": rows,
    }, indent=1))
    print(f"\nwrote {OUT_DIR/'validation_predictions.json'} ({len(rows)} validation variants)")

    print("\nbest 5 predictions:")
    for r in rows[:5]:
        print(f"  {r['mutant']:8} measured={r['measured']:+.5f} "
              f"pred={r['predicted']:+.5f} err={r['error']:+.5f} "
              f"{'IN RANGE' if r['in_range'] else ''}")
    print("worst 5 predictions:")
    for r in rows[-5:]:
        print(f"  {r['mutant']:8} measured={r['measured']:+.5f} "
              f"pred={r['predicted']:+.5f} err={r['error']:+.5f}")

    t315 = [r for r in rows if r["mutant"].upper().startswith("T315")]
    if t315:
        print("\nT315 substitutions in the validation set:")
        for r in t315:
            print(f"  {r['mutant']:8} measured={r['measured']:+.5f} "
                  f"pred={r['predicted']:+.5f} err={r['error']:+.5f} "
                  f"{'IN RANGE' if r['in_range'] else ''}")
    else:
        print("\nT315 is in the TRAINING half for this seed "
              "(position-grouped split; use --seed to move it)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
