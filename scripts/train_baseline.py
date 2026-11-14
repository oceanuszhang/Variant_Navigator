#!/usr/bin/env python3
"""Train and honestly evaluate net-growth prediction on the BCR-ABL DMS screen.

    python scripts/train_baseline.py

Reports an ablation ladder under POSITION-GROUPED cross-validation, against the
replicate noise ceiling. The ladder is the point: it shows which feature block
is doing the work, and whether a model beats the trivial predictor at all.

    global mean        predict the training mean for everything
    position only      distance to the ligand + where in the domain
    substitution only  BLOSUM62 + hydropathy / volume / charge change
    all                both blocks

Held out by RESIDUE POSITION, never by variant: all ~19 substitutions at a
position share their structural features, so a random split puts near-duplicate
rows on both sides of the fold and inflates the score. `--leaky` reproduces that
mistake on purpose so the size of the inflation is visible.
"""

from __future__ import annotations

import argparse
import statistics as st
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pipeline import dms, features  # noqa: E402

try:
    from sklearn.ensemble import HistGradientBoostingRegressor, RandomForestRegressor
    from sklearn.linear_model import RidgeCV
    from sklearn.preprocessing import StandardScaler
except ImportError:  # pragma: no cover
    print("needs scikit-learn (already present on this machine): pip install scikit-learn")
    raise SystemExit(1)


def pearson(xs: list[float], ys: list[float]) -> float:
    n = len(xs)
    if n < 3:
        return float("nan")
    mx, my = st.mean(xs), st.mean(ys)
    sx, sy = st.pstdev(xs), st.pstdev(ys)
    if sx == 0 or sy == 0:
        return 0.0
    return sum((a - mx) * (b - my) for a, b in zip(xs, ys)) / n / (sx * sy)


def spearman(xs: list[float], ys: list[float]) -> float:
    def rank(v: list[float]) -> list[float]:
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    return pearson(rank(xs), rank(ys))


def mae(a: list[float], b: list[float]) -> float:
    return st.mean(abs(x - y) for x, y in zip(a, b))


def make_models() -> dict:
    return {
        "ridge": lambda: RidgeCV(alphas=[0.01, 0.1, 1.0, 10.0, 100.0]),
        "hist_gbm": lambda: HistGradientBoostingRegressor(
            max_iter=400, learning_rate=0.06, max_depth=6,
            min_samples_leaf=20, l2_regularization=1.0, random_state=0,
        ),
        "random_forest": lambda: RandomForestRegressor(
            n_estimators=300, min_samples_leaf=2, n_jobs=-1, random_state=0,
        ),
    }


def evaluate(
    variants: list, *, blocks: str, model_name: str, folds: int, leaky: bool,
    esm_marginals: dict | None = None,
) -> dict:
    """Out-of-fold predictions, pooled, then scored once."""
    X, y, groups, names = features.featurize(
        variants, blocks=blocks, esm_marginals=esm_marginals
    )
    builder = make_models()[model_name] if model_name != "mean" else None

    oof_true: list[float] = []
    oof_pred: list[float] = []

    if leaky:
        # Deliberately wrong: split by row index, ignoring position.
        idx_folds = [[i for i in range(len(y)) if i % folds == f] for f in range(folds)]
    else:
        positions = sorted({g for g in groups})
        idx_folds = [
            [i for i, g in enumerate(groups)
             if positions.index(g) % folds == f]
            for f in range(folds)
        ]

    for f in range(folds):
        test_idx = set(idx_folds[f])
        tr = [i for i in range(len(y)) if i not in test_idx]
        te = sorted(test_idx)
        if not te or not tr:
            continue
        ytr = [y[i] for i in tr]
        if model_name == "mean":
            pred = [st.mean(ytr)] * len(te)
        else:
            Xtr = [X[i] for i in tr]
            Xte = [X[i] for i in te]
            sc = StandardScaler().fit(Xtr)
            m = builder()
            m.fit(sc.transform(Xtr), ytr)
            pred = list(m.predict(sc.transform(Xte)))
        oof_true += [y[i] for i in te]
        oof_pred += pred

    return {
        "blocks": blocks,
        "model": model_name,
        "n": len(oof_true),
        "pearson_r": pearson(oof_pred, oof_true),
        "spearman": spearman(oof_pred, oof_true),
        "mae": mae(oof_pred, oof_true),
        "n_features": len(names),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--leaky", action="store_true",
                    help="split by row instead of position, to show the inflation")
    ap.add_argument("--esm", metavar="MODEL", nargs="?",
                    const="esm2_t33_650M_UR50D", default=None,
                    help="add ESM2 masked-marginal features (needs torch+fair-esm)")
    args = ap.parse_args()

    vs = dms.load_netgrowth()
    nf = dms.noise_floor(vs)
    mar = None
    if args.esm:
        import json, urllib.request
        from pipeline import esm2_score
        if not esm2_score.available():
            print("ESM requested but torch/fair-esm missing; run with the venv python")
            return 1
        seq = json.load(urllib.request.urlopen(
            "https://rest.uniprot.org/uniprotkb/P00519.json?fields=sequence",
            timeout=60))["sequence"]["value"]
        mar = esm2_score.masked_marginals(
            seq, sorted({v.position for v in vs}),
            model_name=args.esm, accession="P00519",
        )
    print(f"BCR-ABL DMS — {len(vs)} variants, {len({v.position for v in vs})} positions")
    print(f"replicate r = {nf['replicate_pearson_r']}  ->  attenuated ceiling "
          f"r = {nf['attenuated_r_ceiling']}   (scoring above this = fitting noise)")
    print(f"replicate MAE = {nf['replicate_mae']}   response sd = {nf['response_sd']}")
    split = "ROW-WISE (LEAKY, for contrast)" if args.leaky else "POSITION-GROUPED"
    print(f"\n{args.folds}-fold {split} cross-validation, pooled out-of-fold\n")

    print(f"{'features':16} {'model':14} {'r':>8} {'spearman':>9} {'MAE':>9} {'feats':>6}")
    print("-" * 68)
    rows = []
    plan = [
        ("all", "mean"),
        ("position", "ridge"),
        ("position", "hist_gbm"),
        ("substitution", "ridge"),
        ("substitution", "hist_gbm"),
        ("all", "ridge"),
        ("all", "hist_gbm"),
        ("all", "random_forest"),
    ]
    if mar is not None:
        plan += [
            ("esm", "ridge"),
            ("esm", "hist_gbm"),
            ("substitution+esm", "ridge"),
            ("all+esm", "ridge"),
            ("all+esm", "hist_gbm"),
            ("all+esm", "random_forest"),
        ]
    for blocks, model in plan:
        r = evaluate(vs, blocks=blocks, model_name=model,
                     folds=args.folds, leaky=args.leaky, esm_marginals=mar)
        rows.append(r)
        label = "global mean" if model == "mean" else blocks
        print(f"{label:16} {model:14} {r['pearson_r']:>8.4f} "
              f"{r['spearman']:>9.4f} {r['mae']:>9.5f} {r['n_features']:>6}")

    best = max((r for r in rows if r["model"] != "mean"), key=lambda r: r["pearson_r"])
    ceil = nf["attenuated_r_ceiling"]
    print("-" * 68)
    print(f"best: {best['blocks']} / {best['model']}  r={best['pearson_r']:.4f}")
    print(f"  that is {best['pearson_r'] / ceil * 100:.0f}% of the attainable ceiling "
          f"({ceil}), not of 1.0")
    print(f"  MAE {best['mae']:.5f} vs replicate MAE {nf['replicate_mae']} "
          f"({best['mae'] / nf['replicate_mae']:.1f}x the noise floor)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
