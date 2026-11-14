#!/usr/bin/env python3
"""Mutant-marginal ESM-2 embeddings: one vector per measured variant.

    .venv/bin/python scripts/embed_variants.py

For each of the 4,925 measured BCR-ABL variants, build the mutant sequence, run
it through ESM-2, and take the hidden state AT THE MUTATED RESIDUE. The feature
is the difference from the reference run at the same position:

    feature = h_mut[i] - h_ref[i]        (480-dim for the 35M model)

WHY THE DIFFERENCE, AND WHY AT ONE RESIDUE. Mean-pooling the whole sequence
would be useless here: every variant is the same protein, so pooled vectors are
~identical (99.82% sequence identity). Reading one residue keeps the signal, and
subtracting the reference removes everything about the position that is constant
across its ~19 substitutions, leaving what the substitution itself did.

NUMBERING. The DMS uses ABL1a numbering, which agrees with UniProt P00519
(1130 aa) at offset 0 - verified at six residues, not assumed. ESM-2 has 1024
position embeddings, so the sequence is cropped to residues 1-1022; the kinase
domain (242-512) sits well inside that. Every mutant is checked against the
reference residue before being embedded, and a mismatch raises.

Cost: 4,926 forward passes (one reference + one per variant). Minutes on MPS
with the 35M model. Cached to data/esm/, so this runs once.
"""

from __future__ import annotations

import argparse
import json
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pipeline import dms  # noqa: E402
from pipeline.esm2_score import MAX_LEN, SMALL_MODEL, _device  # noqa: E402

OUT_DIR = Path("data/esm")
UNIPROT = "P00519"


def reference_sequence() -> str:
    cache = OUT_DIR / f"{UNIPROT}.seq"
    if cache.exists():
        return cache.read_text().strip()
    url = f"https://rest.uniprot.org/uniprotkb/{UNIPROT}.json?fields=sequence"
    seq = json.load(urllib.request.urlopen(url, timeout=60))["sequence"]["value"]
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    cache.write_text(seq)
    return seq


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--model", default=SMALL_MODEL)
    ap.add_argument("--batch", type=int, default=16)
    ap.add_argument("--limit", type=int, default=0, help="debug: only N variants")
    args = ap.parse_args()

    import numpy as np
    import torch
    import esm as esm_lib

    seq = reference_sequence()
    print(f"{UNIPROT}: {len(seq)} aa")
    if len(seq) != 1130:
        print(f"  WARNING: expected 1130 aa for P00519, got {len(seq)}")

    variants = dms.load_netgrowth()
    if args.limit:
        variants = variants[: args.limit]
    lo, hi = 1, min(len(seq), MAX_LEN)
    crop = seq[lo - 1 : hi]
    print(f"crop {lo}-{hi} ({len(crop)} aa); kinase domain 242-512 is inside")

    # Validate every variant against the reference BEFORE spending any compute.
    bad = [
        v for v in variants
        if not (lo <= v.position <= hi) or seq[v.position - 1].upper() != v.ref_aa.upper()
    ]
    if bad:
        for v in bad[:5]:
            got = seq[v.position - 1] if v.position <= len(seq) else "?"
            print(f"  MISMATCH {v.mutant}: reference has {got} at {v.position}")
        raise SystemExit(f"{len(bad)} variants disagree with the reference - refusing to embed")
    print(f"{len(variants)} variants, all validated against the reference")

    dev = _device()
    print(f"loading {args.model} on {dev}")
    model, alphabet = getattr(esm_lib.pretrained, args.model)()
    model = model.to(dev).eval()
    bc = alphabet.get_batch_converter()
    repr_layer = model.num_layers

    def embed(seqs: list[tuple[str, str]]) -> "np.ndarray":
        _l, _s, toks = bc(seqs)
        with torch.no_grad():
            out = model(toks.to(dev), repr_layers=[repr_layer])
        return out["representations"][repr_layer].float().cpu().numpy()

    # Reference pass: hidden state at every position, once.
    ref = embed([("ref", crop)])[0]          # (len+2, dim)
    dim = ref.shape[-1]
    print(f"reference embedded, hidden dim = {dim}")

    X = np.zeros((len(variants), dim), dtype=np.float32)
    y = np.zeros(len(variants), dtype=np.float32)
    groups = np.zeros(len(variants), dtype=np.int32)
    rep1 = np.full(len(variants), np.nan, dtype=np.float32)
    rep2 = np.full(len(variants), np.nan, dtype=np.float32)
    names: list[str] = []

    for start in range(0, len(variants), args.batch):
        chunk = variants[start : start + args.batch]
        seqs = []
        for v in chunk:
            idx = v.position - lo
            seqs.append((v.mutant, crop[:idx] + v.alt_aa.upper() + crop[idx + 1 :]))
        reps = embed(seqs)
        for row, v in enumerate(chunk):
            tpos = (v.position - lo) + 1          # +1 for BOS
            i = start + row
            X[i] = reps[row, tpos] - ref[tpos]    # mutant-marginal difference
            y[i] = v.response
            groups[i] = v.position
            if v.rep1 is not None:
                rep1[i] = v.rep1
            if v.rep2 is not None:
                rep2[i] = v.rep2
            names.append(v.mutant)
        done = min(start + args.batch, len(variants))
        print(f"  {done}/{len(variants)} embedded", end="\r", flush=True)
    print()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    # A --limit run must never overwrite a full one: tag partial outputs.
    suffix = "_mutdiff" if not args.limit else f"_mutdiff_n{len(variants)}"
    stem = OUT_DIR / f"{UNIPROT}_{args.model}{suffix}"
    np.savez_compressed(
        f"{stem}.npz", X=X, y=y, groups=groups, rep1=rep1, rep2=rep2,
        names=np.array(names), crop=np.array([lo, hi]), model=np.array([args.model]),
    )
    print(f"wrote {stem}.npz  X={X.shape}  y={y.shape}")
    nz = float(np.abs(X).mean())
    print(f"mean |feature| = {nz:.5f} (0 would mean the difference carries nothing)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
