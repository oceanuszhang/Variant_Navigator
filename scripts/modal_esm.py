#!/usr/bin/env python3
"""ESM-2 embeddings on Modal GPU — the 650M model, which is impractical locally.

    /opt/anaconda3/bin/modal run scripts/modal_esm.py                  # 650M
    /opt/anaconda3/bin/modal run scripts/modal_esm.py --model esm2_t12_35M_UR50D

WHY MODAL AND NOT THE LAPTOP. The 35M model runs fine on MPS and is what the
README's numbers use. The 650M model is a 2.5 GB download and far slower per
pass; on conference wifi it stalled twice locally. On Modal the weights land in
a persistent Volume once, and the forward passes run on an A10G, so the larger
model becomes a sensible default rather than an overnight job.

WHAT IT COMPUTES. The same mutant-marginal feature as
`scripts/embed_variants.py`, so results are directly comparable:

    feature = h_mut[i] - h_ref[i]     at the mutated residue i

Numbering is validated against the reference sequence before any GPU time is
spent — ABL1a == UniProt P00519 at offset 0, verified at six residues. ESM-2
holds 1024 position embeddings, so the sequence is cropped to 1-1022, which
contains the kinase domain (242-512).
"""

from __future__ import annotations

import modal

APP_NAME = "variantnavigator-esm"
SMALL = "esm2_t12_35M_UR50D"
LARGE = "esm2_t33_650M_UR50D"
MAX_LEN = 1022

image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install("torch==2.5.1", "fair-esm==2.0.0", "numpy==2.1.3")
    .env({"TORCH_HOME": "/weights"})
)
weights = modal.Volume.from_name("vn-esm-weights", create_if_missing=True)
app = modal.App(APP_NAME)


@app.function(image=image, gpu="A10G", volumes={"/weights": weights},
              timeout=3600, retries=1)
def embed_mutants(
    crop: str, jobs: list[tuple[int, str]], lo: int, model_name: str, batch: int = 32
) -> list[list[float]]:
    """Mutant-marginal difference vectors for (position, alt_aa) pairs.

    `jobs` positions are in the ORIGINAL sequence numbering; `lo` is the crop
    start so the index into `crop` is position - lo.
    """
    import esm as esm_lib
    import torch

    model, alphabet = getattr(esm_lib.pretrained, model_name)()
    model = model.cuda().eval()
    weights.commit()  # persist the downloaded checkpoint for the next run
    bc = alphabet.get_batch_converter()
    layer = model.num_layers

    def run(seqs):
        _l, _s, toks = bc(seqs)
        with torch.no_grad():
            out = model(toks.cuda(), repr_layers=[layer])
        return out["representations"][layer].float().cpu().numpy()

    ref = run([("ref", crop)])[0]
    out: list[list[float]] = []
    for start in range(0, len(jobs), batch):
        chunk = jobs[start : start + batch]
        seqs = []
        for pos, alt in chunk:
            i = pos - lo
            seqs.append((f"p{pos}{alt}", crop[:i] + alt + crop[i + 1 :]))
        reps = run(seqs)
        for row, (pos, _alt) in enumerate(chunk):
            t = (pos - lo) + 1
            out.append((reps[row, t] - ref[t]).tolist())
        print(f"  {min(start + batch, len(jobs))}/{len(jobs)}", flush=True)
    return out


@app.local_entrypoint()
def main(model: str = LARGE, limit: int = 0, batch: int = 32):
    import json
    import sys
    import urllib.request
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    import numpy as np

    from pipeline import dms

    out_dir = Path("data/esm")
    out_dir.mkdir(parents=True, exist_ok=True)
    seq_cache = out_dir / "P00519.seq"
    if seq_cache.exists():
        seq = seq_cache.read_text().strip()
    else:
        url = "https://rest.uniprot.org/uniprotkb/P00519.json?fields=sequence"
        seq = json.load(urllib.request.urlopen(url, timeout=60))["sequence"]["value"]
        seq_cache.write_text(seq)
    print(f"P00519: {len(seq)} aa")

    variants = dms.load_netgrowth()
    if limit:
        variants = variants[:limit]
    lo, hi = 1, min(len(seq), MAX_LEN)
    crop = seq[lo - 1 : hi]

    bad = [v for v in variants
           if not (lo <= v.position <= hi) or seq[v.position - 1].upper() != v.ref_aa.upper()]
    if bad:
        raise SystemExit(f"{len(bad)} variants disagree with the reference; refusing to spend GPU")
    print(f"{len(variants)} variants validated; crop {lo}-{hi}")

    jobs = [(v.position, v.alt_aa.upper()) for v in variants]
    print(f"submitting {len(jobs)} mutants to Modal ({model}, A10G)…")
    vecs = embed_mutants.remote(crop, jobs, lo, model, batch)

    X = np.asarray(vecs, dtype=np.float32)
    y = np.array([v.response for v in variants], dtype=np.float32)
    groups = np.array([v.position for v in variants], dtype=np.int32)
    rep1 = np.array([np.nan if v.rep1 is None else v.rep1 for v in variants], dtype=np.float32)
    rep2 = np.array([np.nan if v.rep2 is None else v.rep2 for v in variants], dtype=np.float32)
    names = np.array([v.mutant for v in variants])

    stem = out_dir / f"P00519_{model}_mutdiff"
    np.savez_compressed(f"{stem}.npz", X=X, y=y, groups=groups, rep1=rep1,
                        rep2=rep2, names=names, crop=np.array([lo, hi]),
                        model=np.array([model]))
    print(f"wrote {stem}.npz  X={X.shape}")
    print(f"mean |feature| = {float(np.abs(X).mean()):.5f}")
