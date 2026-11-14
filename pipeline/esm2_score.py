"""ESM2 variant scoring by masked marginals.

WHY MASKED MARGINALS AND NOT SEQUENCE EMBEDDINGS. Every variant in the BCR-ABL
screen is the same protein; any two mutant sequences are 99.82% identical. A
mean-pooled whole-sequence embedding therefore returns an almost identical vector
for all 4,925 variants and a model fit on it has nothing to separate them. The
standard fix, from ESM-1v, is to mask the varied position and read the model's
own distribution over residues there:

    score(ref -> alt at i) = log P(alt | x masked at i) - log P(ref | x masked at i)

One forward pass per *position*, not per variant, so 270 passes cover all 4,925
variants. The score varies per substitution, which is what the geometry features
could not do.

LENGTH LIMIT, and it bites here. ESM2 is trained with 1024 position embeddings
while ABL1 P00519 is 1130 aa, so the full sequence cannot be fed in. We crop to
residues 1-1022, which still contains the whole kinase domain (242-512) with
~240 residues of leading context and the crop is recorded in the cache so a
score can never be silently compared across different crops.

WHAT THIS DOES AND DOES NOT MEASURE. ESM scores evolutionary likelihood, not
drug resistance. Expect it to capture whether a substitution breaks the kinase -
which matters, because Ba/F3 cells here die without functional BCR-ABL - and to
say nothing about whether imatinib still binds. Thr315 is both highly conserved
and the most famous imatinib-escape residue, so on the resistance axis alone ESM
should look wrong. It is a complement to pocket geometry, not a replacement.

Needs `torch` and `fair-esm`. Both are optional: `available()` is False without
them and every caller skips rather than crashing.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

CACHE_DIR = Path("data/esm")
SMALL_MODEL = "esm2_t12_35M_UR50D"   # ~135 MB, benchmarked, the default
LARGE_MODEL = "esm2_t33_650M_UR50D"  # ~2.5 GB download, better but slow to get

# Default to the small model on purpose. The 650M weights are a 2.5 GB download
# that will stall a demo on conference wifi, and the numbers in the README were
# measured with the 35M model. Override deliberately:
#   VN_ESM_MODEL=esm2_t33_650M_UR50D
DEFAULT_MODEL = os.getenv("VN_ESM_MODEL") or SMALL_MODEL
MAX_LEN = 1022  # 1024 minus the BOS/EOS tokens


@dataclass
class VariantScore:
    variant: str
    delta_ll: float          # log P(alt) - log P(ref), masked at the position
    ll_ref: float
    ll_alt: float
    model: str
    position: int

    @property
    def __dict__(self) -> dict[str, Any]:  # type: ignore[override]
        return {
            "variant": self.variant,
            "delta_ll": self.delta_ll,
            "ll_ref": self.ll_ref,
            "ll_alt": self.ll_alt,
            "model": self.model,
            "position": self.position,
        }


class ESMUnavailable(RuntimeError):
    """torch / fair-esm missing, or no cached scores to fall back on."""


def available() -> bool:
    try:
        import esm  # noqa: F401
        import torch  # noqa: F401
    except ImportError:
        return False
    return True


def _device(prefer: str | None = None) -> str:
    import torch
    if prefer:
        return prefer
    if torch.backends.mps.is_available():
        return "mps"
    if torch.cuda.is_available():
        return "cuda"
    return "cpu"


def _cache_path(model: str, accession: str, lo: int, hi: int) -> Path:
    return CACHE_DIR / f"{accession}_{model}_{lo}-{hi}_masked_marginals.json"


def masked_marginals(
    sequence: str,
    positions: list[int],
    *,
    model_name: str = DEFAULT_MODEL,
    accession: str = "seq",
    device: str | None = None,
    batch_positions: int = 8,
    verbose: bool = True,
) -> dict[str, Any]:
    """log-softmax over the 20 residues at each requested position.

    Returns {"log_probs": {position: {aa: logp}}, ...} in 1-indexed numbering of
    the ORIGINAL sequence. Cached to disk, keyed by model and crop.
    """
    if not available():
        raise ESMUnavailable("torch / fair-esm not installed")
    import torch
    import esm as esm_lib

    lo, hi = 1, min(len(sequence), MAX_LEN)
    needed = [p for p in positions if lo <= p <= hi]
    dropped = sorted(set(positions) - set(needed))
    cache = _cache_path(model_name, accession, lo, hi)
    if cache.exists():
        data = json.loads(cache.read_text())
        if set(str(p) for p in needed) <= set(data.get("log_probs", {})):
            if verbose:
                print(f"  esm: reusing {cache}")
            return data

    crop = sequence[lo - 1:hi]
    dev = _device(device)
    if verbose:
        print(f"  esm: loading {model_name} on {dev} "
              f"(crop {lo}-{hi} of {len(sequence)} aa, {len(needed)} positions)")

    loader = getattr(esm_lib.pretrained, model_name)
    model, alphabet = loader()
    model = model.to(dev).eval()
    bc = alphabet.get_batch_converter()
    aas = "ACDEFGHIKLMNPQRSTVWY"
    tok = {a: alphabet.get_idx(a) for a in aas}

    out: dict[str, dict[str, float]] = {}
    with torch.no_grad():
        for start in range(0, len(needed), batch_positions):
            chunk = needed[start:start + batch_positions]
            seqs = []
            for p in chunk:
                idx = p - lo                      # 0-based into the crop
                masked = crop[:idx] + "<mask>" + crop[idx + 1:]
                seqs.append((f"p{p}", masked))
            _lbl, _strs, toks = bc(seqs)
            toks = toks.to(dev)
            logits = model(toks)["logits"]
            lp = torch.log_softmax(logits, dim=-1)
            for row, p in enumerate(chunk):
                # +1 for the prepended BOS token
                tpos = (p - lo) + 1
                out[str(p)] = {a: float(lp[row, tpos, tok[a]]) for a in aas}
            if verbose:
                done = min(start + batch_positions, len(needed))
                print(f"  esm: {done}/{len(needed)} positions", end="\r", flush=True)
    if verbose:
        print()

    data = {
        "model": model_name,
        "accession": accession,
        "crop": [lo, hi],
        "sequence_length": len(sequence),
        "dropped_positions": dropped,
        "log_probs": out,
    }
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache.write_text(json.dumps(data))
    if verbose:
        print(f"  esm: wrote {cache}")
    return data


def delta_ll(
    marginals: dict[str, Any], position: int, ref_aa: str, alt_aa: str
) -> float | None:
    """log P(alt) - log P(ref) at a masked position. None if not scored."""
    row = (marginals.get("log_probs") or {}).get(str(position))
    if not row:
        return None
    r, a = row.get(ref_aa.upper()), row.get(alt_aa.upper())
    if r is None or a is None:
        return None
    return a - r


def score_variant(
    wt_sequence: str,
    mt_sequence: str,
    variant: str,
    *,
    model: str | None = None,
    accession: str = "seq",
) -> VariantScore:
    """Score one substitution given as e.g. "T315I".

    `model=None` resolves DEFAULT_MODEL at call time, so overriding the module
    attribute works. Binding it as a default argument would freeze it at import
    and silently trigger the 2.5 GB download the override was meant to avoid.
    """
    import re
    m = re.match(r"^([A-Z])(\d+)([A-Z])$", variant.strip().upper())
    if not m:
        raise ESMUnavailable(f"score_variant handles single substitutions, got {variant!r}")
    ref, pos, alt = m.group(1), int(m.group(2)), m.group(3)
    if pos <= len(wt_sequence) and wt_sequence[pos - 1].upper() != ref:
        raise ESMUnavailable(
            f"{variant}: sequence has {wt_sequence[pos-1]} at {pos}, not {ref}"
        )
    mar = masked_marginals(
        wt_sequence, [pos], model_name=model or DEFAULT_MODEL,
        accession=accession, verbose=False,
    )
    row = (mar.get("log_probs") or {}).get(str(pos))
    if not row:
        raise ESMUnavailable(f"{variant}: position {pos} outside the scored crop")
    return VariantScore(
        variant=variant.upper(),
        delta_ll=row[alt] - row[ref],
        ll_ref=row[ref],
        ll_alt=row[alt],
        model=model or DEFAULT_MODEL,
        position=pos,
    )


def _main() -> int:
    import urllib.request
    acc = "P00519"
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.json?fields=sequence"
    seq = json.load(urllib.request.urlopen(url, timeout=60))["sequence"]["value"]
    model = os.getenv("VN_ESM_MODEL", DEFAULT_MODEL)
    print(f"ESM2 masked marginals — {acc} ({len(seq)} aa), model {model}")
    mar = masked_marginals(
        seq, list(range(242, 513)), model_name=model, accession=acc
    )
    print(f"\nscored positions: {len(mar['log_probs'])}  crop {mar['crop']}")
    print("\nsanity: known imatinib-resistance substitutions")
    for v in ("T315I", "Y253H", "E255K", "M351T", "H396R", "G250E"):
        d = delta_ll(mar, int(v[1:-1]), v[0], v[-1])
        print(f"  {v:7} delta_ll = {d:+.3f}" if d is not None else f"  {v}: n/a")
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
