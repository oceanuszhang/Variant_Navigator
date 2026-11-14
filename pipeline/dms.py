"""BCR-ABL deep mutational scanning: the public dataset that unblocks training.

WHY THIS REPLACED THE TRANSCRIBED CSVs. The original plan needed two files
hand-transcribed from published figures (`egfr_ic50.csv`, `rice_gr50.csv`), which
no API could supply and which blocked everything downstream of feature
extraction. This dataset removes that blocker entirely: it is public, MIT
licensed, programmatically fetchable, and ~1000x larger.

    Inam H.*, Tomaszkiewicz M.*, Reynolds J.A., Yang Z., Leighow S.M.,
    Pritchard J.R. "Characterizing Variants of Uncertain Drug Response (VUDRs)
    Using Quantitative Measurements at Clinical Exposures." (*equal
    contribution). github.com/pritchardlabatpsu/abl_dms  (MIT)

4,925 single-amino-acid variants of the BCR-ABL kinase domain (residues
242-512), each with a net growth rate in Ba/F3 cells at 1200 nM imatinib, in two
replicates. Ba/F3 is the same assay system as the repo's Experiment 1, so the
schema carries over.

NUMBERING, verified rather than assumed. The DMS uses ABL1a numbering, where the
gatekeeper is Thr315. That agrees with UniProt P00519 (1130 aa) at offset 0, and
with PDB 1IEP and 2HYY at offset 0. It does NOT agree with 1OPJ, which needs
+19. `pipeline.pocket.validate_references()` is the guard.

WHAT THE RESPONSE ACTUALLY MEASURES, and why this matters for modelling. Net
growth under drug is not a resistance score. Ba/F3 cells here depend on BCR-ABL
to survive, so the readout confounds two mechanisms:

    kinase still functional?   a loss-of-function variant grows slowly whatever
                               the drug does
    imatinib still bound?      a variant that blocks imatinib grows fast

Those pull in opposite directions, and a single feature cannot separate them.
This is the central modelling fact about this dataset: see README "What the
measured data taught us".
"""

from __future__ import annotations

import csv
import os
import statistics as st
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterator

RAW_BASE = (
    "https://raw.githubusercontent.com/pritchardlabatpsu/abl_dms/master/Figures"
)
# 1200 nM imatinib net-growth table (manuscript Fig. 2D).
NETGROWTH_URL = f"{RAW_BASE}/Figure_2D/high_imatinib_netgrowth_data.csv"
# Per-variant IC50 estimates fitted across 300/600/1200 nM (manuscript Fig. 3D).
IC50_URL = f"{RAW_BASE}/Figure_3D/all_mutants_with_IC50_estimates.csv"

CACHE_DIR = Path("data/datasets")
NETGROWTH_CSV = CACHE_DIR / "abl_dms_netgrowth_1200nM.csv"
IC50_CSV = CACHE_DIR / "abl_dms_ic50_estimates.csv"

UNIPROT = "P00519"          # ABL1, 1130 aa, agrees with DMS numbering at offset 0
DOMAIN = (242, 512)         # kinase domain tiled by the screen
DRUG = "imatinib"
DOSE_NM = 1200
TIMEOUT = 120

# The response column. `netgrowth_mean_corrected` is the replicate mean after
# the authors' correction; lower = more susceptible at this dose.
RESPONSE = "netgrowth_mean_corrected"
REPLICATES = ("netgr_corr_rep1.high", "netgr_corr_rep2.high")


class DMSError(RuntimeError):
    """Dataset missing, unreadable, or failing a sanity check."""


@dataclass(frozen=True)
class Variant:
    """One measured single-amino-acid variant."""

    mutant: str          # e.g. "T315I", in ABL1a numbering
    position: int        # ABL1a == UniProt P00519 numbering
    ref_aa: str
    alt_aa: str
    response: float      # net growth at 1200 nM imatinib; higher = more escape
    rep1: float | None
    rep2: float | None

    @property
    def replicate_spread(self) -> float | None:
        """|rep1 - rep2|: this variant's own measurement uncertainty."""
        if self.rep1 is None or self.rep2 is None:
            return None
        return abs(self.rep1 - self.rep2)


def _download(url: str, dest: Path) -> Path:
    """Fetch once and cache. Offline runs reuse the cached copy."""
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        with urllib.request.urlopen(url, timeout=TIMEOUT) as resp:
            data = resp.read()
    except Exception as exc:  # noqa: BLE001
        raise DMSError(f"could not fetch {url}: {exc}") from exc
    if len(data) < 1024:
        raise DMSError(f"{url} returned only {len(data)} bytes")
    dest.write_bytes(data)
    return dest


def _as_float(raw: str | None) -> float | None:
    if raw is None:
        return None
    raw = raw.strip()
    if not raw or raw.upper() in {"NA", "NAN", "NULL", "INF", "-INF"}:
        return None
    try:
        v = float(raw)
    except ValueError:
        return None
    return None if v != v else v  # drop NaN


def load_netgrowth(*, path: Path | None = None) -> list[Variant]:
    """Load the 1200 nM net-growth table, fetching and caching if needed.

    Rows without a usable response are dropped rather than imputed: a missing
    measurement is not a zero, and imputing one would quietly invent data.
    """
    src = path or _download(NETGROWTH_URL, NETGROWTH_CSV)
    out: list[Variant] = []
    with src.open(newline="") as fh:
        for row in csv.DictReader(fh):
            y = _as_float(row.get(RESPONSE))
            pos = _as_float(row.get("protein_start"))
            if y is None or pos is None:
                continue
            ref, alt = (row.get("ref_aa") or "").strip(), (row.get("alt_aa") or "").strip()
            if len(ref) != 1 or len(alt) != 1 or ref == alt:
                continue  # synonymous or malformed
            out.append(Variant(
                mutant=(row.get("mutant") or f"{ref}{int(pos)}{alt}").strip(),
                position=int(pos),
                ref_aa=ref,
                alt_aa=alt,
                response=y,
                rep1=_as_float(row.get(REPLICATES[0])),
                rep2=_as_float(row.get(REPLICATES[1])),
            ))
    if not out:
        raise DMSError(f"no usable rows parsed from {src}")
    return out


def noise_floor(variants: list[Variant]) -> dict[str, Any]:
    """Replicate agreement — the error floor any model should be judged against.

    A model whose test error is below this is fitting measurement noise, not
    biology, so this is reported alongside every score. The attenuated ceiling
    is the correlation a *perfect* predictor of the true signal could reach
    given noisy labels: sqrt of the replicate correlation.
    """
    pairs = [(v.rep1, v.rep2) for v in variants if v.rep1 is not None and v.rep2 is not None]
    if len(pairs) < 3:
        raise DMSError("too few replicated rows to estimate a noise floor")
    a = [x for x, _ in pairs]
    b = [y for _, y in pairs]
    r = _pearson(a, b)
    return {
        "n_replicated": len(pairs),
        "replicate_pearson_r": round(r, 4),
        "replicate_mae": round(st.mean(abs(x - y) for x, y in pairs), 5),
        # Spearman-Brown style: a perfect model cannot beat this against the
        # mean of two noisy replicates.
        "attenuated_r_ceiling": round(min(1.0, (2 * r / (1 + r)) ** 0.5), 4),
        "response_sd": round(st.pstdev([v.response for v in variants]), 5),
    }


def _pearson(xs: list[float], ys: list[float]) -> float:
    n = len(xs)
    mx, my = st.mean(xs), st.mean(ys)
    sx, sy = st.pstdev(xs), st.pstdev(ys)
    if sx == 0 or sy == 0:
        return 0.0
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / n / (sx * sy)


def split_by_position(
    variants: list[Variant], *, fold: int = 0, folds: int = 5
) -> tuple[list[Variant], list[Variant]]:
    """Group-split held out by RESIDUE POSITION, not by variant.

    This is load-bearing. Structural features (distance to the ligand, burial,
    conservation) are shared by all ~19 substitutions at a position, so a random
    variant-level split puts near-duplicate rows on both sides and inflates the
    score. Holding out whole positions measures what we actually want: does this
    generalise to a residue the model has never seen?

    Deterministic, so a reported number can be reproduced.
    """
    if not 0 <= fold < folds:
        raise ValueError(f"fold {fold} out of range for {folds} folds")
    positions = sorted({v.position for v in variants})
    test_pos = {p for i, p in enumerate(positions) if i % folds == fold}
    train = [v for v in variants if v.position not in test_pos]
    test = [v for v in variants if v.position in test_pos]
    return train, test


def iter_folds(
    variants: list[Variant], *, folds: int = 5
) -> Iterator[tuple[int, list[Variant], list[Variant]]]:
    for f in range(folds):
        tr, te = split_by_position(variants, fold=f, folds=folds)
        yield f, tr, te


def provenance() -> dict[str, Any]:
    """Everything needed to cite and re-fetch this dataset."""
    return {
        "dataset": "BCR-ABL kinase-domain DMS, imatinib",
        "source_repo": "https://github.com/pritchardlabatpsu/abl_dms",
        "license": "MIT (Pritchard Lab, Penn State)",
        "citation": (
            "Inam H.*, Tomaszkiewicz M.*, Reynolds J.A., Yang Z., Leighow S.M., "
            "Pritchard J.R. Characterizing Variants of Uncertain Drug Response "
            "(VUDRs) Using Quantitative Measurements at Clinical Exposures."
        ),
        "netgrowth_url": NETGROWTH_URL,
        "ic50_url": IC50_URL,
        "assay": "Ba/F3 net growth rate, 6 days, TileSeq + duplex sequencing",
        "dose_nM": DOSE_NM,
        "drug": DRUG,
        "uniprot": UNIPROT,
        "numbering": "ABL1a; == UniProt P00519 and PDB 1IEP/2HYY at offset 0",
        "domain_residues": list(DOMAIN),
        "response_column": RESPONSE,
        "cached": {
            "netgrowth": str(NETGROWTH_CSV) if NETGROWTH_CSV.exists() else None,
            "ic50": str(IC50_CSV) if IC50_CSV.exists() else None,
        },
    }


def available() -> bool:
    """True when the data is cached, or fetchable if the network is up."""
    if NETGROWTH_CSV.exists() and NETGROWTH_CSV.stat().st_size > 0:
        return True
    return not os.getenv("VN_OFFLINE")


def _main() -> int:
    print("BCR-ABL DMS dataset")
    print("=" * 68)
    vs = load_netgrowth()
    print(f"variants           {len(vs)}")
    positions = sorted({v.position for v in vs})
    print(f"positions          {len(positions)}  ({positions[0]}-{positions[-1]})")
    nf = noise_floor(vs)
    print(f"\nnoise floor (the number to beat, and not to beat by too much)")
    for k, v in nf.items():
        print(f"  {k:22} {v}")
    tr, te = split_by_position(vs, fold=0)
    print(f"\nposition-grouped split, fold 0: train={len(tr)} test={len(te)}")
    print(f"  leakage check: shared positions = "
          f"{len({v.position for v in tr} & {v.position for v in te})}")
    hi = sorted(vs, key=lambda v: -v.response)[:5]
    print("\nmost imatinib-escaping variants measured:")
    for v in hi:
        print(f"  {v.mutant:8} netgrowth={v.response:+.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
