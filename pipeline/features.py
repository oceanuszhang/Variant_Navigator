"""Features for predicting variant effect on drug response.

Two kinds, and the split matters more than any individual feature:

    POSITION-LEVEL   distance to the ligand, and where in the domain the
                     residue sits. Shared by all ~19 substitutions at a
                     position.
    SUBSTITUTION      what the residue changed *into*: BLOSUM62, and the
                     change in hydropathy, volume and charge.

59% of the response spread in the BCR-ABL screen is *within* a single position
(see README), so position-level features are mathematically incapable of
explaining most of the variance no matter how good they are. The substitution
block has to carry it. That is why both blocks exist and why
`scripts/train_baseline.py` ablates them separately.

Everything here is stdlib — the property tables are embedded rather than pulled
from biopython, which keeps the install empty.
"""

from __future__ import annotations

from typing import Any

from pipeline.pocket import PocketError, pocket_distance_for

AA = "ARNDCQEGHILKMFPSTWYV"

# Kyte-Doolittle hydropathy. Positive = hydrophobic.
HYDROPATHY = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

# Side-chain volume, cubic angstroms (Zamyatnin).
VOLUME = {
    "A": 88.6, "R": 173.4, "N": 114.1, "D": 111.1, "C": 108.5, "Q": 143.8,
    "E": 138.4, "G": 60.1, "H": 153.2, "I": 166.7, "L": 166.7, "K": 168.6,
    "M": 162.9, "F": 189.9, "P": 112.7, "S": 89.0, "T": 116.1, "W": 227.8,
    "Y": 193.6, "V": 140.0,
}

CHARGE = {"D": -1.0, "E": -1.0, "K": 1.0, "R": 1.0, "H": 0.1}
POLAR = set("RNDQEHKSTY")
AROMATIC = set("FWYH")

# BLOSUM62, rows/cols in AA order above.
_B62_ROWS = """
 4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0
-1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3
-2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3
-2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3
 0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
-1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2
-1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2
 0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3
-2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3
-1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3
-1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1
-1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2
-1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1
-2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1
-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2
 1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2
 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0
-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3
-2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1
 0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4
"""
BLOSUM62: dict[tuple[str, str], int] = {}
for _i, _line in enumerate(_B62_ROWS.strip().splitlines()):
    for _j, _v in enumerate(_line.split()):
        BLOSUM62[(AA[_i], AA[_j])] = int(_v)


def blosum(ref: str, alt: str) -> int:
    """BLOSUM62 substitution score; 0 for anything non-standard."""
    return BLOSUM62.get((ref.upper(), alt.upper()), 0)


# --- feature blocks -------------------------------------------------------

SUBSTITUTION_NAMES = [
    "blosum62",
    "hyd_ref", "hyd_alt", "hyd_delta", "hyd_absdelta",
    "vol_ref", "vol_alt", "vol_delta", "vol_absdelta",
    "chg_ref", "chg_alt", "chg_delta", "chg_absdelta",
    "polar_ref", "polar_alt", "polar_gain",
    "arom_ref", "arom_alt", "arom_loss",
    "alt_is_pro", "alt_is_gly", "ref_is_gly", "alt_is_cys",
]

POSITION_NAMES = ["distance_angstrom", "contact_score", "domain_frac", "has_distance"]

# ESM2 masked marginals. `delta_ll` and `ll_alt` vary per substitution; `ll_ref`
# and `entropy` are position-level conservation measures.
ESM_NAMES = ["esm_delta_ll", "esm_ll_ref", "esm_ll_alt", "esm_entropy", "has_esm"]

FEATURE_NAMES = SUBSTITUTION_NAMES + POSITION_NAMES


def substitution_features(ref: str, alt: str) -> list[float]:
    """What the residue became. Varies between substitutions at one position."""
    ref, alt = ref.upper(), alt.upper()
    hr, ha = HYDROPATHY.get(ref, 0.0), HYDROPATHY.get(alt, 0.0)
    vr, va = VOLUME.get(ref, 0.0), VOLUME.get(alt, 0.0)
    cr, ca = CHARGE.get(ref, 0.0), CHARGE.get(alt, 0.0)
    pr, pa = float(ref in POLAR), float(alt in POLAR)
    ar, aa_ = float(ref in AROMATIC), float(alt in AROMATIC)
    return [
        float(blosum(ref, alt)),
        hr, ha, ha - hr, abs(ha - hr),
        vr, va, va - vr, abs(va - vr),
        cr, ca, ca - cr, abs(ca - cr),
        pr, pa, max(0.0, pa - pr),
        ar, aa_, max(0.0, ar - aa_),
        float(alt == "P"), float(alt == "G"), float(ref == "G"), float(alt == "C"),
    ]


def position_features(
    position: int,
    *,
    compound: str,
    ref: str | None = None,
    domain: tuple[int, int] = (242, 512),
    distance_cache: dict[int, float] | None = None,
) -> list[float]:
    """Where the residue sits. Constant across substitutions at one position.

    `has_distance` flags rows where the structure has no such residue, so a
    missing distance is a declared gap rather than a silent zero.
    """
    from pipeline.baseline import contact_score

    d: float | None = None
    if distance_cache is not None and position in distance_cache:
        d = distance_cache[position]
    else:
        try:
            d = pocket_distance_for(compound, position, expect_residue=ref)["distance_angstrom"]
        except (PocketError, KeyError):
            d = None
        if distance_cache is not None and d is not None:
            distance_cache[position] = d

    lo, hi = domain
    frac = (position - lo) / max(1, hi - lo)
    if d is None:
        # Median-ish filler, flagged by has_distance=0 so a model can learn to
        # distrust it instead of reading it as a real measurement.
        return [10.0, 0.0, frac, 0.0]
    return [d, contact_score(d), frac, 1.0]


def esm_features(
    position: int, ref: str, alt: str, marginals: dict[str, Any] | None
) -> list[float]:
    """ESM2 masked-marginal features for one substitution.

    `has_esm` flags an unscored position so a model can distrust the filler
    instead of reading zeros as a confident prediction.
    """
    import math

    row = ((marginals or {}).get("log_probs") or {}).get(str(position))
    if not row:
        return [0.0, 0.0, 0.0, 0.0, 0.0]
    lr = row.get(ref.upper())
    la = row.get(alt.upper())
    if lr is None or la is None:
        return [0.0, 0.0, 0.0, 0.0, 0.0]
    # Entropy over the 20 standard residues: high = tolerant position.
    ent = -sum(math.exp(v) * v for v in row.values())
    return [la - lr, lr, la, ent, 1.0]


def featurize(
    variants: list[Any],
    *,
    compound: str = "imatinib",
    domain: tuple[int, int] = (242, 512),
    blocks: str = "all",
    esm_marginals: dict[str, Any] | None = None,
) -> tuple[list[list[float]], list[float], list[int], list[str]]:
    """Build the design matrix.

    `blocks` selects an ablation: "substitution", "position", or "all". Returns
    (X, y, groups, feature_names) where `groups` is the residue position, for
    group-aware cross-validation.
    """
    valid = {"all", "substitution", "position", "esm", "all+esm",
             "substitution+esm", "position+esm"}
    if blocks not in valid:
        raise ValueError(f"unknown blocks={blocks!r}, expected one of {sorted(valid)}")
    want_sub = blocks in ("all", "substitution", "all+esm", "substitution+esm")
    want_pos = blocks in ("all", "position", "all+esm", "position+esm")
    want_esm = blocks.endswith("esm")
    if want_esm and esm_marginals is None:
        raise ValueError(f"blocks={blocks!r} needs esm_marginals")
    cache: dict[int, float] = {}
    X: list[list[float]] = []
    y: list[float] = []
    groups: list[int] = []
    for v in variants:
        row: list[float] = []
        if want_sub:
            row += substitution_features(v.ref_aa, v.alt_aa)
        if want_pos:
            row += position_features(
                v.position, compound=compound, ref=v.ref_aa,
                domain=domain, distance_cache=cache,
            )
        if want_esm:
            row += esm_features(v.position, v.ref_aa, v.alt_aa, esm_marginals)
        X.append(row)
        y.append(v.response)
        groups.append(v.position)
    names: list[str] = []
    if want_sub:
        names += SUBSTITUTION_NAMES
    if want_pos:
        names += POSITION_NAMES
    if want_esm:
        names += ESM_NAMES
    return X, y, groups, names
