"""Apply protein variant notation to a wildtype sequence.

Supports the notations we need for the demos:
    T790M        single substitution (1-letter)
    p.Thr790Met  3-letter, with or without the p. prefix
    L858R/T790M  multiple substitutions, slash- or comma-separated
    del746_750   in-frame deletion (EGFR exon 19 style)
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Iterable

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
AA = set(THREE_TO_ONE.values())

_SUB_1 = re.compile(r"^([A-Z])(\d+)([A-Z])$")
_SUB_3 = re.compile(r"^([A-Za-z]{3})(\d+)([A-Za-z]{3})$")
_DEL = re.compile(r"^(?:[A-Z]?(\d+))(?:_[A-Z]?(\d+))?del$|^del(\d+)(?:_(\d+))?$")


class VariantError(ValueError):
    """Variant notation that cannot be parsed, or disagrees with the sequence."""


@dataclass(frozen=True)
class Substitution:
    wt: str
    position: int
    mt: str

    def __str__(self) -> str:
        return f"{self.wt}{self.position}{self.mt}"


@dataclass(frozen=True)
class Deletion:
    start: int
    end: int

    def __str__(self) -> str:
        return f"{self.start}_{self.end}del" if self.end != self.start else f"{self.start}del"


@dataclass
class Mutant:
    wildtype: str
    sequence: str
    variant: str
    substitutions: list[Substitution] = field(default_factory=list)
    deletions: list[Deletion] = field(default_factory=list)

    @property
    def positions(self) -> list[int]:
        return sorted(
            [s.position for s in self.substitutions]
            + [p for d in self.deletions for p in range(d.start, d.end + 1)]
        )

    def window(self, flank: int = 15) -> str:
        """Sequence slice around the first affected residue, for quick eyeballing."""
        if not self.positions:
            return ""
        p = self.positions[0]
        lo, hi = max(0, p - 1 - flank), min(len(self.sequence), p + flank)
        return self.sequence[lo:hi]


def _norm_aa(token: str) -> str:
    t = token.upper()
    if len(t) == 1:
        if t not in AA:
            raise VariantError(f"{token!r} is not an amino acid")
        return t
    if t in THREE_TO_ONE:
        return THREE_TO_ONE[t]
    raise VariantError(f"{token!r} is not an amino acid")


def parse_variant(notation: str) -> list[Substitution | Deletion]:
    """Parse one or more edits out of a variant string."""
    text = notation.strip()
    if not text:
        raise VariantError("empty variant notation")

    edits: list[Substitution | Deletion] = []
    for raw in re.split(r"[/,;+]| and ", text):
        token = raw.strip()
        if not token:
            continue
        if token.lower().startswith("p."):
            token = token[2:]
        token = token.strip()

        if m := _SUB_1.match(token.upper()):
            edits.append(
                Substitution(_norm_aa(m.group(1)), int(m.group(2)), _norm_aa(m.group(3)))
            )
            continue
        if m := _SUB_3.match(token):
            edits.append(
                Substitution(_norm_aa(m.group(1)), int(m.group(2)), _norm_aa(m.group(3)))
            )
            continue
        if m := _DEL.match(token):
            start = m.group(1) or m.group(3)
            end = m.group(2) or m.group(4) or start
            edits.append(Deletion(int(start), int(end)))
            continue
        raise VariantError(f"cannot parse variant {raw.strip()!r}")

    if not edits:
        raise VariantError(f"no edits found in {notation!r}")
    return edits


def apply_variant(
    wildtype: str,
    notation: str,
    *,
    strict: bool = True,
) -> Mutant:
    """Return the mutant sequence for `notation` applied to `wildtype`.

    strict=True (default) raises if the stated wildtype residue does not match
    the sequence — that mismatch usually means a wrong isoform or offset, which
    would silently poison every downstream score.
    """
    wt = wildtype.strip().upper().replace("\n", "")
    if not wt or set(wt) - AA - {"X", "U", "B", "Z", "O"}:
        bad = sorted(set(wt) - AA - {"X", "U", "B", "Z", "O"})
        raise VariantError(f"wildtype sequence has non-residue characters: {bad}")

    edits = parse_variant(notation)
    chars = list(wt)
    subs: list[Substitution] = []
    dels: list[Deletion] = []

    for edit in edits:
        if isinstance(edit, Substitution):
            if not 1 <= edit.position <= len(chars):
                raise VariantError(
                    f"{edit} position outside sequence of length {len(chars)}"
                )
            observed = chars[edit.position - 1]
            if observed != edit.wt:
                msg = (
                    f"{edit}: sequence has {observed} at position {edit.position}, "
                    f"not {edit.wt}"
                )
                if strict:
                    raise VariantError(msg)
            chars[edit.position - 1] = edit.mt
            subs.append(edit)
        else:
            if not 1 <= edit.start <= edit.end <= len(chars):
                raise VariantError(f"{edit} outside sequence of length {len(chars)}")
            dels.append(edit)

    # Deletions applied last, high-to-low, so earlier positions stay 1-indexed
    # against the wildtype numbering the user typed.
    for d in sorted(dels, key=lambda d: d.start, reverse=True):
        del chars[d.start - 1 : d.end]

    return Mutant(
        wildtype=wt,
        sequence="".join(chars),
        variant=notation.strip(),
        substitutions=subs,
        deletions=dels,
    )


def mutate_positions(wildtype: str, edits: Iterable[Substitution]) -> str:
    """Low-level helper: apply pre-parsed substitutions, no validation."""
    chars = list(wildtype)
    for e in edits:
        chars[e.position - 1] = e.mt
    return "".join(chars)


if __name__ == "__main__":
    import sys

    from pipeline.fetch import fetch_protein

    gene = sys.argv[1] if len(sys.argv) > 1 else "EGFR"
    variant = sys.argv[2] if len(sys.argv) > 2 else "T790M"
    protein = fetch_protein(gene)
    m = apply_variant(protein.sequence, variant)
    print(f"{protein.accession} {gene} {variant}")
    print(f"wt  len={len(m.wildtype)}  mut len={len(m.sequence)}")
    print(f"edits: {[str(s) for s in m.substitutions] + [str(d) for d in m.deletions]}")
    print(f"window: {m.window()}")
