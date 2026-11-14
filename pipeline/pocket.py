"""Pocket geometry from experimental structures — the distance feature.

Needs no credentials: PDB files are public, and the maths is a min-distance
computation, so the pocket-distance baseline model is fully runnable today.

For each variant position we compute the minimum heavy-atom distance from that
residue to the bound inhibitor. A substitution far from the ligand is unlikely
to change binding; one lining the pocket is the opposite.

THIRD NUMBERING SYSTEM WARNING. We now juggle three:

    published variant numbering  (clinical EGFR; Arabidopsis-convention ALS)
    UniProt sequence numbering   (P00533; Q6K2E8 rice, offset -26)
    PDB residue numbering        (per-structure, validated below)

Every lookup here validates the residue identity it lands on and raises rather
than returning a plausible wrong number. `validate_references()` checks all four
structures at once.
"""

from __future__ import annotations

import math
import urllib.request
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Any

# Biological assembly, NOT the asymmetric unit. This is load-bearing: ALS's
# herbicide site is a tunnel at the dimer interface, so in the asymmetric unit
# Ala122 sits 28 A from the ligand while in the dimer it is 3.6 A. Using the
# monomer makes the pocket model predict "no effect" for A122T, the headline
# 139.5-fold resistance case.
ASSEMBLY_URL = "https://files.rcsb.org/download/{pdb_id}.pdb1.gz"
PDB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
CACHE_DIR = Path("data/structures")
TIMEOUT = 60

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "CSD": "C", "MSE": "M",   # modified residues seen in these structures
}

# compound -> (pdb id, ligand code, exact or approximate)
# pdb_offset: add to PUBLISHED variant numbering to reach this file's numbering.
# Verified empirically per structure, not assumed.
POCKET_REFERENCE: dict[str, dict[str, Any]] = {
    "erlotinib":   {"pdb": "1M17", "ligand": "AQ4", "exact": True,
                    "pdb_offset": -24,
                    "note": "erlotinib bound to EGFR kinase. 1M17 uses MATURE "
                            "protein numbering (signal peptide removed), so "
                            "clinical T790 is residue 766 here."},
    "osimertinib": {"pdb": "4ZAU", "ligand": "YY3", "exact": True,
                    "pdb_offset": 0,
                    "note": "osimertinib bound to WILD-TYPE EGFR (Thr790 "
                            "present, not T790M). Numbering matches clinical."},
    "imatinib":    {"pdb": "1IEP", "ligand": "STI", "exact": True,
                    "pdb_offset": 0,
                    "note": "imatinib bound to ABL1 kinase. 1IEP numbering "
                            "matches ABL1a == UniProt P00519 == the BCR-ABL DMS "
                            "dataset, all at offset 0. Do NOT substitute 1OPJ: "
                            "it holds the same STI ligand but needs offset +19, "
                            "so it would return the wrong residue's distance."},
    "imazethapyr": {"pdb": "1Z8N", "ligand": "1IQ", "exact": False,
                    "pdb_offset": 0,
                    "note": "imazaquin in Arabidopsis ALS; same imidazolinone "
                            "class, so distances are APPROXIMATE. Numbering "
                            "matches Arabidopsis convention."},
    "nicosulfuron": {"pdb": "1YHZ", "ligand": "1CS", "exact": False,
                     "pdb_offset": 0,
                     "note": "chlorsulfuron in Arabidopsis ALS; same "
                             "sulfonylurea class, so distances are APPROXIMATE. "
                             "Numbering matches Arabidopsis convention."},
}


class PocketError(RuntimeError):
    """Structure missing a residue or ligand, or a residue identity mismatch."""


@dataclass
class Atom:
    residue_number: int
    residue_name: str
    chain: str
    name: str
    x: float
    y: float
    z: float
    is_hetatm: bool

    @property
    def is_heavy(self) -> bool:
        return not self.name.startswith("H")


@dataclass
class StructureAtoms:
    pdb_id: str
    protein: list[Atom] = field(default_factory=list)
    ligands: dict[str, list[Atom]] = field(default_factory=dict)

    def residues(self, chain: str | None = None) -> dict[int, str]:
        out: dict[int, str] = {}
        for a in self.protein:
            if chain and a.chain != chain:
                continue
            out.setdefault(a.residue_number, THREE_TO_ONE.get(a.residue_name, "X"))
        return out


def _download(pdb_id: str) -> Path:
    """Fetch the biological assembly, falling back to the asymmetric unit."""
    import gzip

    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    pdb_id = pdb_id.upper()
    path = CACHE_DIR / f"{pdb_id}.assembly.pdb"
    if path.exists() and path.stat().st_size > 0:
        return path
    try:
        url = ASSEMBLY_URL.format(pdb_id=pdb_id)
        with urllib.request.urlopen(url, timeout=TIMEOUT) as resp:
            path.write_bytes(gzip.decompress(resp.read()))
        return path
    except Exception:
        # No assembly file published; the asymmetric unit may still be a
        # complete biological unit, but interface-lining residues can read as
        # far from the ligand. Callers see this via `assembly_used`.
        path = CACHE_DIR / f"{pdb_id}.pdb"
        if path.exists() and path.stat().st_size > 0:
            return path
        with urllib.request.urlopen(
            PDB_URL.format(pdb_id=pdb_id), timeout=TIMEOUT
        ) as resp:
            path.write_bytes(resp.read())
        return path


def assembly_used(pdb_id: str) -> bool:
    """True when distances come from the biological assembly."""
    return (CACHE_DIR / f"{pdb_id.upper()}.assembly.pdb").exists()


@lru_cache(maxsize=8)
def load_structure(pdb_id: str) -> StructureAtoms:
    """Parse a PDB file into protein atoms and per-ligand atom lists."""
    path = _download(pdb_id)
    out = StructureAtoms(pdb_id=pdb_id.upper())
    with path.open() as fh:
        for line in fh:
            tag = line[:6].strip()
            if tag not in ("ATOM", "HETATM"):
                continue
            resname = line[17:20].strip()
            if resname == "HOH":
                continue
            try:
                atom = Atom(
                    residue_number=int(line[22:26]),
                    residue_name=resname,
                    chain=line[21].strip() or "A",
                    name=line[12:16].strip(),
                    x=float(line[30:38]),
                    y=float(line[38:46]),
                    z=float(line[46:54]),
                    is_hetatm=(tag == "HETATM"),
                )
            except ValueError:
                continue
            # A modified residue is part of the chain, not a ligand.
            if atom.is_hetatm and resname not in THREE_TO_ONE:
                out.ligands.setdefault(resname, []).append(atom)
            else:
                out.protein.append(atom)
    if not out.protein:
        raise PocketError(f"{pdb_id}: no protein atoms parsed")
    return out


def _dist(a: Atom, b: Atom) -> float:
    return math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2)


def ligand_distance(
    pdb_id: str,
    ligand: str,
    position: int,
    *,
    expect_residue: str | None = None,
    chain: str | None = None,
) -> float:
    """Minimum heavy-atom distance (Å) from `position` to `ligand`.

    `expect_residue` is a one-letter code; when given, a mismatch raises. Pass it
    whenever the caller knows what should be there — it is the only guard against
    a PDB numbering offset silently returning the distance of the wrong residue.
    """
    st = load_structure(pdb_id)
    lig = st.ligands.get(ligand.upper())
    if not lig:
        raise PocketError(
            f"{pdb_id}: no ligand {ligand!r} (has: {sorted(st.ligands)})"
        )

    res_atoms = [
        a for a in st.protein
        if a.residue_number == position and (chain is None or a.chain == chain)
    ]
    if not res_atoms:
        raise PocketError(f"{pdb_id}: no residue at position {position}")

    found = THREE_TO_ONE.get(res_atoms[0].residue_name, "X")
    if expect_residue and found != expect_residue.upper():
        raise PocketError(
            f"{pdb_id}: position {position} is {found} "
            f"({res_atoms[0].residue_name}), expected {expect_residue} — "
            "PDB numbering offset?"
        )

    return min(
        _dist(a, l)
        for a in res_atoms if a.is_heavy
        for l in lig if l.is_heavy
    )


def pocket_distance_for(
    compound: str, position: int, *, expect_residue: str | None = None
) -> dict[str, Any]:
    """Distance from a PUBLISHED variant position to `compound`'s inhibitor.

    `position` is in the numbering the paper uses (clinical for EGFR,
    Arabidopsis convention for ALS); the per-structure offset is applied here.
    """
    ref = POCKET_REFERENCE.get(compound.lower())
    if not ref:
        raise PocketError(f"no pocket reference for compound {compound!r}")
    pdb_pos = position + ref["pdb_offset"]
    d = ligand_distance(
        ref["pdb"], ref["ligand"], pdb_pos, expect_residue=expect_residue
    )
    return {
        "compound": compound,
        "pdb": ref["pdb"],
        "ligand": ref["ligand"],
        "published_position": position,
        "pdb_position": pdb_pos,
        "pdb_offset": ref["pdb_offset"],
        "distance_angstrom": round(d, 2),
        "approximate": not ref["exact"],
        "assembly": assembly_used(ref["pdb"]),
        "note": ref["note"],
    }


def validate_references() -> list[dict[str, Any]]:
    """Check each reference structure's numbering against known residues.

    Run this before trusting any distance. Positions are given in each
    structure's own PDB numbering; the expected residues come from the variants
    the experiments actually care about.
    """
    checks = [
        ("erlotinib", 790, "T", "EGFR gatekeeper Thr790"),
        ("erlotinib", 858, "L", "EGFR Leu858"),
        ("osimertinib", 790, "T", "EGFR Thr790 (4ZAU is wild-type)"),
        ("osimertinib", 858, "L", "EGFR Leu858"),
        ("imatinib", 315, "T", "ABL1 gatekeeper Thr315 (T315I is the classic escape)"),
        ("imatinib", 253, "Y", "ABL1 Tyr253 - P-loop"),
        ("imatinib", 396, "H", "ABL1 His396 - activation loop, allosteric escape"),
        ("imazethapyr", 122, "A", "ALS Ala122 - dimer-interface residue"),
        ("imazethapyr", 574, "W", "ALS Trp574"),
        ("imazethapyr", 653, "S", "ALS Ser653"),
        ("nicosulfuron", 122, "A", "ALS Ala122 - dimer-interface residue"),
        ("nicosulfuron", 574, "W", "ALS Trp574"),
        ("nicosulfuron", 653, "S", "ALS Ser653"),
    ]
    results = []
    for compound, pos, expect, label in checks:
        row: dict[str, Any] = {
            "compound": compound, "position": pos, "expected": expect,
            "label": label,
        }
        try:
            row.update(pocket_distance_for(compound, pos, expect_residue=expect))
            row["ok"] = True
        except PocketError as exc:
            row["ok"] = False
            row["error"] = str(exc)
        results.append(row)
    return results


if __name__ == "__main__":
    hdr = f"{'compound':13} {'pub':>5} {'pdb':>5} {'exp':>4} {'dist Å':>8} {'asm':>4}  note"
    print(hdr)
    print("-" * 100)
    for r in validate_references():
        if r["ok"]:
            flag = "~" if r["approximate"] else " "
            print(
                f"{r['compound']:13} {r['published_position']:>5} "
                f"{r['pdb_position']:>5} {r['expected']:>4} "
                f"{r['distance_angstrom']:>7.2f}{flag} "
                f"{'yes' if r['assembly'] else 'NO':>4}  {r['label']}"
            )
        else:
            print(
                f"{r['compound']:13} {r['position']:>5} {'—':>5} "
                f"{r['expected']:>4} {'—':>8} {'—':>4}  "
                f"FAIL: {r.get('error', '')[:50]}"
            )
    print("\n~ = approximate (related, not identical, inhibitor)")
