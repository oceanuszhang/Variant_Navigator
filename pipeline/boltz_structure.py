"""Boltz structure prediction for the mutant protein.

STUB — depends on the Modal/GPU workstream. Returns a path to a predicted
structure plus pocket geometry that drug_screen.py consumes.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class Structure:
    path: Path                              # predicted mutant structure (mmCIF/PDB)
    sequence: str
    confidence: float                       # mean pLDDT, 0-100
    pocket_residues: list[int] = field(default_factory=list)
    per_residue_confidence: dict[int, float] = field(default_factory=dict)
    source: str = "boltz"

    @property
    def reliable(self) -> bool:
        return self.confidence >= 70.0


def predict_structure(
    sequence: str,
    *,
    name: str = "mutant",
    out_dir: Path | str = "output/structures",
    pocket_hint: list[int] | None = None,
) -> Structure:
    """Predict the mutant structure. Raises until Boltz is wired up."""
    raise NotImplementedError(
        "Boltz prediction not wired yet. Implement the Modal Boltz call here; "
        "pocket_hint carries UniProt binding-site residues from fetch.py."
    )


def available() -> bool:
    return bool(os.getenv("MODAL_TOKEN_ID") and os.getenv("MODAL_TOKEN_SECRET"))
