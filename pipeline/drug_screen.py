"""Ptarmigan-1 binding screen across a drug panel.

STUB — the Ptarmigan-1 API docs were not available in this working copy (see
README "Blocked / needs you"). The panel loading and ranking logic below is
real and testable; only `_call_ptarmigan` needs filling in.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal

PANEL_DIR = Path("data/drug_panels")
Verdict = Literal["work", "partial", "fail", "n/a"]


@dataclass
class DrugResult:
    name: str
    smiles: str
    wt_affinity: float | None = None        # predicted, kcal/mol (lower = tighter)
    mt_affinity: float | None = None
    delta_affinity: float | None = None     # mt - wt; positive = weaker binding
    verdict: Verdict = "n/a"
    confidence: float = 0.0
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass
class ScreenResult:
    gene: str
    variant: str
    panel: str
    results: list[DrugResult] = field(default_factory=list)

    def ranked(self) -> list[DrugResult]:
        """Best retained binders first; unscored drugs last."""
        return sorted(
            self.results,
            key=lambda r: (r.mt_affinity is None, r.mt_affinity or 0.0),
        )

    def still_works(self) -> list[DrugResult]:
        return [r for r in self.results if r.verdict == "work"]


def load_panel(panel: str) -> list[dict[str, Any]]:
    """Load a drug panel by name ('oncology' / 'fungicides') or path."""
    path = Path(panel)
    if not path.exists():
        path = PANEL_DIR / f"{panel}.json"
    if not path.exists():
        raise FileNotFoundError(f"no drug panel at {panel!r}")
    with path.open() as fh:
        drugs = json.load(fh)
    missing = [d.get("name", "?") for d in drugs if not d.get("smiles")]
    if missing:
        raise ValueError(f"panel {path.name} has entries without SMILES: {missing}")
    return drugs


def classify(delta_affinity: float, *, fail_at: float = 1.5, partial_at: float = 0.75) -> Verdict:
    """Turn an affinity shift into a work/partial/fail call.

    Thresholds are placeholders until we calibrate against the positive
    controls in data/test_cases/ — calibrate before trusting any verdict.
    """
    if delta_affinity >= fail_at:
        return "fail"
    if delta_affinity >= partial_at:
        return "partial"
    return "work"


def _call_ptarmigan(sequence: str, smiles: str, **kwargs: Any) -> dict[str, Any]:
    key = os.getenv("PTARMIGAN_API_KEY")
    if not key:
        raise RuntimeError("PTARMIGAN_API_KEY not set (see .env.example)")
    raise NotImplementedError(
        "Ptarmigan-1 request not implemented: API docs still needed."
    )


def screen_panel(
    gene: str,
    variant: str,
    wt_sequence: str,
    mt_sequence: str,
    panel: str,
    *,
    pocket_residues: list[int] | None = None,
) -> ScreenResult:
    """Score every drug in the panel against wildtype and mutant."""
    drugs = load_panel(panel)
    out = ScreenResult(gene=gene, variant=variant, panel=panel)
    for drug in drugs:
        wt = _call_ptarmigan(wt_sequence, drug["smiles"], pocket=pocket_residues)
        mt = _call_ptarmigan(mt_sequence, drug["smiles"], pocket=pocket_residues)
        delta = mt["affinity"] - wt["affinity"]
        out.results.append(
            DrugResult(
                name=drug["name"],
                smiles=drug["smiles"],
                wt_affinity=wt["affinity"],
                mt_affinity=mt["affinity"],
                delta_affinity=delta,
                verdict=classify(delta),
                confidence=mt.get("confidence", 0.0),
                metadata={k: drug[k] for k in ("target", "class", "approval") if k in drug},
            )
        )
    return out


def available() -> bool:
    return bool(os.getenv("PTARMIGAN_API_KEY"))
