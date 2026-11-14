"""Pocket-distance baseline: the one predictive model that runs today.

No credentials, no GPU, no blocked dataset. It takes the one real feature we can
compute — minimum heavy-atom distance from the variant site to the bound
inhibitor (`pipeline.pocket`) — and turns it into a per-compound call.

WHAT THIS IS NOT. It is not a log2 fold change. Producing the matrix cell that
`figures/prediction_matrix.svg` advertises requires regressing on measured
dose-response data, and `data/datasets/` is still empty (see README). Fitting a
fold change with nothing to fit against would be inventing numbers. So this
module reports a dimensionless *contact score* and a qualitative call, and
`calibrated()` returns False until the datasets land.

WHAT IT IS GOOD FOR. A geometric prior, and a triage filter: a substitution 10 A
from the ligand is very unlikely to change binding directly, and that is worth
knowing before spending GPU budget. It is also a genuine baseline — any learned
model has to beat geometry alone to have earned its complexity.

WHERE IT PROVABLY FAILS. Distance is blind to chemistry. Osimertinib was
designed to tolerate the T790M gatekeeper (it binds C797 covalently), so
geometry says "Thr790 touches the drug, expect an effect" and the clinic says
"osimertinib works". `validate_against_controls()` reports that miss rather than
hiding it. Read the docstring there before quoting any accuracy number.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Literal

from pipeline.pocket import POCKET_REFERENCE, PocketError, pocket_distance_for

TARGETS_PATH = Path("data/targets.json")

Call = Literal["likely_affected", "uncertain", "likely_unaffected", "no_reference"]

# Distance at which a residue stops touching the ligand. Same 5.0 A cutoff the
# figures draw, which is the conventional heavy-atom contact distance.
CONTACT_CUTOFF = 5.0
# Softness of the transition, in angstroms. 1.0 A makes the score fall from
# ~0.88 at 3 A to ~0.12 at 7 A: confident at the extremes, hedged in between.
SOFTNESS = 1.0

AFFECTED_AT = 0.5
UNAFFECTED_BELOW = 0.2


def contact_score(distance_angstrom: float) -> float:
    """Map a distance to a 0-1 contact score, 1 = hard against the ligand.

    A logistic in distance centred on the contact cutoff. Monotone decreasing,
    so ranking by this score is just ranking by proximity — the score adds no
    information beyond the distance, it only puts it on a comparable scale.
    """
    z = (distance_angstrom - CONTACT_CUTOFF) / SOFTNESS
    # Guard the exponential against overflow for far-away residues.
    if z > 60:
        return 0.0
    if z < -60:
        return 1.0
    return 1.0 / (1.0 + math.exp(z))


def _call_from_score(score: float) -> Call:
    if score >= AFFECTED_AT:
        return "likely_affected"
    if score < UNAFFECTED_BELOW:
        return "likely_unaffected"
    return "uncertain"


def covered_compounds() -> list[str]:
    """Compounds with a reference structure, so a distance can be computed."""
    return sorted(POCKET_REFERENCE)


def compounds_for_target(accession: str) -> list[str]:
    """Compounds that have BOTH a reference structure and this target in scope.

    Without this the pipeline would try every reference compound against every
    target and ask an ALS structure for EGFR residue 790 — which raises, gets
    captured as "no_reference", and buries the real rows in noise.

    Reads data/targets.json so the pairing stays in the data file rather than
    being duplicated here. Unknown accession returns [] — the caller then
    reports the stage as not applicable instead of guessing.
    """
    try:
        with (TARGETS_PATH).open() as fh:
            targets = json.load(fh)
    except (OSError, json.JSONDecodeError):
        return []
    acc = accession.strip().upper()
    for entry in targets.values():
        if str(entry.get("uniprot", "")).upper() != acc:
            continue
        return [c for c in entry.get("compounds", []) if c.lower() in POCKET_REFERENCE]
    return []


def predict_cell(
    compound: str, position: int, *, expect_residue: str | None = None
) -> dict[str, Any]:
    """One variant x compound cell of the baseline prediction.

    `position` is in PUBLISHED numbering (clinical for EGFR, Arabidopsis
    convention for ALS); `pipeline.pocket` applies the per-structure offset.
    Pass `expect_residue` whenever the caller knows what should be there — it is
    what turns a numbering slip into an error instead of a wrong answer.
    """
    if compound.lower() not in POCKET_REFERENCE:
        return {
            "compound": compound,
            "published_position": position,
            "call": "no_reference",
            "reason": "no reference structure for this compound",
            "calibrated": False,
        }
    geom = pocket_distance_for(compound, position, expect_residue=expect_residue)
    score = contact_score(geom["distance_angstrom"])
    return {
        **geom,
        "contact_score": round(score, 3),
        "call": _call_from_score(score),
        # Loud on purpose: nothing downstream should read this as a fold change.
        "log2_fold_change": None,
        "calibrated": False,
    }


def predict_panel(
    positions: list[tuple[int, str]], compounds: list[str] | None = None
) -> list[dict[str, Any]]:
    """Baseline calls for every position x compound pair.

    `positions` is [(published_position, expected_one_letter_residue), ...].
    A per-cell failure is captured in the row rather than aborting the sweep, so
    one bad structure cannot take out the whole matrix.
    """
    compounds = compounds or covered_compounds()
    rows: list[dict[str, Any]] = []
    for pos, expect in positions:
        for compound in compounds:
            try:
                rows.append(predict_cell(compound, pos, expect_residue=expect))
            except PocketError as exc:
                rows.append({
                    "compound": compound,
                    "published_position": pos,
                    "call": "no_reference",
                    "reason": str(exc),
                    "calibrated": False,
                })
    return rows


def calibrated() -> bool:
    """False until the measured dose-response data exists to calibrate against.

    Kept as a function so callers gate on it instead of hardcoding optimism.
    """
    return False


# Ground truth for the EGFR controls, from the README benchmark table. Only the
# EGFR cells are listed: the ALS truth values are exactly what the missing
# datasets are supposed to supply, so asserting them here would be circular.
#
# "affected" means the variant degrades that drug's activity in the clinic.
CONTROL_TRUTH: list[tuple[str, int, str, bool, str]] = [
    ("erlotinib",   790, "T", True,
     "T790M is the classic gatekeeper escape from first-generation EGFR TKIs"),
    ("osimertinib", 790, "T", False,
     "osimertinib was designed to tolerate T790M - it binds Cys797 covalently"),
    ("erlotinib",   858, "L", False,
     "L858R is activating and sensitising, not resistance-conferring"),
    ("osimertinib", 858, "L", False,
     "L858R stays osimertinib-sensitive"),
]


def validate_against_controls() -> dict[str, Any]:
    """Score the baseline on the four EGFR control cells whose truth we know.

    The headline number here is deliberately small and deliberately imperfect.
    Four cells is not an evaluation set; it is a smoke test that the geometry is
    wired up and a demonstration of where distance-only reasoning breaks. Do not
    quote the hit rate as model accuracy.
    """
    rows = []
    hits = 0
    for compound, pos, expect, truth_affected, why in CONTROL_TRUTH:
        row: dict[str, Any] = {
            "compound": compound,
            "position": pos,
            "truth_affected": truth_affected,
            "rationale": why,
        }
        try:
            cell = predict_cell(compound, pos, expect_residue=expect)
            predicted_affected = cell["call"] == "likely_affected"
            row.update(
                distance_angstrom=cell["distance_angstrom"],
                contact_score=cell["contact_score"],
                call=cell["call"],
                predicted_affected=predicted_affected,
                correct=(predicted_affected == truth_affected),
                approximate=cell["approximate"],
            )
            hits += int(row["correct"])
        except PocketError as exc:
            row.update(error=str(exc), correct=False)
        rows.append(row)
    return {
        "cells": rows,
        "hits": hits,
        "total": len(CONTROL_TRUTH),
        "misses": [r for r in rows if not r.get("correct")],
        "calibrated": calibrated(),
    }


def _main() -> int:
    print("Pocket-distance baseline - validation on known EGFR controls")
    print("=" * 72)
    report = validate_against_controls()
    hdr = f"{'compound':13} {'pos':>5} {'dist A':>8} {'score':>6} {'predicted':>18} {'truth':>9}  ok"
    print(hdr)
    print("-" * 72)
    for r in report["cells"]:
        if "error" in r:
            print(f"{r['compound']:13} {r['position']:>5}  FAIL {r['error'][:40]}")
            continue
        print(
            f"{r['compound']:13} {r['position']:>5} {r['distance_angstrom']:>8.2f} "
            f"{r['contact_score']:>6.2f} {r['call']:>18} "
            f"{('affected' if r['truth_affected'] else 'unaffected'):>9}  "
            f"{'Y' if r['correct'] else 'N'}"
        )
    print("-" * 72)
    print(f"{report['hits']}/{report['total']} control cells reproduced "
          f"(calibrated={report['calibrated']})")
    for m in report["misses"]:
        print(f"\nMISS  {m['compound']} @ {m['position']}: {m['rationale']}")
        print("      Distance cannot see this: geometry says contact, chemistry "
              "says tolerated.")
    print("\nThis is a geometric prior, not a calibrated efficacy prediction.")
    print("No log2 fold change is produced until data/datasets/ is populated.")
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
