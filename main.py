#!/usr/bin/env python3
"""VariantNavigator — predict which drugs still work against a mutant target.

    python main.py --gene EGFR --variant T790M --mode oncology
    python main.py --test-case egfr_t790m
    python main.py --test-case all

Stages that are not wired up yet (ESM2, Boltz, Ptarmigan-1) are skipped with a
note rather than aborting the run, so the pipeline always produces a report.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

from pipeline import agent_report, boltz_structure, drug_screen, esm2_score
from pipeline.fetch import fetch_protein
from pipeline.mutant import VariantError, apply_variant

TEST_CASE_DIR = Path("data/test_cases")
DEFAULT_PANEL = {"oncology": "oncology", "agriculture": "fungicides"}


def load_env(path: str | Path = ".env") -> None:
    """Minimal .env loader — avoids a python-dotenv dependency."""
    p = Path(path)
    if not p.exists():
        return
    import os

    for line in p.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, _, value = line.partition("=")
        os.environ.setdefault(key.strip(), value.strip().strip("'\""))


def run(
    gene: str,
    variant: str,
    mode: str,
    *,
    panel: str | None = None,
    organism_id: int | None = None,
    expected: dict[str, str] | None = None,
    use_agent: bool = True,
) -> dict[str, Any]:
    panel = panel or DEFAULT_PANEL.get(mode, "oncology")
    print(f"\n{'=' * 66}\n{gene} {variant}  [{mode} / {panel} panel]\n{'=' * 66}")

    # 1 — wildtype sequence + binding annotations
    print("[1/5] UniProt fetch…", end=" ", flush=True)
    protein = fetch_protein(gene, organism_id)
    print(f"{protein.accession} ({protein.organism}), {protein.length} aa")

    # 2 — apply the variant
    print("[2/5] applying variant…", end=" ", flush=True)
    mutant = apply_variant(protein.sequence, variant)
    pocket = protein.binding_positions()
    in_site = any(p in pocket for p in mutant.positions)
    print(
        f"ok  ({', '.join(str(s) for s in mutant.substitutions) or variant}"
        f"{'; in annotated binding site' if in_site else ''})"
    )

    payload: dict[str, Any] = {
        "gene": protein.gene,
        "variant": variant,
        "mode": mode,
        "panel": panel,
        "uniprot": protein.accession,
        "organism": protein.organism,
        "sequence_length": protein.length,
        "binding_residues": pocket,
        "in_binding_site": in_site,
        "mutated_positions": mutant.positions,
        "expected": expected or {},
        "skipped_stages": [],
    }

    # 3 — ESM2 sequence-model score
    print("[3/5] ESM2 scoring…", end=" ", flush=True)
    try:
        score = esm2_score.score_variant(protein.sequence, mutant.sequence, variant)
        payload["esm2"] = score.__dict__
        print(f"Δll={score.delta_ll:+.3f}")
    except NotImplementedError as exc:
        payload["skipped_stages"].append({"stage": "esm2", "reason": str(exc)})
        print("skipped (not wired)")

    # 4 — mutant structure
    print("[4/5] Boltz structure…", end=" ", flush=True)
    structure = None
    try:
        structure = boltz_structure.predict_structure(
            mutant.sequence, name=f"{gene}_{variant}", pocket_hint=pocket
        )
        payload["structure"] = {
            "path": str(structure.path),
            "confidence": structure.confidence,
        }
        print(f"pLDDT={structure.confidence:.1f}")
    except NotImplementedError as exc:
        payload["skipped_stages"].append({"stage": "boltz", "reason": str(exc)})
        print("skipped (not wired)")

    # 5 — drug panel screen
    print("[5/5] drug screen…", end=" ", flush=True)
    try:
        screen = drug_screen.screen_panel(
            protein.gene,
            variant,
            protein.sequence,
            mutant.sequence,
            panel,
            pocket_residues=structure.pocket_residues if structure else pocket,
        )
        payload["screen"] = {
            "panel": screen.panel,
            "results": [r.__dict__ for r in screen.ranked()],
        }
        print(f"{len(screen.results)} drugs, {len(screen.still_works())} retained")
    except (NotImplementedError, RuntimeError) as exc:
        drugs = drug_screen.load_panel(panel)
        payload["skipped_stages"].append({"stage": "drug_screen", "reason": str(exc)})
        payload["screen"] = {"panel": panel, "results": []}
        print(f"skipped (not wired; {len(drugs)} drugs staged)")

    report = agent_report.write_report(payload, use_agent=use_agent)
    print(f"\nreport: {report.path}  [{report.generated_by}]")
    if payload["skipped_stages"]:
        stages = ", ".join(s["stage"] for s in payload["skipped_stages"])
        print(f"stages pending: {stages}")
    return payload


def run_test_case(name: str, *, use_agent: bool = True) -> dict[str, Any]:
    path = TEST_CASE_DIR / f"{name}.json"
    if not path.exists():
        raise FileNotFoundError(f"no test case {name!r} in {TEST_CASE_DIR}")
    case = json.loads(path.read_text())
    return run(
        case["gene"],
        case["variant"],
        case["mode"],
        panel=case.get("panel"),
        organism_id=case.get("organism_id"),
        expected=case.get("expected"),
        use_agent=use_agent,
    )


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gene", help="gene name or UniProt accession, e.g. EGFR / P00533")
    ap.add_argument("--variant", help="variant notation, e.g. T790M")
    ap.add_argument("--mode", choices=["oncology", "agriculture"], default="oncology")
    ap.add_argument("--panel", help="override drug panel (name or path)")
    ap.add_argument("--organism-id", type=int, help="NCBI taxon id, to disambiguate genes")
    ap.add_argument("--test-case", help="run a benchmark case by id, or 'all'")
    ap.add_argument("--no-agent", action="store_true", help="force the local report renderer")
    args = ap.parse_args(argv)

    load_env()
    use_agent = not args.no_agent

    try:
        if args.test_case == "all":
            for path in sorted(TEST_CASE_DIR.glob("*.json")):
                run_test_case(path.stem, use_agent=use_agent)
        elif args.test_case:
            run_test_case(args.test_case, use_agent=use_agent)
        elif args.gene and args.variant:
            run(
                args.gene,
                args.variant,
                args.mode,
                panel=args.panel,
                organism_id=args.organism_id,
                use_agent=use_agent,
            )
        else:
            ap.error("give --gene and --variant, or --test-case")
    except (VariantError, LookupError, FileNotFoundError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
