#!/usr/bin/env python3
"""Cancer variant assay response and structural research workflow.

Offline demo: python main.py --gene ABL1 --variant T315I --offline --no-agent
Online structures: add --boltz to opt in to budget-capped ABL1/imatinib calls.
Optional stages report skips; invalid input fails before prediction.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from dataclasses import asdict
import subprocess
import hashlib
from typing import Any

from pipeline import (
    agent_report,
    baseline,
    predict,
    boltz_structure,
    esm2_score,
)
from pipeline.fetch import fetch_protein, Protein
from pipeline.mutant import VariantError, apply_variant

TEST_CASE_DIR = Path("data/test_cases")
SUPPORTED_TARGETS = {"ABL1": "P00519", "P00519": "P00519", "EGFR": "P00533", "P00533": "P00533"}


def oncology_cases() -> list[Path]:
    """Return active cancer controls, excluding historical other domains."""
    return [p for p in sorted(TEST_CASE_DIR.glob("*.json"))
            if json.loads(p.read_text(encoding="utf-8")).get("mode") == "oncology"]



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
    offset: int = 0,
    use_agent: bool = True,
    offline: bool = False,
    run_boltz: bool = False,
    esm_run_dir: Path | None = None,
) -> dict[str, Any]:
    if mode != "oncology":
        raise ValueError("Only oncology is supported")
    accession = SUPPORTED_TARGETS.get(gene.strip().upper())
    if accession is None:
        raise ValueError("Supported targets: ABL1/P00519 and EGFR/P00533")
    if organism_id not in (None, 9606) or offset != 0:
        raise ValueError("Supported targets use human (9606), offset 0 numbering")
    panel = panel or "oncology"
    if panel != "oncology":
        raise ValueError("Only the oncology panel is supported")
    print(f"\n{'=' * 66}\n{gene} {variant}  [{mode} / {panel} panel]\n{'=' * 66}")

    # 1 — wildtype sequence + binding annotations
    print("[1/6] UniProt fetch…", end=" ", flush=True)
    if offline:
        if accession != "P00519":
            raise ValueError("Offline sequence reference is available only for ABL1")
        seq = Path("data/cancer-exp/ABL1_NM_005157.6_DMS_reference_aa.txt").read_text().strip()
        protein = Protein("P00519", "ABL1", "ABL1 study reference", "Homo sapiens", seq)
    else:
        protein = fetch_protein(accession, organism_id)
    if protein.organism != "Homo sapiens":
        raise ValueError("Only human cancer targets are supported")
    print(f"{protein.accession} ({protein.organism}), {protein.length} aa")

    # 2 — apply the variant
    print("[2/6] applying variant…", end=" ", flush=True)
    mutant = apply_variant(protein.sequence, variant, offset=offset)
    pocket = protein.binding_positions()
    in_site = any(p in pocket for p in mutant.positions)
    renumbered = (
        f" [clinical {variant} -> UniProt "
        f"{', '.join(str(s) for s in mutant.substitutions)}]" if offset else ""
    )
    print(
        f"ok  ({', '.join(str(s) for s in mutant.substitutions) or variant}"
        f"{'; in annotated binding site' if in_site else ''}){renumbered}"
    )

    payload: dict[str, Any] = {
        "gene": protein.gene,
        "variant": variant,
        "mode": mode,
        "panel": panel,
        "uniprot": protein.accession,
        "organism": protein.organism,
        "sequence_length": protein.length,
        "sequence_sha256": hashlib.sha256(protein.sequence.encode("ascii")).hexdigest(),
        "binding_residues": pocket,
        "in_binding_site": in_site,
        "mutated_positions": mutant.positions,
        "numbering_offset": offset,
        "expected": expected or {},
        "skipped_stages": [],
        "schema_version": 1,
        "offline": offline,
        "sequence_source": "518-aa study reference" if offline else "UniProt",
    }

    # 3 — ESM2 sequence-model score
    print("[3/6] ESM2 scoring…", end=" ", flush=True)
    try:
        if offline:
            raise esm2_score.ESMUnavailable("disabled in offline mode")
        if not esm2_score.available():
            raise esm2_score.ESMUnavailable(
                "torch / fair-esm not installed; use .venv/bin/python"
            )
        score = esm2_score.score_variant(
            protein.sequence, mutant.sequence, variant,
            accession=protein.accession,
        )
        payload["esm2"] = score.__dict__
        print(f"Δll={score.delta_ll:+.3f} ({score.model})")
    except (NotImplementedError, esm2_score.ESMUnavailable, OSError) as exc:
        payload["skipped_stages"].append({"stage": "esm2", "reason": str(exc)})
        print(f"skipped ({str(exc)[:60]})")

    # 4 — mutant structure
    print("[4/6] Boltz structure…", end=" ", flush=True)
    structure = None
    try:
        if offline or not run_boltz:
            raise boltz_structure.BoltzUnavailable("disabled; opt in with --boltz (online only)")
        if accession != "P00519":
            raise boltz_structure.BoltzUnavailable("Boltz crop is validated only for ABL1/imatinib")
        if mutant.deletions or any(not 242 <= p <= 512 for p in mutant.positions):
            raise boltz_structure.BoltzUnavailable("Boltz supports ABL1 substitutions at 242-512 only")
        structure = boltz_structure.predict_structure(
            mutant.sequence, name=f"{protein.gene}_{variant}", pocket_hint=pocket, compound="imatinib"
        )
        payload["structure"] = {**asdict(structure), "path": str(structure.path)}
        payload["structure"].pop("sequence", None)
        print(f"pLDDT={structure.confidence:.1f}")
    except (NotImplementedError, boltz_structure.BoltzUnavailable,
            boltz_structure.BoltzTooExpensive, OSError, subprocess.TimeoutExpired, ValueError, KeyError) as exc:
        payload["skipped_stages"].append({"stage": "boltz", "reason": str(exc)})
        print(f"skipped ({exc})")

    # 5 — pocket-distance baseline. Runs on public structures with no API key,
    # so this is the stage that keeps the demo predictive while Steps 3, 4 and 6
    # wait on credentials. Published numbering = UniProt positions minus the
    # offset we applied above.
    print("[5/6] pocket baseline…", end=" ", flush=True)
    try:
        if offline:
            raise RuntimeError("reference-structure retrieval disabled in offline mode")
        compounds = baseline.compounds_for_target(protein.accession)
        if not compounds:
            raise RuntimeError(
                f"no reference-structure compounds for {protein.accession}"
            )
        published = [
            (s.position - offset, s.wt) for s in mutant.substitutions
        ]
        rows = baseline.predict_panel(published, compounds)
        payload["baseline"] = {
            "calibrated": baseline.calibrated(),
            "contact_cutoff_angstrom": baseline.CONTACT_CUTOFF,
            "cells": rows,
        }
        scored = [r for r in rows if r.get("call") != "no_reference"]
        affected = [r for r in scored if r["call"] == "likely_affected"]
        if scored:
            print(
                f"{len(scored)} cells; "
                f"{len(affected)} likely affected "
                f"({', '.join(sorted({r['compound'] for r in affected})) or 'none'})"
            )
        else:
            print("no cells scored (no residue match in reference structures)")
    except (RuntimeError, OSError) as exc:
        payload["skipped_stages"].append({"stage": "baseline", "reason": str(exc)})
        print(f"skipped ({exc})")

    # 5b — trained net-growth model. Only ABL1/imatinib has measured data, so
    # Only ABL1/imatinib has response data.
    if protein.accession.upper() == "P00519":
        print("[5b/6] trained model…", end=" ", flush=True)
        try:
            subs = mutant.substitutions
            if len(subs) != 1 or mutant.deletions:
                raise RuntimeError("trained model covers single substitutions only")
            s = subs[0]
            pred = predict.predict_variant(s.position, s.wt, s.mt, allow_inference=not offline and esm_run_dir is None)
            payload["prediction"] = pred
            try:
                payload["validation"] = predict.held_out(str(s))
            except (OSError, ValueError, KeyError, TypeError) as exc:
                payload["skipped_stages"].append({"stage": "validation", "reason": str(exc)})
            val = pred.get("predicted_netgrowth", pred.get("measured_netgrowth"))
            kind = "predicted" if pred["predicted"] else "MEASURED in screen"
            print(f"netgrowth={val:+.5f} ({kind})")
        except (RuntimeError, OSError, ValueError) as exc:
            payload["skipped_stages"].append({"stage": "model", "reason": str(exc)})
            print(f"skipped ({exc})")

    else:
        payload["skipped_stages"].append({"stage": "model", "reason": "Response model covers ABL1/imatinib only"})

    if esm_run_dir is not None:
        try:
            if offline or accession != "P00519":
                raise ValueError("Explicit ESM inference requires online ABL1 input")
            from pipeline import embedding_model
            payload["esm_prediction"] = embedding_model.predict_variant(esm_run_dir, variant)
        except (RuntimeError, ValueError, OSError, ImportError, KeyError) as exc:
            payload["skipped_stages"].append({"stage": "esm_inference", "reason": str(exc)})

    reason = "Ptarmigan supports wild-type catalog targets, not mutant binding prediction"
    payload["screen"] = {"panel": panel, "results": [], "status": "unsupported", "reason": reason}
    payload["skipped_stages"].append({"stage": "drug_screen", "reason": reason})
    payload["status"] = "partial" if payload["skipped_stages"] else "complete"

    report = agent_report.write_report(payload, use_agent=use_agent)
    print(f"\nreport: {report.path}  [{report.generated_by}]")
    if payload["skipped_stages"]:
        stages = ", ".join(s["stage"] for s in payload["skipped_stages"])
        print(f"stages pending: {stages}")
    return payload


def run_test_case(name: str, *, use_agent: bool = True, offline: bool = False, run_boltz: bool = False, esm_run_dir: Path | None = None) -> dict[str, Any]:
    path = TEST_CASE_DIR / f"{name}.json"
    if not path.exists():
        raise FileNotFoundError(f"no test case {name!r} in {TEST_CASE_DIR}")
    case = json.loads(path.read_text())
    # An explicit accession beats gene-name search: several of these targets have
    # multiple alleles in UniProt and only some carry the residue we need.
    return run(
        case.get("uniprot") or case["gene"],
        case["variant"],
        case["mode"],
        panel=case.get("panel"),
        organism_id=case.get("organism_id"),
        expected=case.get("expected"),
        offset=case.get("numbering_offset", 0),
        use_agent=use_agent, offline=offline, run_boltz=run_boltz, esm_run_dir=esm_run_dir,
    )


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gene", help="gene name or UniProt accession, e.g. EGFR / P00533")
    ap.add_argument("--variant", help="variant notation, e.g. T790M")
    ap.add_argument(
        "--mode",
        choices=["oncology"],
        default="oncology",
    )
    ap.add_argument("--panel", help="override drug panel (name or path)")
    ap.add_argument("--organism-id", type=int, help="NCBI taxon id, to disambiguate genes")
    ap.add_argument(
        "--offset",
        type=int,
        default=0,
        help="supported human targets require offset 0",
    )
    ap.add_argument("--test-case", help="run a benchmark case by id, or 'all'")
    ap.add_argument("--no-agent", action="store_true", help="force the local report renderer")
    ap.add_argument("--offline", action="store_true", help="ABL1 study reference and cached assay; no network or inference")
    ap.add_argument("--boltz", action="store_true", help="opt in to budget-capped Boltz ABL1/imatinib calls")
    ap.add_argument("--esm-run-dir", type=Path, help="trusted saved ESM model directory; opt in to one Modal embedding")
    args = ap.parse_args(argv)

    load_env()
    use_agent = not args.no_agent

    try:
        if args.test_case == "all":
            failures = []
            for path in oncology_cases():
                try:
                    run_test_case(path.stem, use_agent=use_agent, offline=args.offline, run_boltz=args.boltz, esm_run_dir=args.esm_run_dir)
                except (VariantError, LookupError, ValueError, OSError) as exc:
                    failures.append((path.stem, exc))
                    print(f"  FAILED: {exc}")
            print(f"\n{'=' * 66}")
            total = len(oncology_cases())
            print(f"{total - len(failures)}/{total} cases ran")
            for name, exc in failures:
                print(f"  {name}: {exc}")
            return 1 if failures else 0
        elif args.test_case:
            run_test_case(args.test_case, use_agent=use_agent, offline=args.offline, run_boltz=args.boltz, esm_run_dir=args.esm_run_dir)
        elif args.gene and args.variant:
            run(
                args.gene,
                args.variant,
                args.mode,
                panel=args.panel,
                organism_id=args.organism_id,
                offset=args.offset,
                use_agent=use_agent, offline=args.offline, run_boltz=args.boltz, esm_run_dir=args.esm_run_dir,
            )
        else:
            ap.error("give --gene and --variant, or --test-case")
    except (VariantError, LookupError, FileNotFoundError, ValueError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
