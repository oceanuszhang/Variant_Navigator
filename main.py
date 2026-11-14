#!/usr/bin/env python3
"""Run the minimal BCR-ABL imatinib growth-prediction experiment."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

<<<<<<< Updated upstream
from pipeline import agent_report, boltz_structure, drug_screen, esm2_score
from pipeline.fetch import fetch_protein
from pipeline.mutant import VariantError, apply_variant
=======
from pipeline.train_heads import train_and_evaluate
>>>>>>> Stashed changes

ROOT = Path(__file__).resolve().parent
DEFAULT_SEQUENCE = ROOT / "data/cancer-exp/ABL1_NM_005157.6_DMS_reference_aa.txt"
DEFAULT_DATA = ROOT / "data/datasets/abl_dms_netgrowth_1200nM.csv"
DEFAULT_OUTPUT = ROOT / "output/cancer-exp"


def generate_embeddings(
    sequence_path: Path,
    data_path: Path,
    embeddings_path: Path,
    batch_size: int,
) -> None:
    command = [
        sys.executable,
        "-m",
        "modal",
        "run",
        str(ROOT / "pipeline/modal_esm.py"),
        "--sequence-path",
        str(sequence_path),
        "--csv-path",
        str(data_path),
        "--output-path",
        str(embeddings_path),
        "--batch-size",
        str(batch_size),
    ]
    subprocess.run(command, cwd=ROOT, check=True)


<<<<<<< Updated upstream
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
) -> dict[str, Any]:
    panel = panel or DEFAULT_PANEL.get(mode, "oncology")
    print(f"\n{'=' * 66}\n{gene} {variant}  [{mode} / {panel} panel]\n{'=' * 66}")

    # 1 — wildtype sequence + binding annotations
    print("[1/5] UniProt fetch…", end=" ", flush=True)
    protein = fetch_protein(gene, organism_id)
    print(f"{protein.accession} ({protein.organism}), {protein.length} aa")

    # 2 — apply the variant
    print("[2/5] applying variant…", end=" ", flush=True)
    mutant = apply_variant(protein.sequence, variant, offset=offset)
    pocket = protein.binding_positions()
    in_site = any(p in pocket for p in mutant.positions)
    renumbered = (
        f" [clinical {variant} -> UniProt "
        f"{', '.join(str(s) for s in mutant.substitutions)}]" if offset else ""
=======
def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sequence", type=Path, default=DEFAULT_SEQUENCE)
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument(
        "--recompute-embeddings",
        action="store_true",
        help="rerun ESM-2 even when the cached embedding file exists",
>>>>>>> Stashed changes
    )
    parser.add_argument(
        "--train-only",
        action="store_true",
        help="train from an existing embeddings.npz without calling Modal",
    )
    args = parser.parse_args()

    sequence_path = args.sequence.resolve()
    data_path = args.data.resolve()
    output_dir = args.output_dir.resolve()
    embeddings_path = output_dir / "embeddings.npz"
    output_dir.mkdir(parents=True, exist_ok=True)

<<<<<<< Updated upstream
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
        use_agent=use_agent,
    )


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gene", help="gene name or UniProt accession, e.g. EGFR / P00533")
    ap.add_argument("--variant", help="variant notation, e.g. T790M")
    ap.add_argument(
        "--mode",
        choices=["oncology", "agriculture", "antimicrobial"],
        default="oncology",
    )
    ap.add_argument("--panel", help="override drug panel (name or path)")
    ap.add_argument("--organism-id", type=int, help="NCBI taxon id, to disambiguate genes")
    ap.add_argument(
        "--offset",
        type=int,
        default=0,
        help="shift variant numbering into UniProt numbering (rpoB needs +6)",
    )
    ap.add_argument("--test-case", help="run a benchmark case by id, or 'all'")
    ap.add_argument("--no-agent", action="store_true", help="force the local report renderer")
    args = ap.parse_args(argv)

    load_env()
    use_agent = not args.no_agent

    try:
        if args.test_case == "all":
            failures = []
            for path in sorted(TEST_CASE_DIR.glob("*.json")):
                try:
                    run_test_case(path.stem, use_agent=use_agent)
                except (VariantError, LookupError, ValueError) as exc:
                    failures.append((path.stem, exc))
                    print(f"  FAILED: {exc}")
            print(f"\n{'=' * 66}")
            total = len(list(TEST_CASE_DIR.glob('*.json')))
            print(f"{total - len(failures)}/{total} cases ran")
            for name, exc in failures:
                print(f"  {name}: {exc}")
            return 1 if failures else 0
        elif args.test_case:
            run_test_case(args.test_case, use_agent=use_agent)
        elif args.gene and args.variant:
            run(
                args.gene,
                args.variant,
                args.mode,
                panel=args.panel,
                organism_id=args.organism_id,
                offset=args.offset,
                use_agent=use_agent,
            )
        else:
            ap.error("give --gene and --variant, or --test-case")
    except (VariantError, LookupError, FileNotFoundError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
=======
    if not sequence_path.exists():
        parser.error(f"reference sequence not found: {sequence_path}")
    if not data_path.exists():
        parser.error(f"experiment CSV not found: {data_path}")

    if args.train_only and not embeddings_path.exists():
        parser.error(f"--train-only requires {embeddings_path}")

    if not args.train_only and (args.recompute_embeddings or not embeddings_path.exists()):
        print("[1/2] Generating frozen ESM-2 embeddings on Modal", flush=True)
        generate_embeddings(sequence_path, data_path, embeddings_path, args.batch_size)
    else:
        print(f"[1/2] Reusing {embeddings_path}", flush=True)

    print("[2/2] Training Ridge and XGBoost prediction heads", flush=True)
    metrics = train_and_evaluate(embeddings_path, output_dir)

    for model_name in ("ridge", "xgboost"):
        result = metrics["models"][model_name]
        print(
            f"{model_name}: "
            f"{result['in_range_percent']:.2f}% in replicate range; "
            f"MAE={result['mean_absolute_error']:.6f}"
        )
    print(f"Results written to {output_dir}")
>>>>>>> Stashed changes
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
