"""Fixed four-complex ABL1/imatinib contrast; estimate-only unless --run.

This deliberately small, selected panel tests contrasts, not general accuracy.
Binding confidence is not affinity, and its delta is never assigned energy units.
"""
from __future__ import annotations
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import subprocess
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pipeline import boltz_structure as bs, dms
from pipeline.mutant import apply_variant

PANEL = ("WILDTYPE", "T315I", "T315A", "H396P")
OUT = Path("output/structures/boltz_benchmark")
REFERENCE = Path("data/cancer-exp/ABL1_NM_005157.6_DMS_reference_aa.txt")


def plan() -> list[dict]:
    """Validate variants and freeze identical settings before submitting."""
    reference = REFERENCE.read_text().strip()
    assay = {v.mutant: v for v in dms.load_netgrowth(path=dms.NETGROWTH_CSV)}
    jobs = []
    for variant in PANEL:
        sequence = reference if variant == "WILDTYPE" else apply_variant(reference, variant).sequence
        body = bs.build_input(bs.crop_sequence(sequence), compound="imatinib", num_samples=5)
        measured = assay.get(variant)
        if variant != "WILDTYPE" and measured is None:
            raise ValueError(f"Missing measurement for {variant}")
        jobs.append({"variant": variant, "sequence": sequence, "body": body,
                     "idempotency_key": bs.idempotency_key(body),
                     "measured_netgrowth": measured.response if measured else None,
                     "replicate_1": measured.rep1 if measured else None,
                     "replicate_2": measured.rep2 if measured else None})
    return jobs


def summarize(results: dict) -> dict:
    """Report prespecified mutant contrasts without converting confidence to efficacy."""
    rows = []
    wt = results.get("WILDTYPE", {}).get("binding_confidence")
    for variant in PANEL:
        row = dict(results.get(variant, {"error": "not completed"}))
        row["variant"] = variant
        binding = row.get("binding_confidence")
        row["delta_binding_confidence_vs_wt"] = binding-wt if binding is not None and wt is not None else None
        rows.append(row)
    control = results.get("T315A", {}).get("binding_confidence")
    contrasts = []
    for variant in ("T315I", "H396P"):
        binding = results.get(variant, {}).get("binding_confidence")
        contrasts.append({"variant": variant, "comparator": "T315A",
                          "expected_direction": "lower binding confidence for higher measured net growth",
                          "binding_confidence_lower_than_control": binding < control if binding is not None and control is not None else None})
    return {"rows": rows, "prespecified_contrasts": contrasts,
            "n_evaluable_contrasts": sum(c["binding_confidence_lower_than_control"] is not None for c in contrasts),
            "n_supported_directions": sum(c["binding_confidence_lower_than_control"] is True for c in contrasts),
            "caveats": ["Three selected mutants; no general accuracy estimate or calibrated work/fail threshold.",
                        "Wild-type assay response is unavailable in the substitution CSV; no invented baseline measurement.",
                        "Binding confidence is not affinity, drug efficacy, or a calibrated escape probability.",
                        "Five structural samples do not provide replicate binding-response uncertainty.",
                        "Net growth combines protein function and response to imatinib."]}


def markdown(summary: dict) -> str:
    def number(value, places=5):
        return "unavailable" if value is None else f"{value:.{places}f}"
    lines = ["# Boltz ABL1/imatinib contrast benchmark", "",
             "271-aa crop (242-512), five samples, boltz-2.1, identical unconstrained settings.", "",
             "| Variant | Measured net growth | Binding confidence | Delta vs WT | Mean pLDDT +/- sample SD |",
             "| --- | --- | --- | --- | --- |"]
    for row in summary["rows"]:
        lines.append(f"| {row['variant']} | {number(row.get('measured_netgrowth'))} | {number(row.get('binding_confidence'))} | {number(row.get('delta_binding_confidence_vs_wt'))} | {number(row.get('confidence'), 2)} +/- {number(row.get('confidence_sd'), 2)} |")
    lines += ["", "## Prespecified contrasts", ""]
    for contrast in summary["prespecified_contrasts"]:
        lines.append(f"- {contrast['variant']} lower binding confidence than T315A: {contrast['binding_confidence_lower_than_control']}")
    lines += ["", "## Limits", ""] + [f"- {item}" for item in summary["caveats"]]
    return "\n".join(lines) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--max-total-usd", type=float, default=0.50)
    args = ap.parse_args()
    jobs = plan()
    OUT.mkdir(parents=True, exist_ok=True)
    estimates = []
    for job in jobs:
        estimate = bs.estimate_cost(job["body"])
        estimates.append({"variant": job["variant"], "estimate_usd": float(estimate["estimated_cost_usd"]),
                          "idempotency_key": job["idempotency_key"]})
    total = sum(row["estimate_usd"] for row in estimates)
    manifest = {"model": bs.MODEL, "crop": [242, 512], "num_samples": 5,
                "compound": "imatinib", "constraints": "none", "estimated_total_usd": total,
                "reference_sha256": hashlib.sha256(REFERENCE.read_bytes()).hexdigest(),
                "dataset_sha256": hashlib.sha256(dms.NETGROWTH_CSV.read_bytes()).hexdigest(), "jobs": estimates}
    (OUT / "estimate.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2), flush=True)
    if total > args.max_total_usd:
        raise ValueError(f"Estimate ${total:.2f} exceeds total cap ${args.max_total_usd:.2f}")
    if not args.run:
        print("Estimate only; no prediction submitted.")
        return 0
    # Never relabel old results with a new estimate/source manifest.
    (OUT / "summary.json").unlink(missing_ok=True)
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    results = {}
    for job in jobs:
        variant = job["variant"]
        print(f"Running {variant}...", flush=True)
        try:
            st = bs.predict_structure(job["sequence"], name=f"benchmark_{variant}",
                                      out_dir=OUT / "runs", compound="imatinib", num_samples=5,
                                      max_usd=min(0.25, next(e["estimate_usd"] for e in estimates if e["variant"] == variant)))
            row = asdict(st)
            row.pop("sequence")
            row["path"] = str(st.path)
            raw_path = next(parent / "run.json" for parent in st.path.parents if (parent / "run.json").is_file())
            raw = json.loads(raw_path.read_text())
            row["binding_metrics"] = raw.get("output", {}).get("binding_metrics")
            row["sample_metrics"] = [s.get("metrics", {}) for s in raw.get("output", {}).get("all_sample_results", [])]
            row.update({k: job[k] for k in ("measured_netgrowth", "replicate_1", "replicate_2", "idempotency_key")})
            results[variant] = row
            print(f"{variant}: binding={st.binding_confidence}, pLDDT={st.confidence:.2f}", flush=True)
        except (bs.BoltzUnavailable, bs.BoltzTooExpensive, OSError, ValueError, subprocess.SubprocessError) as exc:
            results[variant] = {"error": str(exc)}
            print(f"{variant}: failed: {exc}", flush=True)
        (OUT / "results.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    summary = summarize(results)
    summary["manifest_sha256"] = hashlib.sha256((OUT / "manifest.json").read_bytes()).hexdigest()
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (OUT / "report.md").write_text(markdown(summary), encoding="utf-8")
    print(markdown(summary))
    return int(any("error" in row for row in results.values()))


if __name__ == "__main__":
    raise SystemExit(main())
