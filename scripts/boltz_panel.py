#!/usr/bin/env python3
"""Boltz-2 structure + binding across a 20-variant BCR-ABL panel, plus wild-type.

    python3 scripts/boltz_panel.py                 # estimate only, spends nothing
    python3 scripts/boltz_panel.py --run           # submit (asks for the total first)

WHY A WILD-TYPE RUN IS IN THE PANEL. Ptarmigan-1 cannot score mutants at all,
so the original mutant-minus-wildtype affinity screen is unsupported there.
Boltz accepts both mutant and wild-type complexes. Comparing their binding
confidence is possible, but this is not a measured affinity shift or a calibrated
efficacy prediction. Keep wild type as a reference for descriptive comparisons.

HOW THE PANEL WAS CHOSEN. Two requirements, both load-bearing:

1. Every variant is MEASURED in abl_dms_netgrowth_1200nM.csv, so each
   prediction has ground truth to be scored against rather than admired.
2. The panel spans the full measured range (-0.062 to +0.049), not just the
   famous resistant end. A panel of 20 resistance mutations would all sit at
   one extreme and could not show discrimination -- every prediction would
   look good for the wrong reason. So it carries 14 resistant, 3 near-neutral
   and 3 kinase-dead variants.

It deliberately includes all five variants that the 480-dim mutant-marginal
embeddings got WORST (G250S, G250A, G250E, L248C, L273N). Those are measured at
about +0.04 and were predicted at about -0.023, i.e. the population mean: the
embeddings put the most resistant variants in the middle of the distribution
because local sequence context does not know imatinib stopped binding. Pocket
geometry is the obvious thing that might. Whether Boltz recovers them is the
most informative result this panel can produce, so they are all in it.

Cost: 21 runs x $0.10 = $2.10 at the kinase crop with 5 samples. Each call
re-runs the free estimate first and refuses to exceed the per-call cap.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pipeline import boltz_structure as bs

DMS = Path("data/datasets/abl_dms_netgrowth_1200nM.csv")
SEQ = Path("data/cancer-exp/ABL1_NM_005157.6_DMS_reference_aa.txt")
OUT = Path("output/structures/panel_results.json")

# ATP-site landmarks in P00519 numbering. Wild-type structural annotation,
# identical for every variant, so it biases the pose without telling the model
# anything about this variant's resistance.
POCKET = [271, 286, 315, 318, 381, 393]

PANEL: list[tuple[str, str]] = [
    # --- resistant: the canonical imatinib-escape set -------------------
    ("T315I", "gatekeeper; the classic escape"),
    ("Y253H", "P-loop"),
    ("E255K", "P-loop"),
    ("E255V", "P-loop"),
    ("Q252H", "P-loop"),
    ("M244V", "P-loop edge"),
    ("F317L", "pocket contact"),
    ("H396P", "activation loop"),
    ("F486S", "C-terminal lobe"),
    # --- the five the embeddings got worst -----------------------------
    ("G250E", "P-loop; embeddings put it at the population mean"),
    ("G250S", "P-loop; embeddings put it at the population mean"),
    ("G250A", "P-loop; embeddings put it at the population mean"),
    ("L248C", "P-loop; embeddings put it at the population mean"),
    ("L273N", "embeddings put it at the population mean"),
    # --- near-neutral controls -----------------------------------------
    ("T315A", "same gatekeeper residue as T315I but NOT resistant"),
    ("V299L", "near-neutral"),
    ("L327A", "near-neutral"),
    # --- kinase-dead / sensitive controls ------------------------------
    ("L471Q", "kinase-dead end of the range"),
    ("P484Q", "kinase-dead end of the range"),
    ("H490T", "most negative measured variant"),
]


def wildtype() -> str:
    return "".join(
        x for x in SEQ.read_text().split() if not x.startswith(">")
    )


def measured() -> dict[str, float]:
    out: dict[str, float] = {}
    for r in csv.DictReader(DMS.open()):
        if not r["protein_start"]:
            continue
        key = f"{r['ref_aa']}{r['protein_start']}{r['alt_aa']}"
        raw = r.get("netgrowth_mean_corrected") or r.get("netgrowth_mean")
        try:
            out[key] = float(raw)
        except (TypeError, ValueError):
            continue
    return out


def apply(seq: str, variant: str) -> str:
    """Apply a single substitution in P00519 numbering, validating the ref."""
    ref, pos, alt = variant[0], int(variant[1:-1]), variant[-1]
    if seq[pos - 1] != ref:
        raise ValueError(
            f"{variant}: P00519 position {pos} holds {seq[pos - 1]!r}, "
            f"not {ref!r}. Do not adjust the offset to make this pass."
        )
    return seq[:pos - 1] + alt + seq[pos:]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="store_true",
                    help="actually submit; without this, estimate only")
    ap.add_argument("--yes", action="store_true",
                    help="skip the interactive confirmation of the total")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--compound", default="imatinib")
    args = ap.parse_args()

    if not bs.available():
        print("boltz-api is not authenticated; run `boltz-api auth login`")
        return 1

    wt = wildtype()
    meas = measured()
    panel = PANEL[:args.limit] if args.limit else PANEL

    # Validate every variant BEFORE spending anything -- the same discipline
    # embed_variants.py uses, for the same reason.
    jobs: list[tuple[str, str, str]] = [("WILDTYPE", wt, "reference for the delta")]
    for variant, note in panel:
        try:
            jobs.append((variant, apply(wt, variant), note))
        except (ValueError, IndexError) as exc:
            print(f"REFUSING to run: {exc}")
            return 1
        if variant not in meas:
            print(f"REFUSING to run: {variant} is not measured in the DMS set")
            return 1

    # Free preflight for the whole panel.
    total = 0.0
    for name, seq, _ in jobs:
        st = bs.predict_structure(
            seq, name=f"ABL1_{name}", pocket_hint=POCKET,
            compound=args.compound, dry_run=True,
        )
        total += st.cost_usd or 0.0
    print(f"panel: {len(jobs)} runs ({len(jobs) - 1} variants + wild-type)")
    print(f"crop {bs.CROP_START}-{bs.CROP_END}, {bs.NUM_SAMPLES} samples, "
          f"binding on, model {bs.MODEL}")
    print(f"ESTIMATED TOTAL: ${total:.2f}")

    if not args.run:
        print("\nestimate only; nothing submitted. Re-run with --run to spend.")
        return 0

    if not args.yes:
        reply = input(f"\nSubmit {len(jobs)} runs for ~${total:.2f}? [y/N] ")
        if reply.strip().lower() not in {"y", "yes"}:
            print("aborted; nothing submitted.")
            return 0

    results: dict[str, dict] = {}
    if OUT.exists():
        results = json.loads(OUT.read_text())

    spent = 0.0
    for i, (name, seq, note) in enumerate(jobs, 1):
        if name in results and "error" not in results[name]:
            print(f"[{i}/{len(jobs)}] {name}: cached, skipping")
            continue
        print(f"[{i}/{len(jobs)}] {name}…", end=" ", flush=True)
        try:
            st = bs.predict_structure(
                seq, name=f"ABL1_{name}", pocket_hint=POCKET,
                compound=args.compound,
            )
        except (bs.BoltzUnavailable, bs.BoltzTooExpensive) as exc:
            print(f"FAILED: {exc}")
            # A run that completed but failed to PARSE was still billed, so
            # count it. Reporting $0.00 after a parse failure understates real
            # spend, which is the one number that must never be optimistic.
            billed = bs.recorded_cost(f"ABL1_{name}")
            spent += billed
            results[name] = {"error": str(exc), "note": note,
                             "billed_usd": billed or None}
            OUT.parent.mkdir(parents=True, exist_ok=True)
            OUT.write_text(json.dumps(results, indent=2))
            continue
        spent += st.cost_usd or 0.0
        results[name] = {
            "note": note,
            "measured_netgrowth": meas.get(name),
            "confidence": st.confidence,
            "confidence_sd": st.confidence_sd,
            "binding_confidence": st.binding_confidence,
            "ptm": st.ptm,
            "iptm": st.iptm,
            "crop": list(st.crop) if st.crop else None,
            "num_samples": st.num_samples,
            "cost_usd": st.cost_usd,
            "prediction_id": st.prediction_id,
            "path": str(st.path),
        }
        OUT.parent.mkdir(parents=True, exist_ok=True)
        OUT.write_text(json.dumps(results, indent=2))
        print(f"pLDDT={st.confidence:.1f}"
              + (f"±{st.confidence_sd:.1f}" if st.confidence_sd else "")
              + (f"  binding={st.binding_confidence:.3f}"
                 if st.binding_confidence is not None else ""))

    print(f"\nspent this run: ${spent:.2f}  ->  {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
