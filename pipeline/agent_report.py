"""Deterministic cancer research reports, with an optional unimplemented agent hook."""

from __future__ import annotations

import json
import os
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

REPORT_DIR = Path("output/reports")


@dataclass
class Report:
    path: Path
    markdown: str
    generated_by: str


def _slug(*parts: str) -> str:
    return "_".join(p.lower().replace(" ", "-") for p in parts if p)


def render_markdown(payload: dict[str, Any]) -> str:
    """Deterministic fallback report — no LLM, always works."""
    gene = payload.get("gene", "?")
    variant = payload.get("variant", "?")
    mode = payload.get("mode", "?")
    screen = payload.get("screen") or {}
    results = screen.get("results") or []
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    lines = [
        f"# {gene} {variant} — variant response research report",
        "",
        f"*Mode:* {mode}  ·  *Generated:* {stamp}",
        "",
        "## Variant",
        f"- Gene / protein: **{gene}** ({payload.get('uniprot', 'n/a')})",
        f"- Organism: {payload.get('organism', 'n/a')}",
        f"- Variant: **{variant}**",
        f"- In annotated binding site: {payload.get('in_binding_site', 'unknown')}",
    ]

    esm = payload.get("esm2")
    if esm:
        lines += [
            "",
            "## Sequence-model evidence (ESM2)",
            f"- Δ log-likelihood: **{esm.get('delta_ll'):+.3f}**",
            f"- Model: {esm.get('model', 'n/a')}",
        ]

    pred = payload.get("prediction") or {}
    if pred:
        val = pred.get("predicted_netgrowth", pred.get("measured_netgrowth"))
        lines += ["", "## Net growth under imatinib (1200 nM)"]
        if pred.get("predicted"):
            lines += [
                f"- Predicted net growth: **{val:+.5f}**",
                f"- Percentile vs. the 4,925-variant screen: "
                f"{pred.get('percentile_vs_screen')}%",
                f"- Model: {pred.get('source')}",
                "",
                f"*Validation:* {pred.get('validation_status', 'unavailable')}. "
                + pred.get("caveat", "Research use only."),
            ]
        else:
            lines += [
                f"- **Measured** net growth: **{val:+.5f}**",
                f"- Replicate spread: {pred.get('replicate_spread')}",
                f"- Source: {pred.get('source')}",
                "",
                "*This variant is in the screen, so a real measurement is reported "
                "instead of a prediction.*",
            ]

    esm_prediction = payload.get("esm_prediction")
    if esm_prediction:
        lines += ["", "## Persisted ESM model estimate",
                  f"- Estimate: **{esm_prediction['predicted_netgrowth']:+.5f}**",
                  f"- Head: {esm_prediction['head']}; model: {esm_prediction['model']}",
                  f"- Position membership: {esm_prediction['split_membership']}",
                  f"- Model SHA256: {esm_prediction['model_sha256']}",
                  f"- Caveat: {esm_prediction['caveat']}"]
    validation = payload.get("validation")
    if validation:
        lines += ["", "## Saved held-out model estimate",
                  f"- Estimate: **{validation['predicted']:+.5f}**",
                  f"- Model/head: {validation['model']} / {validation['head']}",
                  f"- Split: {validation['split']}; seed {validation['seed']}",
                  f"- Artifact SHA256: {validation['artifact_sha256']}",
                  f"- Provenance limit: {validation['provenance_limit']}"]
    structure = payload.get("structure")
    if structure:
        lines += ["", "## Boltz structure (ABL1/imatinib)"]
        for key in ("confidence", "confidence_sd", "binding_confidence", "ptm", "iptm", "crop", "num_samples", "cost_usd", "prediction_id", "path"):
            lines.append(f"- {key}: {structure.get(key)}")
        lines.append("Binding confidence is not affinity or drug efficacy.")
    skipped = payload.get("skipped_stages", [])
    if skipped:
        lines += ["", "## Skipped or unsupported stages"]
        lines.extend(f"- {item['stage']}: {item['reason']}" for item in skipped)

    base = payload.get("baseline") or {}
    cells = [c for c in (base.get("cells") or []) if c.get("call") != "no_reference"]
    if cells:
        cutoff = base.get("contact_cutoff_angstrom", 5.0)
        lines += [
            "",
            "## Pocket-distance baseline",
            "",
            "Minimum heavy-atom distance from the variant site to the bound",
            f"inhibitor in a public reference structure. Under {cutoff:g} Å the site",
            "touches the drug.",
            "",
            "| Compound | Position | PDB | Distance (Å) | Contact score | Call |",
            "| --- | --- | --- | --- | --- | --- |",
        ]
        for c in cells:
            approx = "~" if c.get("approximate") else ""
            lines.append(
                f"| {c.get('compound')} | {c.get('published_position')} "
                f"| {c.get('pdb')} | {c.get('distance_angstrom'):.2f}{approx} "
                f"| {c.get('contact_score'):.2f} "
                f"| {c.get('call', 'n/a').replace('_', ' ')} |"
            )
        if any(c.get("approximate") for c in cells):
            lines += [
                "",
                "~ = the reference structure holds a related, not identical, "
                "inhibitor, so the distance is approximate.",
            ]
        if not base.get("calibrated", False):
            lines += [
                "",
                "**Uncalibrated, and measurably weak.** This is a geometric prior, "
                "not an efficacy prediction. Against 4,665 measured BCR-ABL "
                "variants, distance alone scores Pearson r = 0.069 — under 0.5% "
                "of the variance. It is blind to any mechanism that does not "
                "touch the ligand: H396 escapes imatinib allosterically from "
                "15.5 Å away. Use it to triage, not to decide.",
            ]

    if results:
        lines += [
            "",
            "## Drug panel",
            "",
            "| Drug | Class | Δ affinity | Verdict |",
            "| --- | --- | --- | --- |",
        ]
        for r in results:
            d = r.get("delta_affinity")
            lines.append(
                f"| {r.get('name')} | {(r.get('metadata') or {}).get('class', '')} "
                f"| {f'{d:+.2f}' if isinstance(d, (int, float)) else 'n/a'} "
                f"| {r.get('verdict', 'n/a')} |"
            )
        kept = [r["name"] for r in results if r.get("verdict") == "work"]
        lines += ["", f"**Predicted to retain activity:** {', '.join(kept) or 'none'}"]
    elif cells:
        lines += [
            "",
            "## Drug panel",
            "",
            "_Ptarmigan cannot screen mutant sequences. The "
            "pocket-distance baseline above is a geometric descriptor only._",
        ]
    else:
        lines += ["", "## Drug panel", "", "_No screen results available._"]

    lines += [
        "",
        "## Caveats",
        "- Computational prediction only; not a clinical recommendation.",
        "- Verdict thresholds are uncalibrated until the benchmark controls pass.",
        "",
    ]
    return "\n".join(lines)


def write_report(payload: dict[str, Any], *, use_agent: bool = True) -> Report:
    """Write a report to output/reports/, via Keenable if available."""
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    name = _slug(payload.get("gene", "gene"), payload.get("variant", "variant"))
    stamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S-%f")

    markdown, generated_by = None, "fallback-renderer"
    if use_agent and available():
        try:
            markdown = _call_keenable(payload)
            generated_by = "keenable-agent"
        except (NotImplementedError, RuntimeError) as exc:
            print(f"  [agent_report] falling back to local renderer: {exc}")
    if markdown is None:
        markdown = render_markdown(payload)

    path = REPORT_DIR / f"{name}_{stamp}.md"
    path.write_text(markdown, encoding="utf-8")
    (REPORT_DIR / f"{name}_{stamp}.json").write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    return Report(path=path, markdown=markdown, generated_by=generated_by)


def _call_keenable(payload: dict[str, Any]) -> str:
    key = os.getenv("KEENABLE_API_KEY")
    if not key:
        raise RuntimeError("KEENABLE_API_KEY not set (see .env.example)")
    raise NotImplementedError("Keenable request not implemented: API docs still needed.")


def available() -> bool:
    return bool(os.getenv("KEENABLE_API_KEY"))
