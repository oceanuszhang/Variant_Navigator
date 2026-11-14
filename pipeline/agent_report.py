"""Keenable agent: synthesise pipeline outputs into a clinician/agronomist report.

STUB — Keenable API docs were not available in this working copy. The markdown
fallback writer below is real, so the pipeline produces a readable report even
with no agent credentials.
"""

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
        f"# {gene} {variant} — drug susceptibility report",
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
    else:
        lines += ["", "## Drug panel", "", "_No screen results available._"]

    lines += [
        "",
        "## Caveats",
        "- Computational prediction only; not a clinical or agronomic recommendation.",
        "- Verdict thresholds are uncalibrated until the benchmark controls pass.",
        "",
    ]
    return "\n".join(lines)


def write_report(payload: dict[str, Any], *, use_agent: bool = True) -> Report:
    """Write a report to output/reports/, via Keenable if available."""
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    name = _slug(payload.get("gene", "gene"), payload.get("variant", "variant"))
    stamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S")

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
    path.write_text(markdown)
    (REPORT_DIR / f"{name}_{stamp}.json").write_text(json.dumps(payload, indent=2, default=str))
    return Report(path=path, markdown=markdown, generated_by=generated_by)


def _call_keenable(payload: dict[str, Any]) -> str:
    key = os.getenv("KEENABLE_API_KEY")
    if not key:
        raise RuntimeError("KEENABLE_API_KEY not set (see .env.example)")
    raise NotImplementedError("Keenable request not implemented: API docs still needed.")


def available() -> bool:
    return bool(os.getenv("KEENABLE_API_KEY"))
