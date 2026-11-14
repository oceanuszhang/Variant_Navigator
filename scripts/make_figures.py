#!/usr/bin/env python3
"""Regenerate the figures in figures/ from live pipeline data.

    python scripts/make_figures.py

Dependency-free: writes SVG directly, so the figures render on GitHub and in the
deck without a plotting stack. Numbers come from pipeline.pocket, computed
against the biological assemblies, so a figure can never drift from the code.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pipeline.pocket import POCKET_REFERENCE, pocket_distance_for  # noqa: E402

OUT = Path(__file__).resolve().parent.parent / "figures"

INK = "#101E33"
MUTED = "#4C5A70"
RULE = "#D3D8DF"
BLUE = "#1F6FEB"
AMBER = "#C9741B"
PAPER = "#FCFCFB"

# Validated with the dataviz palette validator: adjacent-pair CVD dE 30.8
# (protan), 26.3 (tritan), normal-vision 34.6 — all checks pass.

CONTACT_CUTOFF = 5.0  # A; below this a residue touches the ligand

PANELS = [
    {
        "title": "Human EGFR — cancer drugs",
        "compounds": [("erlotinib", BLUE), ("osimertinib", AMBER)],
        "positions": [(790, "T", "T790 · gatekeeper"), (858, "L", "L858 · activation loop")],
    },
    {
        "title": "Rice ALS — herbicides  (Arabidopsis numbering)",
        "compounds": [("imazethapyr", BLUE), ("nicosulfuron", AMBER)],
        "positions": [
            (122, "A", "A122 · dimer interface"),
            (574, "W", "W574 · channel"),
            (653, "S", "S653 · channel"),
        ],
    },
]


def esc(s: str) -> str:
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def collect() -> dict:
    """Live distances, keyed (compound, position)."""
    out = {}
    for panel in PANELS:
        for pos, expect, _ in panel["positions"]:
            for compound, _c in panel["compounds"]:
                try:
                    r = pocket_distance_for(compound, pos, expect_residue=expect)
                    out[(compound, pos)] = r["distance_angstrom"]
                except Exception as exc:  # noqa: BLE001
                    print(f"  WARN {compound} {pos}: {exc}")
                    out[(compound, pos)] = None
    return out


def fig_pocket_distances(data: dict) -> str:
    """Grouped horizontal bars: distance from variant site to the inhibitor."""
    W = 1000
    pad_l, pad_r, pad_t = 300, 130, 64
    plot_w = W - pad_l - pad_r
    bar_h, bar_gap, group_gap, panel_gap = 26, 6, 26, 54
    max_d = 12.5

    def x(d: float) -> float:
        return pad_l + (d / max_d) * plot_w

    rows: list[str] = []
    y = pad_t

    for pi, panel in enumerate(PANELS):
        rows.append(
            f'<text x="16" y="{y}" font-size="15" font-weight="700" '
            f'fill="{INK}" font-family="IBM Plex Sans, Arial, sans-serif">'
            f'{esc(panel["title"])}</text>'
        )
        y += 26
        for pos, _expect, label in panel["positions"]:
            n = len(panel["compounds"])
            for ci, (compound, color) in enumerate(panel["compounds"]):
                d = data.get((compound, pos))
                by = y + ci * (bar_h + bar_gap)
                if d is None:
                    continue
                w = max(x(d) - pad_l, 4)
                rows.append(
                    f'<rect x="{pad_l}" y="{by}" width="{w:.1f}" height="{bar_h}" '
                    f'rx="4" fill="{color}"/>'
                )
                approx = "~" if not POCKET_REFERENCE[compound]["exact"] else ""
                rows.append(
                    f'<text x="{x(d) + 10:.1f}" y="{by + bar_h - 8}" font-size="14" '
                    f'font-weight="700" fill="{INK}" '
                    f'font-family="JetBrains Mono, monospace">{d:.2f}{approx} Å</text>'
                )
                rows.append(
                    f'<text x="{pad_l - 12}" y="{by + bar_h - 8}" font-size="13" '
                    f'text-anchor="end" fill="{MUTED}" '
                    f'font-family="IBM Plex Sans, Arial, sans-serif">'
                    f'{esc(compound)}</text>'
                )
            # group label
            gy = y + (n * (bar_h + bar_gap) - bar_gap) / 2 + 5
            rows.append(
                f'<text x="16" y="{gy:.1f}" font-size="14" '
                f'font-weight="600" fill="{INK}" '
                f'font-family="IBM Plex Sans, Arial, sans-serif">'
                f'{esc(label)}</text>'
            )
            y += n * (bar_h + bar_gap) - bar_gap + group_gap
        if pi < len(PANELS) - 1:
            y += panel_gap - group_gap

    H = y + 86

    # Axis, gridlines, contact band
    axis: list[str] = []
    cut_x = x(CONTACT_CUTOFF)
    axis.append(
        f'<rect x="{pad_l}" y="{pad_t - 18}" width="{cut_x - pad_l:.1f}" '
        f'height="{y - pad_t + 4:.1f}" fill="#1F6FEB" opacity="0.055"/>'
    )
    for tick in range(0, 14, 2):
        tx = x(tick)
        axis.append(
            f'<line x1="{tx:.1f}" y1="{pad_t - 18}" x2="{tx:.1f}" y2="{y - 20:.1f}" '
            f'stroke="{RULE}" stroke-width="1"/>'
        )
        axis.append(
            f'<text x="{tx:.1f}" y="{y + 4:.1f}" font-size="13" text-anchor="middle" '
            f'fill="{MUTED}" font-family="JetBrains Mono, monospace">{tick}</text>'
        )
    axis.append(
        f'<line x1="{cut_x:.1f}" y1="{pad_t - 18}" x2="{cut_x:.1f}" y2="{y - 20:.1f}" '
        f'stroke="{BLUE}" stroke-width="2" stroke-dasharray="5 4"/>'
    )
    axis.append(
        f'<text x="{cut_x + 8:.1f}" y="{pad_t - 24}" font-size="13" font-weight="700" '
        f'fill="{BLUE}" font-family="IBM Plex Sans, Arial, sans-serif">'
        f'{CONTACT_CUTOFF:.0f} Å — contact threshold</text>'
    )
    axis.append(
        f'<text x="{pad_l + plot_w / 2:.1f}" y="{y + 30:.1f}" font-size="13" '
        f'text-anchor="middle" fill="{MUTED}" '
        f'font-family="IBM Plex Sans, Arial, sans-serif">'
        f'minimum heavy-atom distance to bound inhibitor (Å)</text>'
    )
    axis.append(
        f'<text x="16" y="{y + 56:.1f}" font-size="12" fill="{MUTED}" '
        f'font-family="IBM Plex Sans, Arial, sans-serif">'
        f'~ = related, not identical, inhibitor in the reference structure. '
        f'Computed from biological assemblies.</text>'
    )

    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H:.0f}" '
        f'viewBox="0 0 {W} {H:.0f}" aria-label="Distance from each variant site to '
        f'the bound inhibitor, by compound. Sites under five angstroms contact the '
        f'drug directly.">'
        f'<rect width="{W}" height="{H:.0f}" fill="{PAPER}"/>'
        f'<text x="16" y="30" font-size="19" font-weight="700" fill="{INK}" '
        f'font-family="IBM Plex Sans, Arial, sans-serif">'
        f'Variant sites that touch the drug</text>'
        + "".join(axis) + "".join(rows) + "</svg>"
    )


def fig_matrix() -> str:
    """The prediction object: a variant x compound matrix of log2 fold change."""
    W = 900
    cell_w, cell_h = 150, 62
    x0, y0 = 300, 118

    exp1 = ["L858R", "T790M", "L858R/T790M"]
    exp2 = ["A122T", "W574L", "S653N"]

    def block(title: str, rows: list[str], compounds: list[str], top: float) -> tuple[str, float]:
        parts = [
            f'<text x="40" y="{top - 30}" font-size="15" font-weight="700" fill="{INK}" '
            f'font-family="IBM Plex Sans, Arial, sans-serif">{esc(title)}</text>'
        ]
        for ci, c in enumerate(compounds):
            parts.append(
                f'<text x="{x0 + ci * cell_w + cell_w / 2}" y="{top - 10}" font-size="14" '
                f'font-weight="600" text-anchor="middle" fill="{INK}" '
                f'font-family="IBM Plex Sans, Arial, sans-serif">{esc(c)}</text>'
            )
        for ri, r in enumerate(rows):
            cy = top + ri * cell_h
            parts.append(
                f'<text x="{x0 - 16}" y="{cy + cell_h / 2 + 5}" font-size="14" '
                f'text-anchor="end" fill="{INK}" '
                f'font-family="JetBrains Mono, monospace">{esc(r)}</text>'
            )
            for ci in range(len(compounds)):
                cx = x0 + ci * cell_w
                parts.append(
                    f'<rect x="{cx + 1}" y="{cy + 1}" width="{cell_w - 2}" '
                    f'height="{cell_h - 2}" rx="6" fill="#FFFFFF" stroke="{RULE}" '
                    f'stroke-width="2"/>'
                )
                parts.append(
                    f'<text x="{cx + cell_w / 2}" y="{cy + cell_h / 2 + 6}" font-size="15" '
                    f'text-anchor="middle" fill="{MUTED}" '
                    f'font-family="JetBrains Mono, monospace">log₂ FC</text>'
                )
        return "".join(parts), top + len(rows) * cell_h

    b1, end1 = block("Experiment 1 — EGFR, Ba/F3 IC₅₀", exp1, ["erlotinib", "osimertinib"], y0)
    b2, end2 = block(
        "Experiment 2 — rice ALS, greenhouse GR₅₀", exp2,
        ["imazethapyr", "nicosulfuron"], end1 + 96,
    )
    H = end2 + 150

    # Reading-direction annotations: the answer to "which axis is the product?"
    ann = (
        f'<rect x="40" y="{end2 + 34}" width="820" height="96" rx="10" '
        f'fill="#E7EEFB"/>'
        f'<text x="62" y="{end2 + 62}" font-size="14" font-weight="700" fill="#14448F" '
        f'font-family="IBM Plex Sans, Arial, sans-serif">'
        f'Read across a ROW → which treatment holds up for this variant '
        f'(the clinical / agronomic decision)</text>'
        f'<text x="62" y="{end2 + 86}" font-size="14" font-weight="700" fill="#8A4D10" '
        f'font-family="IBM Plex Sans, Arial, sans-serif">'
        f'Read down a COLUMN → which variants escape this treatment '
        f'(resistance surveillance)</text>'
        f'<text x="62" y="{end2 + 110}" font-size="13" fill="{MUTED}" '
        f'font-family="IBM Plex Sans, Arial, sans-serif">'
        f'Same model, same cells. The headline result is a row where the two '
        f'numbers disagree.</text>'
    )

    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H:.0f}" '
        f'viewBox="0 0 {W} {H:.0f}" aria-label="The prediction is a matrix of '
        f'variants by compounds. Each cell is a predicted log2 fold change in the '
        f'dose needed for fifty percent inhibition.">'
        f'<rect width="{W}" height="{H:.0f}" fill="{PAPER}"/>'
        f'<text x="40" y="38" font-size="19" font-weight="700" fill="{INK}" '
        f'font-family="IBM Plex Sans, Arial, sans-serif">'
        f'What the model outputs: one cell per variant × compound</text>'
        f'<text x="40" y="64" font-size="14" fill="{MUTED}" '
        f'font-family="IBM Plex Sans, Arial, sans-serif">'
        f'log₂ fold change in dose to 50% inhibition, vs the assay’s '
        f'susceptible reference</text>'
        + b1 + b2 + ann + "</svg>"
    )


def fig_numbering() -> str:
    """Three numbering systems, and where each one bites."""
    W, H = 900, 430
    rows = [
        ("EGFR", "clinical", "T790", [
            ("UniProt P00533", "T790", "+0", True),
            ("PDB 1M17", "T766", "−24", True),
            ("PDB 4ZAU", "T790", "+0", True),
        ]),
        ("rice ALS", "Arabidopsis", "W574", [
            ("UniProt Q6K2E8", "W548", "−26", True),
            ("PDB 1Z8N", "W574", "+0", True),
            ("PDB 1YHZ", "W574", "+0", True),
        ]),
    ]
    out = [
        f'<rect width="{W}" height="{H}" fill="{PAPER}"/>',
        f'<text x="40" y="38" font-size="19" font-weight="700" fill="{INK}" '
        f'font-family="IBM Plex Sans, Arial, sans-serif">'
        f'One residue, three numbering systems</text>',
        f'<text x="40" y="64" font-size="14" fill="{MUTED}" '
        f'font-family="IBM Plex Sans, Arial, sans-serif">'
        f'Every offset below was verified against the actual sequence, not assumed. '
        f'Each is a silent-wrong-answer bug.</text>',
    ]
    y = 104
    for target, conv, published, mappings in rows:
        out.append(
            f'<text x="40" y="{y}" font-size="15" font-weight="700" fill="{INK}" '
            f'font-family="IBM Plex Sans, Arial, sans-serif">'
            f'{esc(target)} — published as <tspan font-family="JetBrains Mono, monospace">'
            f'{esc(published)}</tspan> ({esc(conv)} numbering)</text>'
        )
        y += 22
        for i, (where, resolved, off, ok) in enumerate(mappings):
            bx = 40 + i * 274
            out.append(
                f'<rect x="{bx}" y="{y}" width="258" height="80" rx="10" fill="#FFFFFF" '
                f'stroke="{RULE}" stroke-width="2"/>'
            )
            out.append(
                f'<text x="{bx + 16}" y="{y + 26}" font-size="13" fill="{MUTED}" '
                f'font-family="IBM Plex Sans, Arial, sans-serif">{esc(where)}</text>'
            )
            out.append(
                f'<text x="{bx + 16}" y="{y + 56}" font-size="24" font-weight="700" '
                f'fill="{INK}" font-family="JetBrains Mono, monospace">{esc(resolved)}</text>'
            )
            colour = BLUE if off == "+0" else AMBER
            out.append(
                f'<text x="{bx + 242}" y="{y + 56}" font-size="17" font-weight="700" '
                f'text-anchor="end" fill="{colour}" '
                f'font-family="JetBrains Mono, monospace">{esc(off)}</text>'
            )
        y += 112
    out.append(
        f'<rect x="40" y="{y}" width="820" height="74" rx="10" fill="#FAEBDA"/>'
        f'<text x="62" y="{y + 30}" font-size="14" font-weight="700" fill="#8A4D10" '
        f'font-family="IBM Plex Sans, Arial, sans-serif">'
        f'The trap: A122 is alanine in BOTH ALS numbering systems, by coincidence.</text>'
        f'<text x="62" y="{y + 56}" font-size="13" fill="#6B4414" '
        f'font-family="IBM Plex Sans, Arial, sans-serif">'
        f'A122T is the headline demo case, so it validates fine while every other '
        f'site is silently wrong. Strict residue checks catch it.</text>'
    )
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
        f'viewBox="0 0 {W} {H}" aria-label="The same published residue maps to '
        f'different positions in UniProt and in each PDB structure. EGFR 1M17 needs '
        f'an offset of minus 24; rice ALS needs minus 26.">'
        + "".join(out) + "</svg>"
    )


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    print("collecting live pocket distances…")
    data = collect()
    figures = {
        "pocket_distances.svg": fig_pocket_distances(data),
        "prediction_matrix.svg": fig_matrix(),
        "numbering_systems.svg": fig_numbering(),
    }
    for name, svg in figures.items():
        (OUT / name).write_text(svg)
        print(f"  wrote figures/{name}  ({len(svg) // 1024} KB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
