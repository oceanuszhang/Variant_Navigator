#!/usr/bin/env python3
"""Predicted vs measured net growth — the demo scatter.

    .venv/bin/python scripts/plot_predictions.py

Reads output/model/validation_predictions.json (written by
scripts/train_embeddings.py) and writes figures/predicted_vs_measured.svg.

Dependency-free SVG, like the rest of figures/, so it renders on GitHub and in
a deck with no plotting stack.

Design notes, from the dataviz procedure:
  - form: scatter, because the job is agreement between two continuous
    quantities. The y=x diagonal is the reference the eye actually needs.
  - the replicate band is drawn as a shaded ribbon around the diagonal: it is
    the width within which a prediction is as good as redoing the experiment,
    and the honest target to be judged against.
  - two categories (ordinary / highlighted), so a legend is present; the
    highlighted few are also directly labelled, so identity is never carried
    by colour alone.
  - palette #1F6FEB / #C9741B validated: adjacent-pair CVD dE 30.8 protan,
    26.3 tritan, 34.6 normal vision, contrast >= 3:1 on this surface. All pass.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

SRC = Path("output/model/validation_predictions.json")
OUT = Path("figures/predicted_vs_measured.svg")

INK = "#101E33"
MUTED = "#4C5A70"
RULE = "#D3D8DF"
BLUE = "#1F6FEB"
AMBER = "#C9741B"
PAPER = "#FCFCFB"


def esc(s: str) -> str:
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def main() -> int:
    if not SRC.exists():
        print(f"missing {SRC} — run scripts/train_embeddings.py first")
        return 1
    d = json.loads(SRC.read_text())
    rows = d["rows"]
    head = d.get("head", "?")
    model = d.get("model", "?")
    m = d.get("metrics", {}).get(head, {})

    # Highlight: any T315 substitution, plus the best and worst prediction.
    by_err = sorted(rows, key=lambda r: r["abs_error"])
    hi_names = {by_err[0]["mutant"], by_err[-1]["mutant"]}
    hi_names |= {r["mutant"] for r in rows if r["mutant"].upper().startswith("T315")}
    # Mean replicate half-width, for the agreement ribbon.
    widths = [
        abs(r["rep_high"] - r["rep_low"]) / 2
        for r in rows if r.get("rep_high") is not None and r.get("rep_low") is not None
    ]
    half = sum(widths) / len(widths) if widths else 0.0

    vals = [r["measured"] for r in rows] + [r["predicted"] for r in rows]
    lo, hi = min(vals), max(vals)
    pad = (hi - lo) * 0.06 or 0.01
    lo, hi = lo - pad, hi + pad

    W, H = 780, 700
    pad_l, pad_r, pad_t, pad_b = 92, 28, 112, 92
    pw, ph = W - pad_l - pad_r, H - pad_t - pad_b

    def sx(v: float) -> float:
        return pad_l + (v - lo) / (hi - lo) * pw

    def sy(v: float) -> float:
        return pad_t + ph - (v - lo) / (hi - lo) * ph

    p: list[str] = [
        f'<rect width="{W}" height="{H}" fill="{PAPER}"/>',
        f'<text x="24" y="34" font-size="19" font-weight="700" fill="{INK}" '
        f'font-family="IBM Plex Sans, Arial, sans-serif">'
        f'Predicted vs measured net growth under imatinib</text>',
        f'<text x="24" y="58" font-size="13" fill="{MUTED}" '
        f'font-family="IBM Plex Sans, Arial, sans-serif">'
        f'{len(rows)} held-out variants · positions never seen in training · '
        f'{esc(head)} on {esc(model)} embeddings</text>',
        f'<text x="24" y="80" font-size="13" font-weight="600" fill="{INK}" '
        f'font-family="JetBrains Mono, monospace">'
        f'Spearman {m.get("spearman", 0):.3f} · MAE {m.get("mae", 0):.5f} · '
        f'in-range {m.get("in_range_pct", 0):.1f}%</text>',
    ]

    # Agreement ribbon: +/- mean replicate half-width around y=x.
    if half > 0:
        band = [
            (sx(lo), sy(lo + half)), (sx(hi), sy(hi + half)),
            (sx(hi), sy(hi - half)), (sx(lo), sy(lo - half)),
        ]
        pts = " ".join(f"{x:.1f},{y:.1f}" for x, y in band)
        p.append(f'<polygon points="{pts}" fill="{BLUE}" opacity="0.08"/>')

    # Grid + ticks
    steps = 6
    for i in range(steps + 1):
        v = lo + (hi - lo) * i / steps
        x, y = sx(v), sy(v)
        p.append(f'<line x1="{x:.1f}" y1="{pad_t}" x2="{x:.1f}" y2="{pad_t+ph}" '
                 f'stroke="{RULE}" stroke-width="1"/>')
        p.append(f'<line x1="{pad_l}" y1="{y:.1f}" x2="{pad_l+pw}" y2="{y:.1f}" '
                 f'stroke="{RULE}" stroke-width="1"/>')
        p.append(f'<text x="{x:.1f}" y="{pad_t+ph+20}" font-size="11" '
                 f'text-anchor="middle" fill="{MUTED}" '
                 f'font-family="JetBrains Mono, monospace">{v:+.3f}</text>')
        p.append(f'<text x="{pad_l-10}" y="{y+4:.1f}" font-size="11" '
                 f'text-anchor="end" fill="{MUTED}" '
                 f'font-family="JetBrains Mono, monospace">{v:+.3f}</text>')

    # y = x
    p.append(f'<line x1="{sx(lo):.1f}" y1="{sy(lo):.1f}" x2="{sx(hi):.1f}" '
             f'y2="{sy(hi):.1f}" stroke="{MUTED}" stroke-width="2" '
             f'stroke-dasharray="6 4"/>')
    p.append(f'<text x="{sx(hi)-6:.1f}" y="{sy(hi)+18:.1f}" font-size="11" '
             f'text-anchor="end" fill="{MUTED}" '
             f'font-family="IBM Plex Sans, Arial, sans-serif">perfect agreement</text>')

    # Points: ordinary first, highlighted on top with a surface ring.
    for r in rows:
        if r["mutant"] in hi_names:
            continue
        p.append(f'<circle cx="{sx(r["measured"]):.1f}" cy="{sy(r["predicted"]):.1f}" '
                 f'r="2.6" fill="{BLUE}" opacity="0.4"/>')
    for r in rows:
        if r["mutant"] not in hi_names:
            continue
        cx, cy = sx(r["measured"]), sy(r["predicted"])
        p.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="6.5" fill="{AMBER}" '
                 f'stroke="{PAPER}" stroke-width="2"/>')
        p.append(f'<text x="{cx+11:.1f}" y="{cy+4:.1f}" font-size="12" '
                 f'font-weight="700" fill="{INK}" '
                 f'font-family="JetBrains Mono, monospace">{esc(r["mutant"])}</text>')

    # Axis titles
    p.append(f'<text x="{pad_l+pw/2:.1f}" y="{H-34}" font-size="13" '
             f'text-anchor="middle" fill="{MUTED}" '
             f'font-family="IBM Plex Sans, Arial, sans-serif">'
             f'measured net growth (netgrowth_mean_corrected)</text>')
    p.append(f'<text x="22" y="{pad_t+ph/2:.1f}" font-size="13" fill="{MUTED}" '
             f'text-anchor="middle" font-family="IBM Plex Sans, Arial, sans-serif" '
             f'transform="rotate(-90 22 {pad_t+ph/2:.1f})">predicted net growth</text>')

    # Legend
    ly = H - 14
    p.append(f'<circle cx="{pad_l+6}" cy="{ly-4}" r="4" fill="{BLUE}" opacity="0.6"/>')
    p.append(f'<text x="{pad_l+18}" y="{ly}" font-size="12" fill="{MUTED}" '
             f'font-family="IBM Plex Sans, Arial, sans-serif">held-out variant</text>')
    p.append(f'<circle cx="{pad_l+170}" cy="{ly-4}" r="5" fill="{AMBER}"/>')
    p.append(f'<text x="{pad_l+182}" y="{ly}" font-size="12" fill="{MUTED}" '
             f'font-family="IBM Plex Sans, Arial, sans-serif">'
             f'T315 · best · worst</text>')
    if half > 0:
        p.append(f'<rect x="{pad_l+340}" y="{ly-11}" width="16" height="10" '
                 f'fill="{BLUE}" opacity="0.14"/>')
        p.append(f'<text x="{pad_l+364}" y="{ly}" font-size="12" fill="{MUTED}" '
                 f'font-family="IBM Plex Sans, Arial, sans-serif">'
                 f'replicate agreement band (±{half:.4f})</text>')

    svg = (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
        f'viewBox="0 0 {W} {H}" aria-label="Scatter of predicted against '
        f'measured net growth for held-out variants. The dashed diagonal is '
        f'perfect agreement and the shaded band is the width within which a '
        f'prediction matches a repeat measurement.">'
        + "".join(p) + "</svg>"
    )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(svg)
    print(f"wrote {OUT} ({len(svg)//1024} KB, {len(rows)} points, "
          f"{len(hi_names)} highlighted)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
