# VariantNavigator — working notes for Claude

## Web lookups

Use the Keenable tools (`search_web_pages`, `fetch_page_content`) instead of the
built-in web search for anything that needs the web: recent releases, prices,
docs, company facts.

## Interpreters — this matters, there are three

| Use | Interpreter |
| --- | --- |
| The pipeline, stdlib only | `python3` (3.11.4, python.org framework build) |
| ESM2 + the trained model | `.venv/bin/python` (torch, fair-esm, numpy, sklearn) |
| `modal` only | `/opt/anaconda3/bin/python3` (3.13.9) |

`python main.py …` works on the stdlib interpreter: stages needing torch skip
with a note rather than crashing. Use `.venv/bin/python main.py …` to include
ESM2.

## Non-negotiables in this repo

- **Never guess a numbering offset.** Four conventions are in play (clinical
  EGFR, Arabidopsis ALS, ABL1a/1b, rpoB +6). Verify against the real sequence or
  structure and let `strict=True` raise. A guessed offset silently corrupts
  every row — it nearly corrupted all 4,925 ABL1 variants, and the check caught
  it when the guess (+19) turned out to be 0.
- **Score against the noise ceiling, not 1.0.** Replicate *r* is 0.916, so the
  attainable ceiling is 0.978. Report "% of attainable".
- **Split by position, never by variant.** All ~19 substitutions at a residue
  share their structural features; a random split leaks and inflates.
  `--leaky` exists to demonstrate the inflation, not to be used.
- **A measurement beats a prediction.** `predict.py` returns the real measured
  value for the 4,925 variants in the screen; it only predicts for the other
  224.
- **Evidence never feeds prediction.** `evidence.py` runs after a prediction is
  fixed. T790M resistance is in thousands of papers, so a pipeline allowed to
  read them is looking up the answer, not deriving it.

## Boltz costs — measured, not guessed

`estimate-cost` is free (no resource, no GPU). Run it before every submission;
`boltz_structure.predict_structure()` does this for you and refuses to exceed
`VN_BOLTZ_MAX_USD` (default $0.25/call).

Measured 2026-11-14 with `estimate-cost`, model `boltz-2.1`:

| request | cost |
| --- | --- |
| full-length P00519 (1130 aa) + imatinib, 1 sample | $0.20 |
| kinase crop 242–512 (271 aa) + imatinib, 1 sample | $0.05 |
| kinase crop, **5 samples** | **$0.10** |
| kinase crop, 3 samples | $0.10 |
| kinase crop, 10 samples | $0.20 |

Three counter-intuitive findings, all of which change the right design:

- **Binding metrics are free.** Structure-only and structure+affinity cost the
  same, so always request `binding`.
- **5 samples cost the same as 3** and hit the per-unit floor ($0.020). Use 5
  and report the spread; 1 sample gives no spread to report.
- **There is a $0.05 floor.** The published "0–256 aa → $0.025" tier is
  unreachable — a 200-aa crop still bills $0.05. So there is **no cost reason
  to shave the crop**, and shaving it is what breaks position coverage: the DMS
  data spans 242–512 while 1IEP only resolves 225–498, so cropping to the
  structure would silently drop 14 measured positions to save nothing.

`admin:usage` returns 401 under the OAuth session (it wants an admin API key),
so the CLI cannot show a Boltz balance. Budget from the estimates, not a balance.

## Secrets

Real keys go in `.env`, which is gitignored. **Never** in `.env.example`, which
is committed — a key was pasted there once and caught before it reached a
commit. If a key does land in a tracked file, treat it as disclosed and rotate.
