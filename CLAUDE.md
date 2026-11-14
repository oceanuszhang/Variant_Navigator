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

## Secrets

Real keys go in `.env`, which is gitignored. **Never** in `.env.example`, which
is committed — a key was pasted there once and caught before it reached a
commit. If a key does land in a tracked file, treat it as disclosed and rotate.
