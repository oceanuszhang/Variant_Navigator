# VariantNavigator

One pipeline that predicts **which drugs still work against a disease-causing mutation** — run against two domains with the same code:

- **Oncology:** EGFR variants vs. approved kinase inhibitors
- **Agriculture:** CYP51 variants in wheat pathogens vs. registered fungicides

## Pipeline

```
gene + variant + mode
  → fetch.py            UniProt: wildtype sequence + annotated binding site
  → mutant.py           apply the variant, validated against the sequence
  → esm2_score.py       ESM2 Δ log-likelihood (Modal)
  → boltz_structure.py  mutant structure prediction (Boltz)
  → drug_screen.py      Ptarmigan-1 binding screen across the drug panel
  → agent_report.py     Keenable agent → report in output/reports/
```

## Quickstart

```bash
pip install -r requirements.txt
cp .env.example .env        # then fill in keys

python main.py --gene EGFR --variant T790M --mode oncology
python main.py --test-case all
```

`fetch.py` and `mutant.py` need no credentials — they run on the stdlib. Stages
that are not wired up yet are **skipped with a note** rather than crashing, so
`main.py` always produces a report.

## Status

| Module | State |
| --- | --- |
| `fetch.py` | **Working** — tested against UniProt live |
| `mutant.py` | **Working** — subs (1- and 3-letter), multi-variant, in-frame del |
| `main.py` | **Working** — end-to-end with graceful stage skipping |
| `data/drug_panels/` | **Working** — 11 oncology + 12 fungicide drugs, SMILES from PubChem |
| `esm2_score.py` | Stub — interface fixed, awaiting Modal wiring |
| `boltz_structure.py` | Stub — interface fixed, awaiting Modal wiring |
| `drug_screen.py` | Stub — panel loading + ranking real, Ptarmigan-1 call pending docs |
| `agent_report.py` | Partial — deterministic markdown renderer works, Keenable call pending docs |

## Benchmark controls

The pipeline wins if it recovers these from first principles.

| Case | Expected |
| --- | --- |
| EGFR T790M | erlotinib, gefitinib **fail**; osimertinib, lazertinib **work** |
| EGFR L858R | erlotinib, gefitinib, osimertinib all **work** (negative control for resistance) |
| CYP51 Y137F | tebuconazole **fails**; see caveat below |

L858R matters as much as T790M: it catches a model that calls everything resistant.

Both panels deliberately include **off-target controls** (brigatinib in oncology;
azoxystrobin, pyraclostrobin, fluxapyroxad in fungicides). These hit different
targets entirely, so a pipeline that scores them as strong binders is reporting noise.

### Caveat on the agriculture control

`Y137F` is documented in **_Zymoseptoria tritici_** (wheat septoria), and as
`Y136F` in *Blumeria graminis* (powdery mildew). It is **not** an established
*Puccinia* (wheat rust) variant, and rusts are not primarily controlled by DMI
azoles. Options:

1. Reframe the agriculture demo as **wheat septoria** — keeps the literature support and the Y137F/tebuconazole resistance story intact. *Recommended.*
2. Keep wheat rust and pick a genuine *Puccinia* target plus a matching fungicide class.

The test case file carries this flag. Residue numbering must be confirmed against
whichever organism is chosen before it counts as a positive control.

## Verdict thresholds are not calibrated

`drug_screen.classify()` uses placeholder cutoffs (fail ≥ 1.5, partial ≥ 0.75
kcal/mol). Calibrate against the controls above before any verdict is trusted;
until then, treat output as a pipeline demo, not a prediction.

## Blocked / needs you

- **Ptarmigan-1 and Keenable API keys + docs.** No DNAForge repo is reachable from this account, and neither API has public docs, so both calls are stubs with a fixed interface. Drop the keys in `.env` (gitignored) and send the docs.
- **Repo name.** This is `Variant_Navigator`, not `variant-navigator` — kept as-is to avoid breaking your clone; say the word if you want it renamed.

## Scope note

Output is a computational prediction for research use. It is not a clinical
treatment recommendation or an agronomic spray recommendation.
