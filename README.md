# VariantNavigator

One pipeline that predicts **which drugs still work against a disease-causing mutation** — run against two domains with the same code:

- **Oncology:** EGFR variants vs. approved kinase inhibitors
- **Agriculture:** CYP51 variants in wheat pathogens vs. registered fungicides
- **Antimicrobial:** rpoB variants in *M. tuberculosis* vs. TB drug regimens

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
python main.py --test-case all          # all five controls, 3 domains
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
| `data/drug_panels/` | **Working** — 11 oncology + 12 fungicide + 13 antibiotic drugs, SMILES from PubChem |
| `esm2_score.py` | Stub — interface fixed, awaiting Modal wiring |
| `boltz_structure.py` | Stub — interface fixed, awaiting Modal wiring |
| `drug_screen.py` | Stub — panel loading + ranking real, Ptarmigan-1 call pending docs |
| `agent_report.py` | Partial — deterministic markdown renderer works, Keenable call pending docs |

## Benchmark controls

Five cases across three domains; `python main.py --test-case all` runs them all.

| Domain | Gene | Variant | Expected |
| --- | --- | --- | --- |
| Oncology | EGFR | T790M | erlotinib, gefitinib **fail**; osimertinib, lazertinib **work** |
| Oncology | EGFR | L858R | erlotinib, gefitinib, osimertinib all **work** |
| Agriculture | CYP51 | Y137F | tebuconazole **fails** |
| Antimicrobial | rpoB | S450L | all rifamycins **fail**; off-target drugs **unchanged** |
| Antimicrobial | rpoB | D435V | rifampicin **fails**; rifabutin **disputed** (stretch case) |

L858R matters as much as T790M: it catches a model that calls everything resistant.

Every panel includes **off-target controls** — drugs that hit a different protein
entirely (brigatinib; azoxystrobin, pyraclostrobin, fluxapyroxad; isoniazid,
moxifloxacin, bedaquiline and the rest of the non-rifamycins). They must come back
*unchanged*. A pipeline that predicts a binding shift for isoniazid against an
rpoB variant is reporting noise, and the off-target controls are how you see that.

### Two controls were changed, and one rejected

**Rejected: PBP2a S403L (methicillin fails / ceftaroline works).** Two problems.
Ser403 is PBP2a's *catalytic active-site serine*; S403L would abolish
transpeptidase activity rather than confer resistance, and it is not a documented
resistance substitution. More fundamentally, methicillin resistance in MRSA comes
from **acquiring** mecA/PBP2a at all — wildtype PBP2a already binds beta-lactams
poorly. A pipeline that scores mutant-vs-wildtype delta therefore *cannot*
recover "methicillin fails": the delta is ~zero because the wildtype is already
resistant. The control would fail by construction, not because the science is
wrong. Documented ceftaroline-resistance substitutions are Y446N/E447K (binding
pocket) and N146K/N204K/E239K (allosteric).

**Replaced with: rpoB S450L.** Rifamycins bind RpoB directly, wildtype is
sensitive, and the variant abolishes binding — so mutant-vs-wildtype delta is
exactly the right signal. It is the most frequent rifampicin-resistance
substitution in the RRDR and is catalogued in CARD.

**Agriculture reframed from wheat rust to wheat septoria.** Y137F is documented
in *Zymoseptoria tritici* and as Y136F in *Blumeria graminis* — not in *Puccinia*
(wheat rust), which is not primarily controlled by DMI azoles. Septoria keeps the
literature support and the tebuconazole story intact. Still a valid agricultural
resilience demo: septoria is the most economically damaging wheat disease in
Europe.

### Numbering conventions are a live hazard

Clinical nomenclature and UniProt numbering disagree for rpoB. The TB
literature's **S450L is position 456** in UniProt P9WGY9 — an offset of **+6**.
Applied naively, S450L lands on Thr450 and silently mutates the wrong residue.

```
$ python main.py --gene P9WGY9 --variant S450L --mode antimicrobial
error: S450L: sequence has T at position 450, not S
```

`strict=True` in `apply_variant()` catches this instead of producing a confident
wrong answer. Test cases carry a `numbering_offset`, and `mutant.renumber()`
translates clinical notation into UniProt numbering.

The same hazard bites allele choice: UniProt holds three *Z. tritici* CYP51
entries, and one (B2L2K3) has **Lys** at 137, so Y137F against it fails
validation. Test cases pin an explicit accession rather than trusting gene-name
search.

## A delta alone is not an efficacy prediction

`drug_screen` currently reports mutant-minus-wildtype affinity. That is the right
primary signal, but on its own it misreads any case where the **wildtype is
already a poor binder** — the delta is ~0 and the drug scores as "works" when it
never worked. This is what sinks the MRSA framing above, and it is why the screen
needs to carry *absolute* mutant affinity alongside the delta, with a
baseline-sensitivity gate before any verdict. Worth building in before Step 4
integration rather than after.

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
