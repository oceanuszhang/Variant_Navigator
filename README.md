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
| `main.py` | **Working** — 6 stages, end-to-end with graceful stage skipping |
| `pocket.py` | **Working** — 13/13 reference checks pass on biological assemblies |
| `dms.py` | **Working** — 4,925 measured BCR-ABL variants, fetched and cached |
| `baseline.py` | **Working** — the only predictive stage that needs no credentials |
| `evidence.py` | **Working** — Paperclip returns live PMC citations |
| `data/drug_panels/` | **Working** — 11 oncology + 12 fungicide + 13 antibiotic drugs, SMILES from PubChem |
| `esm2_score.py` | Stub — **blocked on `torch`, not on a key**; see "Is ESM2 wired?" |
| `boltz_structure.py` | Stub — interface fixed, awaiting Modal wiring |
| `drug_screen.py` | Stub — panel loading + ranking real, Ptarmigan-1 call pending docs |
| `agent_report.py` | Partial — deterministic markdown renderer works, Keenable call pending docs |

## Data sources and structure

Every accession, CID and offset below was verified live against the source API,
not copied from the plan doc.

### Targets (`data/targets.json`)

| | Experiment 1 | Experiment 2 |
| --- | --- | --- |
| Target | Human EGFR | Rice ALS1 |
| UniProt | `P00533` (1210 aa) | `Q6K2E8` (644 aa) |
| Organism | *Homo sapiens* | *Oryza sativa*, weedy accessions |
| Pocket structures | PDB `1M17`, `4ZAU` | PDB `1Z8N`, `1YHZ` (*Arabidopsis*) |
| Assay | Ba/F3 cell viability IC₅₀ | fresh shoot biomass GR₅₀, 21 d |
| Reference | wild-type EGFR | IRGA417 susceptible cultivar |
| Compounds | erlotinib, osimertinib | imazethapyr, nicosulfuron |
| Numbering offset | 0 | **−26** |

### ⚠ The ALS numbering offset will silently corrupt Experiment 2

The ALS resistance literature numbers variants in **Arabidopsis** convention
(A122, W574, S653). Rice ALS1 is 26 residues shorter. Checked against the real
sequences:

| Published | *Arabidopsis* `P17597` | Rice `Q6K2E8` | Rice, renumbered −26 |
| --- | --- | --- | --- |
| A122 | A ✓ | **A ✓ (coincidence)** | A96 ✓ |
| W574 | W ✓ | D ✗ | W548 ✓ |
| S653 | S ✓ | past 644-aa end ✗ | S627 ✓ |

**A122 is alanine in both numbering systems by accident.** A122T is the headline
live demo case, so it would validate fine while every other site was quietly
wrong — the worst kind of bug, because the demo passes. W548L and S627I are the
documented rice ALS resistance substitutions, which is the independent check
that −26 is right.

`mutant.apply_variant(..., offset=-26)` handles this, and
`strict=True` refuses rather than guessing when a residue disagrees. Genotypes
stay in published numbering in the CSVs so transcription remains checkable
against the paper; renumbering happens in code.

### Compounds (`data/compounds/compounds.json`)

SMILES resolved from **PubChem PUG REST** (no key, no licence friction):
erlotinib `CID 176870`, osimertinib `CID 71496458`, imazethapyr `CID 54740`,
nicosulfuron `CID 73281`.

### Measured response data — **the blocker is gone: Experiment 3**

The original plan needed two files hand-transcribed from published figures
(`egfr_ic50.csv`, `rice_gr50.csv`). Nothing downstream of feature extraction
could run until they existed, and no API key could supply them.

**That is resolved.** Experiment 3 uses a public, MIT-licensed deep mutational
scanning dataset that is fetched programmatically and is ~1000× larger than
what the transcription would have produced:

| | Experiment 3 |
| --- | --- |
| Target | Human ABL1 (BCR-ABL) kinase domain, residues 242–512 |
| UniProt | `P00519` (1130 aa) |
| Compound | imatinib (PubChem `CID 5291`) |
| Assay | **Ba/F3 net growth rate** at 1200 nM imatinib, 6 d, TileSeq + duplex seq |
| Measured variants | **4,925**, each with 2 replicates |
| Pocket structure | PDB `1IEP` (ligand `STI`), offset `0` |
| Numbering offset | 0 — verified, see below |
| Source | [pritchardlabatpsu/abl_dms](https://github.com/pritchardlabatpsu/abl_dms) (MIT) |

> Inam H.\*, Tomaszkiewicz M.\*, Reynolds J.A., Yang Z., Leighow S.M.,
> Pritchard J.R. *Characterizing Variants of Uncertain Drug Response (VUDRs)
> Using Quantitative Measurements at Clinical Exposures.* (\*equal contribution)

It is the same **Ba/F3 assay system** as Experiment 1, so the schema carries
over. `pipeline/dms.py` fetches, caches and validates it:

```bash
python -m pipeline.dms
```

Richer tables from the same source, not yet wired: per-variant IC₅₀s fitted
across 300/600/1200 nM (`Figure_3D/all_mutants_with_IC50_estimates.csv`) and
resistance calls at three clinical doses (`Figure_4C/`). The IC₅₀ table is the
closer match to the original log₂(dose ratio) target and is the natural next
step.

#### A fourth numbering system, and it has a trap too

ABL1 has two conventions: **ABL1a**, where the gatekeeper is Thr315, and
**ABL1b**, which is 19 residues longer so the same residue is 334. Checked
against the real sequences and structures rather than assumed:

| Source | Convention | Offset to ABL1a |
| --- | --- | --- |
| The DMS dataset | ABL1a | 0 |
| UniProt `P00519` (1130 aa) | ABL1a | **0** — verified at T315/Y253/E255/G250/M351/H396 |
| PDB `1IEP`, `2HYY` | ABL1a | 0 |
| PDB `1OPJ` | ABL1b | **+19** |

`1OPJ` carries the *same* imatinib ligand (`STI`) as `1IEP`, so substituting it
looks harmless and silently returns the wrong residue's distance. The initial
guess here was that UniProt would need `+19`; the sequence check proved it needs
`0`. **Guessing the offset would have corrupted all 4,925 rows.**

### What the measured data taught us

Having real labels immediately falsified the headline baseline. Against all
4,665 variants with both a measurement and a structure:

| Metric | Value |
| --- | --- |
| Pearson *r* (pocket distance vs. net growth) | **+0.069** |
| Replicate *r* (the noise floor) | 0.916 |
| Attenuated *r* ceiling (best a perfect model could score) | **0.978** |
| Replicate MAE | 0.0069 |

**Pocket distance alone is not predictive.** On four hand-picked EGFR control
cells it looked like 3/4; on 4,665 real measurements it explains under 0.5% of
the variance. Both facts are true, and the second one is the one that counts.
This is the clearest argument in the repo for benchmarking against measured data
rather than curated anecdotes.

The reason is visible in the most-resistant variants:

```
F359A  netgrowth=+0.0699  dist= 4.08 A     contacts the drug
H396P  netgrowth=+0.0488  dist=15.54 A     15 A away and still resistant
T315L  netgrowth=+0.0481  dist= 2.88 A     the gatekeeper
```

`H396` sits in the activation loop and escapes imatinib **allosterically**, by
shifting the kinase to a conformation imatinib cannot bind. Geometry cannot see
a mechanism that does not touch the ligand.

Two further facts that constrain any model built on this data:

- **59% of the response spread is *within* a single position.** So a feature
  that is constant across a position — distance, burial, conservation — is
  mathematically incapable of explaining most of the variance. The substitution
  identity has to carry it.
- **The response confounds two mechanisms.** Ba/F3 cells here depend on BCR-ABL
  to survive, so net growth falls both when the kinase is *broken* and when the
  drug *works*. A single scalar cannot separate "kinase dead" from "drug
  evaded", and the two have opposite clinical meaning.

### Is ESM2 wired? No — and it is not blocked on a key

`esm2_score.py` is still a stub that raises `NotImplementedError`. It is blocked
on **`torch`, which is not installed in either interpreter on this machine** —
not on credentials.

When it is wired, the embedding must be used correctly, because the obvious
approach fails here:

- All 4,925 variants are **the same protein**. There is one reference sequence,
  not 4,925, and any two mutant sequences are **99.82% identical**.
- A mean-pooled whole-sequence embedding therefore returns an almost identical
  vector for every variant. A linear model on it has no signal to fit.

What works instead, in increasing cost:

1. **Masked-marginal log-likelihood** (zero-shot, no training):
   `score = log P(alt | mask at i) − log P(ref | mask at i)`. One forward pass
   per *position* (270), one scalar per variant. This is the standard ESM-1v /
   ESM2 variant-effect method.
2. **Per-residue hidden state at the mutated position**, not the pooled mean,
   concatenated with an encoding of ref→alt. Varies per variant, and trainable
   with ridge regression.
3. **Mutant-marginal**: embed each mutant sequence and read the hidden state at
   the mutated site. 4,925 forward passes, tractable with the 35M model.

And a caveat worth stating up front: **ESM scores evolutionary likelihood, not
drug resistance.** It should predict the "is the kinase still functional" axis
well, since that is conservation, and predict drug escape poorly — T315I is both
highly conserved and the most famous imatinib-resistance mutation. Expect ESM to
complement pocket geometry, not to subsume it.

### Repo structure

```
data/
├── targets.json              # accessions, PDB refs, numbering offsets, assays
├── compounds/compounds.json  # SMILES + PubChem CIDs, tagged by experiment
├── datasets/
│   ├── SCHEMA.md             # column spec for the original two-experiment design
│   └── abl_dms_netgrowth_1200nM.csv   # 4,925 measured variants (fetched, MIT)
├── structures/               # PDB biological assemblies (gitignored cache)
├── drug_panels/              # from the earlier classification design (see note)
└── test_cases/               # positive controls (see note)
pipeline/
├── fetch.py        # UniProt -> sequence + binding-site residues     [working]
├── mutant.py       # variant notation -> mutant seq, with renumber() [working]
├── pocket.py       # min heavy-atom distance to the bound ligand     [working]
├── dms.py          # BCR-ABL DMS loader, noise floor, position CV    [working]
├── baseline.py     # distance -> contact score -> call; no keys      [working]
├── evidence.py     # Paperclip retrieval, isolated from prediction   [working]
├── esm2_score.py   # ESM-2 via Modal; blocked on torch, not a key    [stub]
├── boltz_structure.py # Boltz-2 via Tamarind Bio                     [stub]
├── drug_screen.py  # Ptarmigan-1; key present, awaiting API docs     [stub]
└── agent_report.py # report writer; Keenable agent to be wired       [partial]
scripts/
└── make_figures.py # regenerates figures/ from live pocket geometry
figures/            # SVG, written straight from the code that computes them
```

### Figures

Regenerate with `python scripts/make_figures.py`. Each is written directly from
`pipeline.pocket`, so a figure cannot drift from the code.

| | |
| --- | --- |
| [`pocket_distances.svg`](figures/pocket_distances.svg) | variant-site→inhibitor distance, with the 5 Å contact threshold |
| [`prediction_matrix.svg`](figures/prediction_matrix.svg) | the output object: variant × compound, and how to read it |
| [`numbering_systems.svg`](figures/numbering_systems.svg) | one residue across every numbering convention |

### Environment gotchas

Two that cost real debugging time on this machine:

- **`paperclip` was silently dark.** `evidence.py` searched only `$HOME` for an
  interpreter carrying Paperclip's dependencies, but this machine's anaconda is
  in `/opt`. `available()` returned `False` behind a bare
  `ModuleNotFoundError: click`. Fixed by searching `/opt` first.
- **`modal` is installed under the wrong interpreter.** `modal 1.6.4` lives in
  `/opt/anaconda3` (Python 3.13), while `python3` resolves to the python.org
  3.11 framework build, where `import modal` fails. Either run the pipeline as
  `/opt/anaconda3/bin/python3 main.py …` — the whole pipeline imports fine
  there — or `pip install modal` into 3.11.

## Tools

### In use now (working, committed)

| Tool | Step | Function | Key needed |
| --- | --- | --- | --- |
| **UniProt REST API** | 1 | Wildtype sequence + annotated binding-site residues | No |
| **PubChem PUG REST** | 4 | Source of SMILES for all 36 panel drugs | No |
| **Python stdlib** (`urllib`, `re`, `dataclasses`) | 1 | Variant parsing + mutant generation | No |

### Planned (interfaces fixed, awaiting credentials or wiring)

| Tool | Step | Function | Key needed |
| --- | --- | --- | --- |
| **ESM2 (650M)** | 2 | Evolutionary tolerance / functional impact score | via Modal |
| **Modal** | 2, 3 | GPU compute backend for ESM2 + Boltz | `MODAL_TOKEN_ID`, `MODAL_TOKEN_SECRET` |
| **Boltz** | 3 | Mutant structure prediction + pocket RMSD | via Modal |
| **Ptarmigan-1** (Talus Bio) | 4 | Drug–mutant binding affinity screening | `PTARMIGAN_API_KEY` **+ docs** |
| **Keenable** | 5 | Agent orchestration + report generation | `KEENABLE_API_KEY` **+ docs** |

### Ground-truth references (Step 5 cross-referencing)

| Source | Domain | Access | Key needed |
| --- | --- | --- | --- |
| **OncoKB** | Cancer | REST API, registration required | `ONCOKB_API_TOKEN` |
| **ClinVar** | Cancer | NCBI E-utilities | No (`NCBI_API_KEY` raises rate limit) |
| **CARD** | Antimicrobial | Bulk download + local lookup | No (free account for download) |
| **FRAC** | Agriculture | **No API** — published tables, manual curation | n/a |

### Considered and not used

| Tool | Why |
| --- | --- |
| **BioPython** | Spec assigned it Step 1 variant application. Not needed — `mutant.py` does this on the stdlib, and skipping the dependency keeps the demo install trivial. Still listed in `requirements.txt` as optional; drop it if nothing else claims it. |
| **DrugBank / ChEMBL** | Spec's SMILES source. Used PubChem instead: no key, no licence friction, and it resolved all 36 drugs first try. DrugBank requires a licence for bulk/commercial use. |
| **AlphaGenome2** | Contrast only, per the plan — it returns expression predictions, not treatment efficacy. Not a pipeline dependency. |

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

- **Ptarmigan-1 and Keenable API keys + docs.** No DNAForge repo is reachable from this account, and neither API has public docs, so both calls are stubs with a fixed interface. Drop the keys in `.env` (gitignored) and send the docs. **These two are the critical path** — Steps 4 and 5 cannot run without them.
- **OncoKB API token**, if you want Step 5 to cross-reference sensitivity calls programmatically rather than by hand.
- **"Paperclip"** — not in the plan doc and nothing by that name turns up in protein/drug-binding tooling. Send a link or say what it does and I will slot it in; see the note below.
- **Repo name.** This is `Variant_Navigator`, not `variant-navigator` — kept as-is to avoid breaking your clone; say the word if you want it renamed.

## Scope note

Output is a computational prediction for research use. It is not a clinical
treatment recommendation or an agronomic spray recommendation.
