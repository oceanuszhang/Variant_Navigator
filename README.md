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

For the ESM2 and model stages, use the project venv (torch is the only heavy
dependency, and it is optional — without it every caller skips):

```bash
python -m venv .venv
.venv/bin/python -m pip install torch fair-esm numpy scikit-learn

.venv/bin/python -m pipeline.esm2_score                            # cache marginals
.venv/bin/python scripts/train_baseline.py --esm esm2_t12_35M_UR50D  # the ladder
python main.py --gene P00519 --variant G250L --mode oncology       # one variant
```

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
| `esm2_score.py` | **Working** — ESM2 masked marginals, local on MPS, cached |
| `features.py` | **Working** — substitution / position / ESM blocks, ablatable |
| `predict.py` | **Working** — measurement when known, model otherwise |
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

## Sponsor tools — what each one actually did

Recorded honestly, including the one that could not do what we hoped. Every
claim below was verified by running it, not by reading a landing page.

| Tool | Used for | Outcome |
| --- | --- | --- |
| **Talus Bio Ptarmigan-1** | Compound→target engagement screen | ✅ Ran 2 screens, $0.06. **Cannot score mutants** — see below |
| **Modal** | ESM-2 on A10G GPU, persistent weights Volume | ✅ Smoke-tested end to end; makes the 650M model practical |
| **ESM-2** (`fair-esm`) | Variant scoring, two different ways | ✅ The single most valuable feature block (r 0.245 → 0.3315) |
| **Paperclip** | Literature evidence, isolated from prediction | ✅ Returns live PMC citations — was silently broken, fixed |
| **Boltz** | Mutant structure prediction (stage 4) | ⏳ CLI v0.44.0 installed, oauth authenticated, not yet called |
| **Keenable** | Web search / page fetch for agents | ⏳ MCP server ✔ Connected; tools need a Claude Code restart |
| **Genomic Intelligence** | Ground-truth cross-referencing | ⏳ `.env` placeholder added, awaiting key |

### Ptarmigan-1 cannot score mutants, and that is structural

Their documentation is explicit:

> *"Only the human proteins listed there can be screened; there is no way to add
> one."*

Confirmed against the live API: EGFR returns exactly one entry,
`uniprot_accession: P00533` — canonical wild-type, no variant field. So
`drug_screen.py`'s design, `screen_panel(gene, variant, wt_sequence,
mt_sequence, …)` → mutant-minus-wildtype delta, is **unimplementable**. T790M and
wild-type EGFR return identical results, making the delta structurally zero.
The API key was never the blocker; the API shape is.

What it *can* do, it did well — as a **panel-design control**:

| | ABL1 screen | EGFR screen |
| --- | --- | --- |
| **imatinib** | **#1** of 12 (0.5572) | **#12** of 12 (0.4744) |
| **brigatinib** (off-target control) | **#12** (0.4876) | #11 (0.5139) |

Imatinib tops the ABL1 screen and comes last on EGFR; the deliberately planted
off-target control sinks on both. One prediction is genuinely sharp —
neratinib peaks at residue **797C**, the real covalent anchor for EGFR
inhibitors, and erlotinib/gefitinib/icotinib cluster at **744** in the ATP site.

Two caveats kept attached to that data in `output/ptarmigan/`:

- **Narrow dynamic range.** 0.4876–0.5572 (ABL1), 0.4744–0.5850 (EGFR), all near
  0.5, and their docs warn against reading 0.5 as a cutoff. The *ranking*
  carries the signal; the absolute values do not. Both targets also exceed the
  model's 766-residue training length, so both are flagged out-of-distribution.
- **Peak residues are about half wrong.** Only 8/12 (ABL1) and 6/12 (EGFR) land
  inside the kinase domain. `gefitinib → 99W` on ABL1 and `dacomitinib → 535C`
  in the EGFR *extracellular* domain are not credible sites.

### Modal: what it is for here

The 35M model runs fine on local MPS and produced the headline numbers. Modal
earns its place on the **650M** model, a 2.5 GB download that stalled twice on
local wifi: on Modal the weights persist in a Volume and the passes run on an
A10G.

```bash
/opt/anaconda3/bin/modal run scripts/modal_esm.py                      # 650M
/opt/anaconda3/bin/modal run scripts/modal_esm.py --model esm2_t12_35M_UR50D
```

`modal` is installed only under `/opt/anaconda3` (Python 3.13), not the venv —
see "Environment gotchas".

## Two ways to use ESM-2, and why both exist

| | Masked marginals | Mutant-marginal difference |
| --- | --- | --- |
| Feature | `log P(alt) − log P(ref)` at masked *i* | `h_mut[i] − h_ref[i]` |
| Dimensions | 5 | 480 (35M) / 1280 (650M) |
| Forward passes | **271** (one per position) | **4,926** (one per variant) |
| Wired into `main.py` | ✅ yes | research script |
| Script | `pipeline/esm2_score.py` | `scripts/embed_variants.py` |

Both read **one residue**, never a pooled mean. That is the whole trick: all
4,925 variants are the same protein at 99.82% sequence identity, so a
mean-pooled whole-sequence embedding returns a near-identical vector for every
one of them and a model fit on it has nothing to separate. Subtracting the
reference removes whatever is constant across a position's ~19 substitutions,
leaving the substitution's own effect (measured mean `|feature|` = 0.139).

### Evaluation protocol for the embedding heads

```bash
.venv/bin/python scripts/embed_variants.py              # one-off, cached
.venv/bin/python scripts/train_embeddings.py            # ridge + hist_gbm
```

80/20 split **grouped by residue position**, fixed seed, asserting zero shared
positions rather than trusting it. Reported: MAE, Spearman, Pearson, and an
**in-range** metric — the share of validation variants whose prediction lands
between that variant's two corrected replicates.

**Read in-range carefully.** The replicate interval has a mean width of 0.0069
against a response spanning ~0.12, so landing inside it is a demanding test. A
low percentage means the model is not yet as precise as repeating the
experiment, not that it is broken. The width is printed next to it so the two
are never confused.

### The "518-aa reference" trap, avoided

The implementation plan this came from specified a *518-aa ABL1 reference*.
**P00519 is 1130 aa.** Applying position 315 to a 518-residue sequence mutates
the wrong residue in every row. `embed_variants.py` validates every variant
against the reference residue before spending any compute and refuses the whole
run on a mismatch — the difference between a result and 4,925 silently corrupted
ones. This is the fourth numbering trap in this repo and the third that would
have passed silently.

## The model — exactly what it is

### One sentence

A **ridge regression** on 32 features — 23 substitution-chemistry, 4
pocket-geometry and 5 **ESM2 masked-marginal** — predicting Ba/F3 net growth
under 1200 nM imatinib, trained on 4,925 measured BCR-ABL variants and
evaluated by **position-grouped** 5-fold cross-validation.

```bash
# the honest evaluation ladder
.venv/bin/python scripts/train_baseline.py --esm esm2_t12_35M_UR50D

# one variant, end to end
python main.py --gene P00519 --variant G250L --mode oncology
```

### Results

5-fold position-grouped CV, pooled out-of-fold, ESM2-35M:

| features | model | r | Spearman | MAE | n feats |
| --- | --- | --- | --- | --- | --- |
| global mean | — | −0.078 | −0.072 | 0.0207 | — |
| position only | hist_gbm | 0.108 | 0.139 | 0.0206 | 4 |
| substitution only | ridge | 0.233 | 0.243 | 0.0196 | 23 |
| all (no ESM) | random_forest | 0.245 | 0.238 | 0.0196 | 27 |
| **ESM only** | ridge | **0.291** | 0.293 | 0.0192 | **5** |
| substitution + ESM | ridge | 0.331 | 0.342 | 0.0188 | 28 |
| **all + ESM** | **ridge** | **0.3315** | **0.339** | **0.0189** | 32 |

Reported against the **attenuated ceiling of r = 0.978**, not against 1.0: with
two noisy replicates (replicate r = 0.916), no model can exceed that. So
**0.3315 is 34% of attainable**, at MAE 2.7× the replicate noise floor.

### Three findings from the ladder

1. **ESM2 alone, 5 features, beats all 27 hand-engineered features combined**
   (0.291 vs 0.245). It is the single most valuable block.
2. **Pocket geometry contributes almost nothing.** 0.233 → 0.2345 when added to
   substitution chemistry, and 0.3312 → 0.3315 when added on top of ESM. Three
   decimal places of nothing. Distance alone against all 4,665 structurally
   covered variants is **r = 0.069**.
3. **Ridge beats gradient boosting and random forest at every rung.** With 4,925
   rows and ~32 features the relationship is close to linear, and the extra
   capacity only fits noise. XGBoost was deliberately not installed — sklearn's
   `HistGradientBoostingRegressor` was already present and lost anyway.

### Why ESM2 is used as masked marginals, not embeddings

Every variant is the same protein, so any two mutant sequences are **99.82%
identical**. A mean-pooled whole-sequence embedding returns a nearly identical
vector for all 4,925 variants and a model fit on it has nothing to separate.
Instead, per ESM-1v, we mask the position and read the model's own distribution:

```
score(ref → alt at i) = log P(alt | x masked at i) − log P(ref | x masked at i)
```

One forward pass per **position** (271), not per variant, covering all 4,925
variants. Five features come out: `delta_ll` and `ll_alt` vary per substitution,
while `ll_ref` and `entropy` are position-level conservation measures, plus a
`has_esm` flag so unscored positions are a declared gap rather than a silent 0.

**Length limit, which bites here.** ESM2 has 1024 position embeddings and ABL1
P00519 is 1130 aa, so the full sequence will not fit. We crop to residues
**1–1022**, which still contains the whole kinase domain (242–512) with ~240
residues of leading context. The crop is recorded in the cache key, so scores
can never be silently compared across different crops.

### ESM2 measures conservation, not resistance

Worth internalising before trusting it. The masked-marginal scores for known
imatinib-resistance mutations:

```
T315I   delta_ll = -2.901      the classic gatekeeper escape
Y253H   delta_ll = -4.354
M351T   delta_ll = -8.430
H396R   delta_ll = +1.023
G250E   delta_ll = +2.776
```

ESM calls the most famous resistance mutations **unlikely**, because those sites
are conserved. It is scoring evolutionary plausibility, not drug escape. It still
helps a great deal, and the reason is mechanistic: Ba/F3 cells here **depend on
BCR-ABL to survive**, so net growth falls when the kinase is broken and rises
when the drug is evaded. ESM predicts the first axis well and the second not at
all — which is precisely why it is complementary rather than sufficient.

### What is deliberately not claimed

- **No log₂ fold change.** `baseline.calibrated()` returns `False`. The matrix
  cell in `figures/prediction_matrix.svg` needs a dose-ratio regression; the
  IC₅₀ table that would support it is fetched but not yet wired.
- **A measurement beats a prediction.** `predict.py` returns the real measured
  value for the 4,925 variants in the screen and only predicts for the
  remaining 224. Re-predicting known values would flatter the demo.
- **Ranking, not thresholds.** At r = 0.33 the output orders variants; it does
  not support a clinical cutoff, and every report says so.

### Is ESM2 wired? Yes — here is what it needed

It runs locally on Apple Silicon via MPS — **no Modal, no GPU rental, no API
key**. It needed `torch` and `fair-esm`, which were absent from both system
interpreters, so they live in a project venv:

```bash
python -m venv .venv
.venv/bin/python -m pip install torch fair-esm numpy scikit-learn
.venv/bin/python -m pipeline.esm2_score          # 271 positions, seconds on MPS
```

Scores are cached to `data/esm/` keyed by model **and crop**, so a second run is
instant and a stale crop can never be mistaken for a fresh one. `available()`
returns `False` without torch and every caller skips rather than crashing, so
`python main.py` still works on the stdlib-only interpreter.

Two model sizes are wired. `esm2_t12_35M_UR50D` is the fast default used for the
numbers above; `esm2_t33_650M_UR50D` is the better model and the one
`esm2_score.DEFAULT_MODEL` points at.

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
├── features.py     # substitution / position / ESM blocks            [working]
├── esm2_score.py   # ESM2 masked marginals, local on MPS, cached     [working]
├── predict.py      # measurement if known, trained model otherwise   [working]
├── evidence.py     # Paperclip retrieval, isolated from prediction   [working]
├── boltz_structure.py # Boltz-2 via Tamarind Bio                     [stub]
├── drug_screen.py  # Ptarmigan-1; key present, awaiting API docs     [stub]
└── agent_report.py # report writer; Keenable agent to be wired       [partial]
scripts/
├── make_figures.py    # regenerates figures/ from live pocket geometry
└── train_baseline.py  # the ablation ladder under position-grouped CV
figures/            # SVG, written straight from the code that computes them
data/esm/           # cached masked marginals, keyed by model AND crop
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
| **PubChem PUG REST** | 4 | Source of SMILES for all 37 panel drugs | No |
| **RCSB PDB** | — | Biological assemblies for pocket geometry | No |
| **BCR-ABL DMS** (MIT) | — | 4,925 measured variants: the training labels | No |
| **ESM2** (`fair-esm`) | 2 | Masked-marginal variant scores, local on MPS | **No** |
| **scikit-learn** | — | Ridge / HistGBM / RandomForest + the CV harness | No |
| **Paperclip** (CLI) | 5 | Literature evidence, isolated from prediction | No |
| **Python stdlib** (`urllib`, `re`, `dataclasses`) | 1 | Variant parsing + mutant generation | No |

Note what is *absent* from that list: **no Modal, no GPU rental, no inference
API**. ESM2 runs on the laptop. The only remaining credentialed dependency on
the prediction path is Ptarmigan-1, which is optional.

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
