# VariantNavigator

Cancer-only research workflow for exploring mutation effects on drug response.
The supported targets are human **ABL1/BCR-ABL** and **EGFR**. Quantitative
measured response and runtime modeling cover **ABL1 + imatinib only**.
This is a research prototype, not a treatment recommender.

## Implementation status (updated 2026-11-14)

| Component | Status | Verification |
| --- | --- | --- |
| Cancer-only orchestration | Implemented | Offline ABL1 report and service-mocked tests pass |
| Modal authentication | Verified in existing `.venv` | Read-only app listing and actual GPU job succeeded |
| ESM2 35M GPU smoke test | Passed | Four validated variants returned a finite `(4, 480)` feature matrix |
| Saved regression models | Implemented and smoke-tested | Both heads reload with identical held-out predictions |
| Model provenance/loader | Verified | Nine regression tests pass, including checksum, dimension, partial-run, and optional-service guards |
| Single-variant saved ESM inference | Passed live verification | Fresh T315I Modal embedding exactly matches saved batch prediction; separate from measured response |
| Full 4,925-variant GPU embedding run | Passed | Finite `(4925, 480)` matrix; evaluated heads saved and reloaded successfully |
| Boltz | Guarded optional integration | Live authentication/control run pending |
| Literature/agent adapters | Incomplete | Deterministic reports work; evidence remains separate |
| UI | Deferred | No UI work in this phase |

## Run the backend

Python 3.11+ is sufficient for the core workflow. Run from the repository root.

```powershell
python main.py --gene ABL1 --variant T315I --offline --no-agent
python main.py --gene EGFR --variant T790M --no-agent
python main.py --test-case all --no-agent
python -m unittest discover -s tests -v
```

The first command is deterministic and needs no credentials, network, GPU, or
third-party packages. It uses the bundled 518-aa ABL1 study reference and assay
CSV. Offline mode skips ESM, structure retrieval, Boltz, and novel inference.
EGFR offline lookup is not supported. Invalid inputs fail explicitly.

`--test-case all` runs only the oncology fixtures; completion means the workflow
ran, not that any uncalibrated efficacy verdict passed a benchmark. Historical
agricultural and antimicrobial data remain in the repository but are outside
the active workflow. Other targets, modes, panels, nonhuman taxa, and numbering
offsets are rejected by the entry point.

## Workflow and output contract

1. Resolve an explicit human accession (ABL1/P00519 or EGFR/P00533).
2. Validate and apply the published variant at offset 0.
3. Optionally compute an ESM2 sequence score when its dependencies are present.
4. Optionally predict an ABL1/imatinib complex with Boltz (`--boltz`).
5. Compute experimental reference-pocket geometry online.
6. For single ABL1 substitutions at 242-512, return a measured assay result if
   present; otherwise attempt the handcrafted-feature runtime model online.
7. Attach a saved held-out ESM estimate separately when available.
8. Write UTF-8 Markdown and JSON reports in ignored `output/reports/`.

JSON includes `schema_version`, `sequence_source`, `offline`, `status`, and
`skipped_stages`. Optional-service failures are preserved with their reasons.
A partial report is expected when integrations are unavailable; it does not
mean all predictive stages executed. Boltz output retains confidence, sample
spread, binding confidence, crop, cost, and prediction ID. Binding confidence
is not binding affinity or efficacy.

## Measurements, validation, and inference are distinct

- **Measurement:** 4,925 ABL1 variants, Ba/F3 corrected net growth at 1200 nM
  imatinib, with two replicates. Known measurements always take precedence.
- **Saved held-out estimate:** `output/model/validation_predictions.json`
  contains 992 variants across 54 test positions, ESM2 35M / HistGradientBoosting,
  seed 1. MAE 0.01911, Spearman 0.2850, Pearson 0.3124; mean-baseline MAE
  0.02080. Each attached estimate identifies its artifact SHA256, model/head,
  seed, and metrics. The legacy artifact lacks training membership and a
  reference-sequence hash, so those provenance claims cannot be reconstructed.
- **Runtime inference:** `pipeline/predict.py` fits a RandomForest from
  substitution and pocket-position features for unmeasured variants. It does
  not use the ESM embedding head. Its current fit has no attached verified
  validation metrics, and never borrows the ESM metrics. It is limited to
  valid ABL1 substitutions at 242-512 under imatinib.

Net growth combines protein fitness and response under drug. A low value cannot
separate kinase loss of function from drug sensitivity. Replicate range is not
a confidence interval. Pocket distance is a geometric descriptor, not an efficacy
prediction; its reported correlation with measured growth is only 0.069.

Data: Inam et al., *Characterizing Variants of Uncertain Drug Response (VUDRs)
Using Quantitative Measurements at Clinical Exposures*;
[public MIT-licensed source](https://github.com/pritchardlabatpsu/abl_dms).
See `data/datasets/SCHEMA.md` for measurement conventions. Preserve censored
bounds and exclude censored values from exact-error metrics.

## Optional dependencies and sponsor tools

```powershell
python -m venv .venv
.venv\Scripts\python -m pip install torch fair-esm numpy scikit-learn
# Modal GPU scripts only:
.venv\Scripts\python -m pip install modal
```

Copy `.env.example` to `.env` for optional credentials. Never commit keys.
Boltz uses its own `boltz-api` CLI and OAuth session, not Modal tokens.

| Tool | Actual role and integration state |
| --- | --- |
| Modal | GPU embedding script: `scripts/modal_esm.py`; no HTTP backend or deployment yet |
| Boltz | ABL1 kinase crop 242-512 plus imatinib; explicitly enable with `--boltz`; free cost estimate and `VN_BOLTZ_MAX_USD` per-call cap (default $0.25) |
| Ptarmigan | Saved wild-type screens in `output/ptarmigan/`; catalog-only service cannot score mutant sequences; mutant screen explicitly unsupported |
| Paperclip | Separate optional literature module `pipeline/evidence.py`; not automatically invoked or fed into prediction |
| Keenable | Report-agent integration remains unimplemented; local deterministic reports work |
| Nebius / Genomic Intelligence | Not implemented in the workflow |

Boltz availability, budget rejection, timeout, and crop errors are recorded as
skips, allowing reporting to continue. The wrapper is never invoked for EGFR
or variants outside its validated crop. Budget caps are per call, not a total
account spending limit. Live services have not been verified by offline tests. The Boltz panel script
uses the bundled study reference and retries entries previously marked as errors.

## Project layout

- `main.py`: cancer-only orchestration and CLI.
- `pipeline/`: sequence retrieval, variants, assay lookup/modeling, ESM, geometry,
  Boltz, reporting, and separate literature evidence.
- `scripts/`: offline embedding/training, Modal GPU work, Boltz panel, figures.
- `data/`: targets, study reference, assay data, compounds, controls.
- `output/model/`: checked-in legacy validation artifact; new generated output
  remains ignored.
- `tests/`: deterministic service-mocked workflow tests (also pytest-compatible).
- `ui/`: separate prototype; UI work is currently deferred.

## Remaining backend work

- Improve and independently validate the full-dataset ESM model; its current
  held-out performance is modest and it misses resistant T315 substitutions.
  Persistence and provenance are implemented; smoke artifacts are guarded.
- Evaluate the runtime RandomForest reproducibly before attaching accuracy claims.
- Validate Boltz results against measured controls; do not turn confidence into
  drug work/fail calls or apply placeholder affinity thresholds.
- Connect literature evidence after predictions are fixed, and implement an
  optional report rewrite only if needed. Deterministic reporting is supported.
- Add a service/API and deployment when the workflow contract is ready.
- Add more cancer targets and drug conditions only with matching datasets and
  target-specific sequence/structure configuration.

## Modal smoke test and persisted ESM models

Use the existing `.venv` and Modal CLI authentication. GPU dependencies are
installed in the remote image; local training needs NumPy and scikit-learn.
On Windows, set UTF-8 output before invoking Modal:

```powershell
$env:PYTHONUTF8 = '1'
.venv\Scripts\python -m modal run scripts/modal_esm.py --limit 4 --batch 2
.venv\Scripts\python scripts/train_embeddings.py --npz data/esm/P00519_esm2_t12_35M_UR50D_study518_mutdiff_n4.npz --allow-partial --out-dir output/model/smoke --seed 1
```

The Modal default is ESM2 35M with the bundled 518-aa study reference. Limited
runs sample across the dataset and use distinct filenames; they do not replace
full feature artifacts. The script checks residue identity before GPU submission
and checks output dimensions and finiteness before writing.

Training writes `ridge.joblib`, `hist_gbm.joblib`, `manifest.json`, and
`validation_predictions.json` into the chosen output directory. It reloads each
model and verifies identical held-out predictions. The manifest records the
feature recipe, sequence/dataset/feature hashes, reference and crop, model
parameters, software versions, exact split membership, replicate count, and
head-specific metrics. The deployed head is preselected with `--head`, never
chosen by its test performance. Spearman handles tied ranks.

`pipeline.embedding_model.predict_embedding` loads a **trusted local** run and
checks its feature contract, dimensions, checksum, and sklearn version. Joblib
is executable serialization; never load untrusted models. Smoke models require
an explicit `allow_partial=True` override. This function consumes an embedding;
it does not automatically submit GPU jobs. `embedding_model.predict_variant`
validates a trusted run before requesting one Modal embedding. The CLI exposes
this as `--esm-run-dir output/model/abl1-study518-seed1`; the resulting
`esm_prediction` is separate from measured response, with explicit train/test/
unseen position membership. It never labels a training prediction as held out. Four-row smoke metrics have no scientific
validation value. Full training should use a new output directory and omit
`--allow-partial`; keep generated artifacts and credentials out of Git.

### Latest full model run

`output/model/abl1-study518-seed1` contains the new locally generated run:
4,925 variants; 3,933 training rows at 216 positions and 992 test rows at 54
positions; zero shared positions. ESM2 35M, actual bundled study reference
(518 aa), preselected HistGradientBoosting head, seed 1.

| Head | MAE | Spearman | Pearson | Within replicate range |
| --- | --- | --- | --- | --- |
| Ridge | 0.01978 | 0.2168 | 0.2172 | 10.8% |
| HistGradientBoosting | 0.01926 | 0.2681 | 0.2989 | 11.7% |
| Training-mean baseline | 0.02080 | ? | ? | ? |

These are a new run's results, separate from the legacy checked-in artifact.
Do not mix their metrics or predictions. The full manifest records the exact
reference, feature/dataset hashes, package versions, and saved head checksums.
T315I remains a failure case: measured +0.04599 versus held-out estimate
approximately -0.02466. Saving the model does not solve its scientific limits.
Generated full models/features are ignored, not committed; reproduce them with:

```powershell
$env:PYTHONUTF8 = '1'
.venv\Scripts\python -m modal run scripts/modal_esm.py --batch 8
.venv\Scripts\python scripts/train_embeddings.py --npz data/esm/P00519_esm2_t12_35M_UR50D_study518_mutdiff.npz --out-dir output/model/abl1-study518-seed1 --seed 1 --head hist_gbm
.venv\Scripts\python main.py --gene ABL1 --variant T315I --esm-run-dir output/model/abl1-study518-seed1 --no-agent
```

The final command uses online sequence/geometry services and explicitly requests
one Modal model estimate; it keeps measured assay lookup separate. Boltz is
still disabled unless `--boltz` is supplied.
