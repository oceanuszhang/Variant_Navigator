<<<<<<< Updated upstream
# Variant Navigator BCR ABL Imatinib Growth Prediction

This experiment tests whether protein sequence embeddings can predict the quantitative growth of cells carrying different BCR-ABL kinase-domain variants under a fixed imatinib treatment.

## Overall Research Question Big Picture

Can we predict how a protein-coding variant changes treatment response from its effect on the target protein sequence? Here, we test whether ABL1 variants predict the net growth of BCR-ABL-expressing Ba/F3 cells under 1,200 nM imatinib.

## Specific Targeted Approach

Construct each mutant from the study-specific 518-amino-acid ABL1 reference. Run the reference and mutant sequences through frozen ESM-2 35M, extract the embedding at the mutated residue, and calculate the 480-dimensional difference between the mutant and reference embeddings. Use that difference vector to train and compare Ridge and XGBoost regression models.

## Input

- `data/cancer-exp/ABL1_NM_005157.6_DMS_reference_aa.txt`: ABL1 amino-acid reference matching the study's residue numbering.
- `data/cancer-exp/high_imatinib_netgrowth_data.csv`: 4,925 missense variants spanning ABL1 residues 242-512.
- Variant fields: `protein_start`, `ref_aa`, and `alt_aa`.
- Training label: `netgrowth_mean_corrected`.

## What Exactly We Are Predicting

The continuous corrected net-growth rate for each variant under 1,200 nM imatinib. Higher predicted growth indicates greater imatinib resistance; lower predicted growth indicates greater treatment susceptibility.

## Benchmarking Accuracy

Create one reproducible 80/20 training-validation split grouped by `protein_start`, keeping every substitution at a residue entirely within one split. A validation prediction is considered in range when it falls between `netgr_corr_rep1.high` and `netgr_corr_rep2.high` for that variant. Compare Ridge and XGBoost by the percentage of validation predictions in range, with absolute prediction error reported for context.

## Models and Tools

- **ESM-2 35M** (`facebook/esm2_t12_35M_UR50D`): frozen protein representation model.
- **Modal**: remote compute for ESM-2 embedding extraction.
- **Ridge regression**: linear prediction head and primary baseline.
- **XGBoost regression**: nonlinear prediction head.
- **Python, pandas, scikit-learn, and XGBoost**: data preparation, splitting, training, and evaluation.
=======
# Variant Navigator: BCR-ABL Growth Prediction

## Overall Research Question (Big Picture)

Given a protein reference sequence and an amino-acid variant, can a frozen protein language model predict how a cancer cell responds to a treatment?

This experiment tests one well-defined case: BCR-ABL variants in Ba/F3 cells treated with a constant 1,200 nM dose of imatinib for six days.

## Specific Targeted Approach

For each ABL1 variant, construct the mutant protein sequence and use frozen ESM-2 representations to describe how the variant changes its local sequence context. Train two regression heads on those features and predict the experimentally measured corrected net growth.

The active pipeline is intentionally limited to:

1. frozen ESM-2 feature extraction on a Modal GPU;
2. Ridge regression;
3. XGBoost regression;
4. a fixed held-out validation split.

## Input

- Reference protein: `data/cancer-exp/ABL1_NM_005157.6_DMS_reference_aa.txt`
  - 518-amino-acid ABL1 study reference;
  - all measured variant positions are checked against this sequence before GPU work begins.
- Labels: `data/datasets/abl_dms_netgrowth_1200nM.csv`
  - 4,925 amino-acid substitutions;
  - one fixed treatment condition: 1,200 nM imatinib;
  - two corrected replicate measurements per variant;
  - target column: `netgrowth_mean_corrected`.

## What Exactly We're Predicting

The quantitative value of `netgrowth_mean_corrected` for each variant under the fixed imatinib treatment.

- Higher predicted net growth means the variant-bearing cells grew better during treatment.
- Lower predicted net growth means the cells grew worse during treatment.
- This is a continuous regression target, not a resistant/sensitive classification.

Because every row uses the same treatment and dose, the model input is the reference protein plus the amino-acid substitution. It does not yet accept treatment identity or dose as a variable.

## Benchmarking Accuracy

The split is a reproducible 80/20 `GroupShuffleSplit` with random seed 42. Rows are grouped by `protein_start`, so every substitution at a held-out residue is kept out of training.

Primary metric:

- percentage of validation predictions that fall between `netgr_corr_rep1.high` and `netgr_corr_rep2.high`, inclusive.

Context metric:

- mean absolute error against `netgrowth_mean_corrected`.

The same training rows, validation rows, and metrics are used for Ridge and XGBoost.

## Models and Tools

- ESM-2: `facebook/esm2_t12_35M_UR50D`, frozen, 480 dimensions.
- Modal: runs ESM-2 inference on an NVIDIA T4 GPU.
- Ridge: standardized features with `alpha=1.0`.
- XGBoost: squared-error regression with 300 trees, learning rate 0.05, maximum depth 4, minimum child weight 2, row and column subsampling 0.8, L2 regularization 1.0, and seed 42.

## Exact Implementation Steps

1. `main.py` resolves the reference, CSV, and output paths.
2. `pipeline/modal_esm.py` reads the 518-amino-acid sequence and all 4,925 CSV rows.
3. Before requesting a GPU, it verifies each position is in range and that the CSV reference amino acid matches the sequence.
4. Modal builds the ESM-2 environment and starts a T4 GPU function.
5. The reference sequence is embedded once.
6. Each mutant sequence is created by replacing exactly one amino acid; mutants are embedded in batches of 32.
7. The feature for a variant is the 480-dimensional mutant-minus-reference ESM vector at the mutated residue.
8. Features, all reconstructed sequences, variant metadata, target values, and both replicate values are cached in `output/cancer-exp/embeddings.npz`.
9. The pipeline also saves all 4,925 full mutant sequences as FASTA, one 480-dimensional feature file per variant, and a manifest connecting every variant to its files.
10. `pipeline/train_heads.py` creates the fixed position-grouped 80/20 split and asserts that no residue position occurs in both sets.
11. Ridge and XGBoost are fit on the training features and evaluated on the identical validation rows.
12. The pipeline saves metrics, row-level predictions, and both fitted models.

## Run

```bash
cd /Users/prattmac/Desktop/Variant_Navigator
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -m modal setup
.venv/bin/python main.py
```

ESM embeddings are cached. Later runs can retrain only the prediction heads:

```bash
.venv/bin/python main.py --train-only
```

Force a new Modal embedding run with:

```bash
.venv/bin/python main.py --recompute-embeddings
```

On macOS, XGBoost also needs the OpenMP runtime. The standard installation is `brew install libomp`. On this machine, the verified Homebrew Apple Silicon runtime was installed inside `.venv` because the existing Homebrew executable used the Intel prefix.

## Outputs

- `output/cancer-exp/embeddings.npz`: cached frozen ESM features and labels;
- `output/cancer-exp/variant_sequences.fasta`: all 4,925 reconstructed full mutant protein sequences;
- `output/cancer-exp/variant_embeddings/`: one 480-dimensional mutant-minus-reference feature file per variant;
- `output/cancer-exp/variant_manifest.csv`: mapping from every dataset row and variant to its sequence and embedding artifact;
- `output/cancer-exp/metrics.json`: split details and aggregate validation metrics;
- `output/cancer-exp/validation_predictions.csv`: measured values, replicate interval, both predictions, absolute errors, and in-range calls;
- `output/cancer-exp/ridge_model.joblib`: fitted Ridge pipeline;
- `output/cancer-exp/xgboost_model.json`: fitted XGBoost model.

## Current Verification Status

- Local dependencies installed and imports verified.
- All three active Python files compile.
- All 4,925 CSV variants validate against the 518-amino-acid reference.
- Ridge and XGBoost both pass an end-to-end synthetic smoke test, including the grouped split, metrics, predictions, and saved model files.
- The full ESM-2 run completed on a Modal T4 GPU on November 14, 2026.
- The aggregate matrix contains 4,925 rows by 480 features.
- The FASTA contains 4,925 reconstructed sequences, each exactly 518 amino acids long.
- The manifest contains 4,925 data rows and `variant_embeddings/` contains 4,925 individual `.npy` files.
- The grouped split contains 3,981 training rows from 216 positions and 944 validation rows from 54 positions.

## Completed Validation Results

| Model | Predictions within replicate interval | Mean absolute error |
| --- | ---: | ---: |
| Ridge | 88/944 (9.32%) | 0.021912 |
| XGBoost | 114/944 (12.08%) | 0.018927 |

XGBoost performed better on both reported validation metrics. These are held-out-position results: the model never saw another substitution from any validation residue during training.
>>>>>>> Stashed changes
