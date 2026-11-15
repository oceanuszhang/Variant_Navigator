# Variant Navigator: BCR-ABL Growth Direction

## Overall Research Question (Big Picture)

Given an ABL1 amino-acid variant, can frozen protein sequence embeddings predict whether BCR-ABL cell growth increases, decreases, or remains unchanged during a fixed imatinib treatment?

## Specific Targeted Approach

Reconstruct all 4,925 full mutant proteins from the 518-amino-acid study reference. On a Modal T4 GPU, use frozen ESM-2 35M to calculate a 480-dimensional mutant-minus-reference vector at each altered residue. Train Ridge and XGBoost classifiers on those features.

## Input

- Reference: `data/cancer-exp/ABL1_NM_005157.6_DMS_reference_aa.txt`
- Dataset: `data/cancer-exp/high_imatinib_netgrowth_data.csv`
- Treatment: 1,200 nM imatinib for six days
- Experimental measurement: `netgrowth_mean_corrected`
- Replicates: `netgr_corr_rep1.high` and `netgr_corr_rep2.high`

## What Exactly We're Predicting

Three directional classes relative to zero corrected net growth:

- **Decrease:** `netgrowth_mean_corrected < -0.01`
- **Unchanged:** `-0.01 <= netgrowth_mean_corrected <= 0.01`
- **Increase:** `netgrowth_mean_corrected > 0.01`

The neutral band has a total width of 0.02, approximately four times the median replicate interval width. The model no longer needs to reproduce a narrow replicate-specific numerical interval.

## Benchmarking Accuracy

Use one reproducible 80/20 `GroupShuffleSplit` with random seed 42, grouped by `protein_start`. Every substitution at a residue remains entirely in training or validation.

Primary evaluation:

- exact three-class accuracy;
- balanced accuracy, which gives equal weight to decrease, unchanged, and increase;
- macro F1;
- opposite-direction errors.

The majority-class baseline predicts **decrease** for every validation row. Because 70.3% of validation variants decrease, ordinary accuracy alone is misleading. Balanced accuracy and macro F1 are required for comparison.

## Models and Tools

- Frozen `facebook/esm2_t12_35M_UR50D` on a Modal T4 GPU
- Ridge classifier with standardized features, `alpha=1.0`, and balanced class weights
- XGBoost multiclass classifier with 300 trees and balanced training-row weights
- scikit-learn, XGBoost, NumPy, and joblib for local training and evaluation

## Exact Implementation Steps

1. Validate all 4,925 substitutions against the 518-amino-acid reference.
2. Reconstruct every full mutant protein.
3. Embed the reference once and mutant sequences in batches of 32 on Modal.
4. Save the 480-dimensional mutant-minus-reference vector at the altered residue.
5. Cache the aggregate feature matrix, all full sequences, and one embedding file per variant.
6. Convert corrected net growth into decrease, unchanged, or increase using the ±0.01 neutral band.
7. Create the fixed position-grouped train-validation split.
8. Fit the Ridge and XGBoost classifiers with class balancing.
9. Save class predictions, XGBoost class probabilities, confusion matrices, metrics, and trained models.

## Completed Validation Results

Validation class counts: 664 decrease, 148 unchanged, and 132 increase.

| Model | Accuracy | Balanced accuracy | Macro F1 | Opposite-direction errors |
| --- | ---: | ---: | ---: | ---: |
| Majority decrease | 70.34% | 33.33% | 0.275 | 132 (13.98%) |
| Ridge | 41.63% | 38.38% | 0.345 | 191 (20.23%) |
| XGBoost | 65.78% | **42.65%** | **0.434** | 127 (13.45%) |

XGBoost improves class-balanced performance over the majority baseline but still struggles with unchanged and increased-growth variants. Its overall accuracy is lower than the majority baseline because the dataset is strongly imbalanced toward decreased growth.

## Run

```bash
cd /Users/prattmac/Desktop/Variant_Navigator
.venv/bin/python main.py --train-only
```

The cached embeddings avoid another Modal GPU run. Use `.venv/bin/python main.py --recompute-embeddings` only when the sequence features must be regenerated.

## Outputs

- `output/cancer-exp/embeddings.npz`
- `output/cancer-exp/variant_sequences.fasta`
- `output/cancer-exp/variant_embeddings/`
- `output/cancer-exp/variant_manifest.csv`
- `output/cancer-exp/metrics.json`
- `output/cancer-exp/validation_predictions.csv`
- `output/cancer-exp/ridge_model.joblib`
- `output/cancer-exp/xgboost_model.json`
