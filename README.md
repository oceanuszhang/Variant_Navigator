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
