# VariantNavigator: ABL1 Growth Direction

## Team 

Nicholas Pratt - University of California, Merced; Mohamed bin Zayed University of Artificial Intelligence (MBZUAI)
Oceanus Zhang - University of California, Merced
James Song - Emory University
Franziska Tietze - University of California, Berkeley; Mohamed bin Zayed University of Artificial Intelligence (MBZUAI)

## Overall Research Question (Big Picture)

We implement models to predict how cancel cells with mutant variants respond to cancer drug treatment. Given the mutated protein and the name of its mutation we predict if cell growth increases decreases or stays constant, this tells us how the cells respond to drug treatment. More cell growth indicates the drug treatment was not effective.

Given an ABL1 amino-acid variant, can frozen protein-sequence embeddings predict whether BCR-ABL cell growth decreases, remains unchanged, or increases during fixed-dose imatinib treatment?

## Demo

Please download the short mp4 demo video

<video src="assets/variantnavigator_demo.mp4" controls width="100%"></video>

[▶ Watch the 55-second VariantNavigator demo](assets/variantnavigator_demo.mp4)

## Specific Targeted Approach

Reconstruct 4,925 full mutant proteins from the 518-amino-acid study reference. Use frozen ESM-2 35M to represent each mutation as a 480-dimensional mutant-minus-reference vector at the altered residue. Train Ridge and XGBoost direction classifiers.

Two model copies are deliberately kept separate:

- **Evaluation models:** trained on 216 residue positions and evaluated on 54 entirely held-out positions.
- **Deployment models:** refit on all 4,925 rows to power the live dataset demo. Their per-row predictions are not independent accuracy estimates.

## Input

- Reference: `data/cancer-exp/ABL1_NM_005157.6_DMS_reference_aa.txt`
- Dataset: `data/cancer-exp/high_imatinib_netgrowth_data.csv`
- Treatment: 1,200 nM imatinib for six days
- Target measurement: `netgrowth_mean_corrected`
- Replicates: `netgr_corr_rep1.high` and `netgr_corr_rep2.high`
- Demo input: the exact study reference sequence, one measured substitution, and a Ridge or XGBoost head

## What Exactly We're Predicting

- **Decrease:** `netgrowth_mean_corrected < -0.01`
- **Unchanged:** `-0.01 <= netgrowth_mean_corrected <= 0.01`
- **Increase:** `netgrowth_mean_corrected > 0.01`

This is directional response under one fixed experimental treatment. It is not a clinical efficacy prediction.

## Benchmarking Accuracy

Use one reproducible 80/20 `GroupShuffleSplit`, seed 42, grouped by `protein_start`. Every substitution at a residue remains entirely in training or validation.

| Model | Accuracy | Balanced accuracy | Macro F1 | Opposite-direction errors |
| --- | ---: | ---: | ---: | ---: |
| Majority decrease | 70.34% | 33.33% | 0.275 | 132 / 944 |
| Ridge | 41.63% | 38.38% | 0.345 | 191 / 944 |
| XGBoost | 65.78% | **42.65%** | **0.434** | 127 / 944 |

Validation contains 944 variants at 54 held-out sites: 664 decrease, 148 unchanged, and 132 increase. Balanced accuracy is primary because the assay is strongly imbalanced toward decreased growth.

The XGBoost class recalls are 84.19% decrease, 26.35% unchanged, and 17.42% increase. The model beats the majority baseline on balanced accuracy but remains weak on minority directions.

## Models and Tools

- **ESM-2 35M on Modal T4:** frozen sequence feature extraction. GPU work is only needed to recompute embeddings.
- **Ridge and XGBoost:** local three-class prediction heads.
- **Ptarmigan-1 cached control:** ranks compounds for wild-type ABL1. It cannot score submitted mutant sequences, so the UI labels it as a wild-type control.
- **Boltz-2.1:** optional mutant–imatinib structure/binding follow-up through the `boltz-api` CLI. It uses a 242–512 kinase-domain crop, five samples, a free cost preflight, a $0.25 per-call cap, and disk caching. It does not use Modal.
- **Local web UI:** accepts the exact 518-aa reference plus one of the 4,925 measured substitutions and displays deployment output, measured response, and a separate held-out result when available.

## Exact Implementation Steps

### Completed preprocessing and ESM extraction

1. Validate all 4,925 substitutions against the 518-aa study reference.
2. Reconstruct and save all 4,925 full mutant sequences.
3. Embed the reference once with frozen `facebook/esm2_t12_35M_UR50D`.
4. Embed mutant sequences in Modal GPU batches.
5. For each variant, save the 480-dimensional mutant-minus-reference vector at the altered residue.
6. Save the aggregate `embeddings.npz`, FASTA file, manifest, and one `.npy` file per variant.

### Completed local evaluation and deployment rerun

7. Load the saved 4,925 × 480 feature matrix; do not call Modal.
8. Convert corrected net growth to the three classes using the ±0.01 neutral band.
9. Recreate the site-grouped split with seed 42 and verify zero position overlap.
10. Fit balanced Ridge and XGBoost evaluation heads on the training positions.
11. Save held-out predictions, XGBoost probabilities, class recalls, confusion matrices, and metrics under `output/cancer-exp/evaluation/`.
12. Refit separate Ridge and XGBoost deployment heads on all 4,925 rows.
13. Save all-row predictions and a manifest containing sequence, embedding, and model SHA-256 hashes under `output/cancer-exp/deployment/`.
14. Preserve the existing root-level metrics and prediction CSV for backward compatibility.

### Completed main-branch integration

15. Keep the NickP ESM feature definition, split, Ridge logic, and XGBoost logic.
16. Add the main-branch Boltz-2.1 adapter, compound records, and cost/caching safeguards.
17. Add the main-branch cached Ptarmigan ABL1 and EGFR wild-type outputs; do not represent them as mutant predictions.
18. Replace the old RandomForest/regression UI backend with the NickP three-class artifact store.
19. Add UI routes for the reference, prediction, validation, Ptarmigan control, free Boltz estimate, and explicitly confirmed Boltz run.
20. Make the demo distinguish all-data deployment output from independently held-out evaluation output.
21. Verify `E450R`: deployment XGBoost = increase, measured = increase, held-out XGBoost = increase.
22. Run local integration tests and browser/API smoke tests.

No Modal GPU job or paid Boltz prediction was submitted during this integration.

## Run

Reuse saved embeddings and regenerate all evaluation/deployment artifacts:

```bash
cd /Users/prattmac/Desktop/Variant_Navigator
.venv/bin/python main.py --train-only
```

Start the demo:

```bash
cd /Users/prattmac/Desktop/Variant_Navigator
.venv/bin/python ui/server.py --port 8000
```

Open `http://127.0.0.1:8000` and use `E450R` for the clean held-out demonstration.

Run integration checks:

```bash
.venv/bin/python -m unittest -v tests/test_integrated_demo.py
```

Only recompute ESM embeddings when required:

```bash
.venv/bin/python main.py --recompute-embeddings
```

That last command submits GPU work through Modal.

## Outputs

- Saved features: `output/cancer-exp/embeddings.npz`
- Full sequences: `output/cancer-exp/variant_sequences.fasta`
- Per-variant embeddings: `output/cancer-exp/variant_embeddings/`
- Variant manifest: `output/cancer-exp/variant_manifest.csv`
- Held-out evaluation: `output/cancer-exp/evaluation/`
- Full-data demo models and predictions: `output/cancer-exp/deployment/`
- Cached Ptarmigan controls: `output/ptarmigan/`
- Optional Boltz cache: `output/structures/`
