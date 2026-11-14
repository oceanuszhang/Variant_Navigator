# Measured response data — required format

The regression has no training signal until these two files exist. Both are
**transcribed from published figures/tables**; neither is fetchable by API.

## `egfr_ic50.csv` — Experiment 1

Source: EGFR Ba/F3 screen, Figure 2 source data.

| column | type | meaning |
| --- | --- | --- |
| `genotype` | str | full genotype as published, e.g. `L858R`, `L858R/T790M` |
| `sites` | str | semicolon-separated substitutions, for site-level holdout |
| `compound` | str | `erlotinib` or `osimertinib` |
| `ic50_nm` | float | measured IC50, nM |
| `ic50_ref_nm` | float | wild-type reference IC50 in the same assay/plate |
| `log2_fc` | float | `log2(ic50_nm / ic50_ref_nm)` — the regression target |
| `censored` | bool | true where 50% inhibition was never reached |
| `bound` | str | `>` or `<` when censored, else empty |

## `rice_gr50.csv` — Experiment 2

Source: weedy rice ALS study, Figure 4 + Supplementary Table S1.

| column | type | meaning |
| --- | --- | --- |
| `accession` | str | weedy rice accession id |
| `genotype` | str | ALS substitutions in **Arabidopsis numbering** (as published) |
| `sites` | str | semicolon-separated; `S653` substitutions are held out together |
| `compound` | str | `imazethapyr` or `nicosulfuron` |
| `gr50` | float | measured GR50 |
| `gr50_ref` | float | IRGA417 reference GR50 |
| `log2_fc` | float | `log2(gr50 / gr50_ref)` — the regression target |
| `resistance_index` | float | as published, for cross-checking the transcription |
| `censored` | bool | true where 50% reduction was never reached |

## Rules

- **Censored rows keep their bound and are excluded from exact-error metrics.**
  Dropping them silently biases the error downward; recording them as a number
  biases it upward.
- `genotype` stays in published numbering. Renumbering to UniProt happens in
  code via `data/targets.json:numbering_offset`, so the transcription stays
  checkable against the paper.
- Record the reference potency per compound, not just the ratio — a
  cross-compound potency comparison needs it.
