# Liver histone post-translational modifications

## Purpose

This module analyses liver histone H3/H4 post-translational modifications (hPTMs) measured by data-dependent acquisition (DDA) LC–MS/MS. It produces a protein-adjusted PTM-level matrix, evaluates diet-by-day effects with feature-specific Bayesian models, assesses and corrects technical batch variation for profile-level analyses, and examines diet-specific hPTM co-variation.

## Sample preparation and LC–MS/MS analysis
Liver histones were derivatised with trimethylacetic anhydride (TMA) before DDA LC–MS/MS analysis (Kuchaříková et al., 2021).

> Kuchaříková H, Dobrovolná P, Lochmanová G, Zdráhal Z. Trimethylacetic Anhydride-Based Derivatization Facilitates Quantification of Histone Marks at the MS1 Level. *Molecular & Cellular Proteomics*. 2021;20:100114. doi:10.1016/j.mcpro.2021.100114.

Histone samples were processed and analysed by [ProGenTomics](https://www.progentomics.ugent.be/), Ghent, Belgium. ProGenTomics provided the histone sample-preparation workflow and DDA LC–MS/MS acquisition. Downstream computational processing and statistical analyses are documented in this repository.

## Shared peptidoform-to-PTM workflow

Liver data use **exactly the same peptidoform preprocessing, parent-histone adjustment, and PTM inference workflow** as the midgut hPTM dataset. See the section **“Peptidoform preprocessing and PTM summarisation”** in the [Midgut hPTM README](../05_Midgut_hPTMs/README.md#1-peptidoform-preprocessing-and-ptm-summarisation).

Briefly, DDA peptide-ion intensities are processed as:

```text
peptide-ion intensities
→ peptidoforms
→ robust H3/H4 parent-protein summaries
→ protein-adjusted peptidoform intensities
→ robust aggregation of peptidoforms carrying each PTM
→ one protein-adjusted PTM-level value per sample
```


## Script order

1. `01_msqrob2PTM_filtering_normalization.R` — shared DDA peptidoform-to-PTM workflow and liver matrix export.
2. `02_process_matrix.R` — liver-specific PTM retention, missingness filtering, and sample-ID harmonisation.
3. `03_batch_effect_pca_exploration.R` — technical-batch assessment before correction.
4. `04_a-brms_regression.R` — feature-specific Bayesian family comparison and final batch-adjusted diet-by-day models.
5. `04_b_brms_regression_heatmap_results.R` — HC-minus-NC posterior dot plot; update `base_dir` after a new timestamped model run.
6. `05_batch_correction_limma.R` — limma batch correction for profile-level analyses.
7. `06_samples_features_heatmap_zscored.R` — corrected-matrix heatmap.
8. `07_PCA_PERMANOVA.R` — PCA, PERMANOVA, dispersion tests, and PC-score comparisons.
9. `08_a_covariance analysis_corr.R` — NC- and HC-specific Spearman correlation matrices.
10. `08_b_covariance_supp_analysis_Cutler_style.R` — supplementary variability and differential-correlation analyses.

## Technical batch handling

The liver technical batch variable is `trout_liver_batch`.

`03_batch_effect_pca_exploration.R` assesses its relationship with hPTM profiles through PCA panels, Euclidean-distance PERMANOVA, multivariate dispersion testing, marginal PERMANOVA for diet, day, and batch, batch-allocation tests, and a rank check for the batch-adjusted biological design.

Two distinct batch strategies are then used:

- **Feature-specific differential analysis:** `04_a-brms_regression.R` uses the original processed PTM matrix and includes technical batch as a fixed covariate:

  ```r
  response ~ day_categ * diet + liver_hptm_batch
  ```

  Here, `liver_hptm_batch` is the model variable created from `trout_liver_batch`.

- **Profile-level analyses:** `05_batch_correction_limma.R` uses `limma::removeBatchEffect()` to remove `trout_liver_batch` variation while protecting the biological design:

  ```r
  design = model.matrix(~ diet * day_categ)
  ```

  The corrected matrix is used only for heatmaps, PCA/PERMANOVA, and co-variation analyses, not for the feature-specific Bayesian contrasts.

The corrected matrix and associated QC outputs are saved under:

```text
05_batch_correction_limma/
```

including:

```text
H3_H4_8PTM_liver_batch_corrected_limma.rds
```

## Bayesian differential analysis

`04_a-brms_regression.R` follows the same Bayesian workflow as the Midgut hPTM module.

For each retained liver PTM, Gaussian, Student-t, and skew-normal response families are compared. Family selection is based on LOOIC, its uncertainty, convergence metrics, and Bayesian R².

Candidate families were compared for each retained PTM, but the Student-t family was retained for all features in both the liver and midgut final analyses. The final model was therefore:

```r
response ~ day_categ * diet + liver_hptm_batch
family = student()
```

Posterior HC-minus-NC contrasts are calculated separately for each sampling day. Effects are differences on the protein-adjusted log2-intensity scale.

Outputs are written to timestamped folders within:
```text
04_brms_regression/
```

and include selected-family tables, fitted models, posterior predictive checks, estimated marginal means, HC-minus-NC contrasts, LOO summaries, and the posterior diet-effect dot plot.

## Downstream profile and co-variation analyses

The limma-corrected liver matrix is used for:

- `06_samples_features_heatmap_zscored.R` — feature-wise z-scored hPTM heatmap with diet, day, temporal-phase, and diet-by-phase annotations;
- `07_PCA_PERMANOVA.R` — PCA, Euclidean PERMANOVA, dispersion tests, and PC-score comparisons for diet, phase, diet-by-phase, and sampling day;
- `08_a_covariance analysis_corr.R` — separate NC and HC Spearman hPTM correlation matrices;
- `08_b_covariance_supp_analysis_Cutler_style.R` — PTM variability, global correlation-structure comparison, and pair-specific differential-correlation permutation testing.

The covariance-named scripts evaluate **PTM-level Spearman co-variation**, not peptidoform-level covariance.

The supplementary co-variation analysis is conceptually informed by the hPTM covariance/network analyses in Cutler et al. (2025); see the [Cutler et al. reproducibility repository](https://github.com/cutleraging/single-cell-histone-ptm). The present implementation differs by analysing protein-adjusted, batch-corrected PTM features from bulk liver samples rather than single-cell or peptidoform-level data.

## Main outputs

- `01_output_files/TroutLiver_peptidoform_H3H4.csv` — protein-adjusted peptidoform matrix.
- `01_output_files/TroutLiver_hPTM_H3H4.csv` — PTM-level matrix before liver-specific feature retention.
- `02_process_matrix/H3_H4_8PTM_liver.rds` — final processed liver hPTM matrix.
- `03_batch_effect_pca_exploration/` — pre-correction liver batch-assessment outputs.
- `04_brms_regression/` — selected-family Bayesian models and posterior HC-minus-NC results.
- `05_batch_correction_limma/` — corrected matrix and before/after correction QC.
- `06_samples_features_heatmap_zscored/` — corrected-matrix heatmap outputs.
- `07_PCA_PERMANOVA/` — multivariate outputs.
- `08_a_covariance analysis_corr/` and `08_b_differential_covariance_analysis/` — hPTM co-variation and differential-correlation outputs.

## Main R packages

Core packages include `QFeatures`, `msqrob2`, `MsCoreUtils`, `SummarizedExperiment`, `brms`, `loo`, `emmeans`, `tidybayes`, `bayesplot`, `limma`, `vegan`, `Hmisc`, `corrplot`, `ComplexHeatmap`, `ggtree`, `ggplot2`, `dplyr`, `tidyr`, `purrr`, `readr`, and `data.table`.
