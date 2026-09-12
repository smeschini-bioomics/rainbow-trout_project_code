# Midgut histone post-translational modifications

## Purpose

This module processes midgut histone H3/H4 LC–MS/MS measurements acquired by data-dependent acquisition (DDA), from peptidoform-level peptide-ion intensities to protein-adjusted histone post-translational modification (hPTM) features. It then evaluates diet-by-day effects with feature-specific Bayesian regression, assesses technical batch structure, generates batch-corrected profile-level visualisations, and analyses diet-specific hPTM co-variation.

## Histone derivatisation method

Midgut histones were processed by propionylation-based chemical derivatisation before DDA LC–MS/MS analysis, following the histone derivatisation strategy described by Garcia et al. (2007). The hPTM DDA-MS workflow and interpretation were additionally informed by the histone PTM atlas described by Provez et al. (2022).

> Garcia BA, Mollah S, Ueberheide BM, Busby SA, Muratore TL, Shabanowitz J, Hunt DF. Chemical derivatization of histones for facilitated analysis by mass spectrometry. *Nature Protocols*. 2007;2:933–938. doi:10.1038/nprot.2007.106.

> Provez L, Van Puyvelde B, Corveleyn L, et al. An interactive mass spectrometry atlas of histone posttranslational modifications in T-cell acute leukemia. *Scientific Data*. 2022;9:1. doi:10.1038/s41597-022-01736-1.

## Sample preparation and LC–MS/MS analysis

Histone samples were processed and analysed by [ProGenTomics](https://www.progentomics.ugent.be/), Ghent, Belgium. ProGenTomics provided the histone sample-preparation workflow and DDA LC–MS/MS acquisition. Downstream computational processing and statistical analyses are documented in this repository.

## Analytical branches

The workflow has two complementary analytical branches:

- **Feature-specific differential analysis** uses the processed, uncorrected PTM matrix. `gut_hptm_batch` is included as a fixed covariate in each Bayesian model. All samples were processed within the same overall workflow, but centrifugation was performed sequentially in groups of 30 samples; `gut_hptm_batch` represents this technical centrifugation batch.

- **Profile-level visualisation and co-variation analysis** use the limma batch-corrected matrix. Batch correction protects the biological `diet * day_categ` design, so diet, day, and diet-by-day variation are retained.

## Directory layout

## Directory layout

- `01_input_files/` — curated H3/H4 peptide-ion table derived from the Progenesis export, containing raw non-convoluted abundances and formatted for downstream processing with `QFeatures` and `msqrob2`.
- `01_output_files/` — protein-adjusted peptidoform and PTM-level matrices.
- `02_process_matrix/` — processed sample-by-PTM matrix used by the Bayesian models and batch workflow.
- `03_batch_effect_pca_exploration/` — initial batch-structure PCA panels and multivariate batch-assessment report.
- `04_brms_regression/` — candidate-family diagnostics, selected-family table, fitted feature-level Bayesian models, posterior contrasts, and diet-effect dot plot.
- `05_batch_correction_limma/` — limma-corrected matrix, before/after correction diagnostics, and batch-balance outputs.
- `06_samples_features_heatmap_zscored/` — corrected-matrix heatmap, annotation tables, and z-score export.
- `07_PCA_PERMANOVA/` — PCA, PERMANOVA, dispersion tests, and PC-score comparisons.
- `08_a_covariance analysis_corr/` — diet-specific hPTM Spearman correlation matrices.
- `08_b_differential_covariance_analysis/` — supplementary variability and differential-correlation analyses.
- `progenesis_final_output/` — CSV exported from Progenesis QI containing chromatographic and mass-spectrometric information of non-deconvoluted peptide ion raw abundances. At this stage, protein assignments were manually reviewed, peptide identifications were assigned to the corresponding histone family, and ambiguous or non-specific peptides were excluded. Only peptides assigned to histones H3 or H4 were retained. The curated table was then reorganised into the layout required for the `QFeatures`/`msqrob2` workflow and used to generate the files in `01_input_files/` for downstream analysis.

Run scripts from `05_Midgut_hPTMs/` so their relative paths resolve correctly.

## Inputs

### Peptide-ion input

`01_msqrob2PTM_filtering_normalization.R` expects:

```text
01_input_files/20260420_TroutGut_peptides_ions_raw_ab_rank3_msqrob-input.csv
```

The peptide-ion dataset contains samples acquired across two independent MS run dates. Exploratory batch assessment indicated a substantial run-date effect. To avoid carrying this technical variation into downstream hPTM analyses, the workflow retains only samples whose quantitative-column names contain `240816`.

Samples from the other acquisition-date subset are excluded before peptidoform filtering, normalisation, and PTM summarisation.

Downstream processing follows the `msqrob2PTM` framework for post-translational modification and peptidoform-level proteomics analysis [Demeulemeester et al., 2023](https://doi.org/10.1016/j.mcpro.2023.100708).


### Metadata

The module uses:

```text
STPN2309_metadata_all.tsv
```

The metadata file must contain the project sample identifier as its first column and the following analysis-specific variables:

| Use | Required metadata field(s) |
|---|---|
| Map raw hPTM sample labels to project sample IDs | `sampleID_hPTMs_gut` |
| Differential models | `diet`, `day_categ`, `gut_hptm_batch` |
| Initial batch exploration | `sample_id`, `diet`, `day_categ`, `rearing_tank`, `gut_hptm_batch` |
| limma correction and downstream analyses | `sampleID_STPN2309`, `diet`, `day_categ`, `gut_hptm_batch` |
| Additional PCA groupings | `phase`, `diet_phase`, `day_categ` |

The primary analyses retain `NC` and `HC` samples. `NC` is the reference diet, and the ordered sampling days are `day_01`, `day_02`, `day_03`, `day_04`, `day_10`, `day_15`, and `day_22`.

## Script order

1. `01_msqrob2PTM_filtering_normalization.R` — peptidoform preprocessing, protein adjustment, robust PTM summarisation, and matrix export.
2. `02_process_matrix.R` — PTM feature-name standardisation, duplicate-signal collapsing, missingness filtering, and sample-ID harmonisation.
3. `03_batch_effect_pca_exploration.R` — initial PCA and technical-batch assessment before correction.
4. `04_a-brms_regression.R` — Bayesian family comparison and final feature-specific diet-by-day regression.
5. `04_b_brms_regression_heatmap_results.R` — posterior diet-effect dot plot/heatmap. Update `base_dir` when a new timestamped Bayesian-model run is generated.
6. `05_batch_correction_limma.R` — batch correction for downstream global-profile analyses, with biological terms retained.
7. `06_samples_features_heatmap_zscored.R` — clustered heatmap of corrected, feature-wise z-scored hPTM profiles.
8. `07_PCA_PERMANOVA.R` — PCA, PERMANOVA, dispersion tests, PC-score comparisons, and analogous analyses for temporal metadata groupings.
9. `08_a_covariance analysis_corr.R` — NC- and HC-specific Spearman hPTM correlation matrices.
10. `08_b_covariance_supp_analysis_Cutler_style.R` — supplementary co-variation, differential-correlation, and variability analyses.


## 1. Peptidoform preprocessing and PTM summarisation

`01_msqrob2PTM_filtering_normalization.R` uses `QFeatures`, `msqrob2`, and `MsCoreUtils` to process H3/H4 peptidoform intensities.

### Preprocessing

The script:

1. imports the selected peptide-ion intensity columns;
2. creates a peptidoform identifier from peptide sequence and variable-modification annotation;
3. treats zero intensities as missing values rather than quantified zero abundance;
4. retains peptidoforms observed with non-zero intensity in at least 50% of retained samples;
5. log2-transforms intensity values;
6. median-centres samples.

### Protein adjustment and PTM-level summarisation

For each sample, retained H3 and H4 peptidoforms are robustly summarised to estimate parent-histone abundance. This summary is subtracted from each log2 peptidoform intensity, producing a protein-adjusted log2 ratio that represents relative peptidoform usage within H3 or H4.

For each histone protein and modification-site label, all protein-adjusted peptidoforms containing that PTM are then robustly aggregated within sample to create one PTM-level feature. Thus, the final matrix represents relative PTM usage conditional on parent H3/H4 abundance, rather than absolute modified-peptide intensity or PTM occupancy.

The script writes:

```text
01_output_files/TroutGut_peptidoform_H3H4.csv
01_output_files/TroutGut_hPTM_H3H4.csv
```

The resulting values are protein-adjusted, log2-scale relative intensities; they are not absolute histone occupancies.

## 2. PTM matrix processing

`02_process_matrix.R` converts the exported PTM matrix to a sample-by-feature matrix aligned to study metadata.

The script:

1. converts verbose modification labels to compact PTM names, for example `H3K9Ac`;
2. retains H3/H4 PTMs containing acetylation (`Ac`), dimethylation (`Me2`), trimethylation (`Me3`), crotonylation (`Cr`), or lactylation (`La`);
3. collapses PTM columns with identical intensity profiles and joins their names with `|`;
4. removes PTM columns with more than 30% missing values;
5. maps assay sample labels through `sampleID_hPTMs_gut` to the project sample IDs.

The processed matrix is saved as:

```text
02_process_matrix/H3_H4_5PTM_midgut.rds
```

## 3. Initial technical-batch assessment

`03_batch_effect_pca_exploration.R` evaluates technical structure before batch correction.

For PCA visualisation only, PTMs with fewer than two observed values or near-zero variance are removed, remaining missing values are median-imputed feature-wise, and PTMs are centred and scaled. The script creates PC1/PC2 panels with marginal density plots coloured by:

- diet;
- sampling day;
- rearing tank;
- `gut_hptm_batch`.

It then performs a complementary multivariate batch assessment on Euclidean distances calculated from feature-wise standardised PTM values:

- unadjusted batch PERMANOVA;
- batch dispersion testing with `betadisper` and a 9,999-permutation test;
- marginal PERMANOVA for `diet + day_categ + gut_hptm_batch`;
- permutation chi-square tests for batch-by-day and batch-by-`Diet × Day` allocation;
- a rank check for the model matrix `~ gut_hptm_batch + diet * day_categ`.

The script writes a text report explaining whether the biological and batch terms are estimable without exact confounding. This stage is diagnostic only and does not create a corrected matrix.

## 4. Feature-specific Bayesian diet-by-day models

`04_a-brms_regression.R` fits one Bayesian model per PTM feature using `brms`.

### Model formula

For every PTM and every candidate family, the model is:

```r
response ~ day_categ * diet + gut_hptm_batch
```

The response is the protein-adjusted log2 PTM intensity. `NC` and `day_01` are the reference levels. The batch variable is included as a fixed covariate in both family comparison and final feature-specific models.

Each PTM is modelled using its available non-missing observations. Features with fewer than six non-missing measurements, or with fewer than two represented batch levels after filtering, are skipped.

### Candidate families and priors

The script compares:

- Gaussian;
- Student-t;
- skew-normal.

The common priors are:

```r
b         ~ normal(0, 2)
Intercept ~ student_t(3, 0, 5)
sigma     ~ exponential(1)
```

For Student-t models, the degrees-of-freedom parameter additionally uses:

```r
nu ~ gamma(2, 0.2)
```

### Family comparison and final fits

Candidate models use four chains, 1,000 iterations per chain, and 500 warm-up iterations. Each fit is evaluated using LOOIC, Bayesian R², R-hat, effective-sample-size ratio, and the number of observed values.

Final selected-family models use four chains, 4,000 iterations per chain, 2,000 warm-up iterations. Per-feature outputs include model summaries, posterior predictive checks, fixed-effect trace plots, conditional-effect plots, estimated marginal means, and LOO summaries.

### Posterior HC-minus-NC contrasts

For each sampling day, `emmeans` obtains the posterior HC-minus-NC contrast. Posterior draws are summarised as:

- posterior mean (`estimate`);
- 95% credible interval (`lower_95`, `upper_95`);
- posterior probability of a positive contrast (`prob_gt0`).

The consolidated table is:

```text
04_brms_regression/<timestamped_run>/summary_tables/
all_features_contrast_diet_by_day_best_family_posterior.csv
```

These effects are HC-minus-NC differences on the protein-adjusted log2-intensity scale. The model is Bayesian; no frequentist multiple-testing-adjusted p-value is produced by this script.

### Diet-effect dot plot

`04_b_brms_regression_heatmap_results.R` filters the posterior table to `HC - NC` contrasts, clusters PTMs according to their day-specific posterior estimates, and adds a dendrogram.

In the final dot plot:

- fill colour indicates direction (`HC > NC` or `HC < NC`);
- dot size is `|HC − NC|` on the log2 scale;
- green outline indicates a 95% credible interval excluding zero;
- text labels show probability of direction (`Pd`) only when `Pd > 0.95`.

The script has a hard-coded `base_dir` pointing to a timestamped final-model run. Update this path after rerunning script 04a.

## 5. limma batch correction for profile-level analyses

`05_batch_correction_limma.R` is used for global-profile visualisation, multivariate analyses, and co-variation analyses. It is not the matrix used for the feature-specific Bayesian contrast estimates.

The script removes PTMs with fewer than two observed values or near-zero global variance, checks batch balance by diet and day, then uses:

```r
limma::removeBatchEffect(
  x = PTM_matrix,
  batch = gut_hptm_batch,
  design = model.matrix(~ diet * day_categ)
)
```

Thus, the correction removes variation associated with `gut_hptm_batch` while retaining the specified diet, day, and diet-by-day biological design.

Batch-correction QC includes:

- PCA panels before and after correction;
- PC1/PC2 association tests for batch, diet, and day;
- feature-level linear models testing batch effects while adjusting for `diet * day_categ`;
- Benjamini–Hochberg adjusted batch p-values before and after correction;
- correlation of HC-minus-NC PTM differences before versus after correction.

The script saves both the original and corrected matrices, metadata, PCA objects, and QC tables in:

```text
05_batch_correction_limma/H3_H4_5PTM_midgut_batch_corrected_limma.rds
```

## 6. Corrected-matrix heatmap and multivariate structure

### Z-scored heatmap

`06_samples_features_heatmap_zscored.R` uses the limma-corrected matrix. PTMs are z-scored across samples using observed values only. For display, values can be capped at the 99th percentile of absolute z-scores; the chosen cap and number of capped cells are exported.

Samples are clustered from Spearman correlations among z-scored profiles. PTMs are clustered from Pearson correlations among PTM z-score profiles, using complete linkage. The heatmap includes diet, day, temporal phase, and `Diet × Temporal Phase` annotations.

Temporal phase is defined as:

- early: days 1–4;
- late: days 10, 15, and 22.

### PCA, PERMANOVA, and dispersion testing

`07_PCA_PERMANOVA.R` uses the limma-corrected matrix. It z-scores every PTM across samples; missing values are replaced by zero after z-scoring, corresponding to the feature mean, so Euclidean-distance and PCA calculations can be performed.

For diet, the script performs:

- PERMANOVA on Euclidean distances from z-scored PTM profiles, with 9,999 permutations;
- `betadisper` and a 9,999-permutation dispersion test;
- PCA directly on the z-scored matrix;
- PC1 and PC2 Wilcoxon rank-sum tests, with BH adjustment across the two axes.

Additional PCA/PERMANOVA/dispersion analyses are run for:

- `phase`;
- `diet_phase`;
- `day_categ`.

For two-level groupings, PC-score comparisons use Wilcoxon tests. For groupings with more than two levels, the script uses Kruskal–Wallis tests, with BH adjustment across PC1 and PC2.

## 7. PTM co-variation and supplementary differential-correlation analysis

The scripts named `covariance` perform **Spearman correlation-based co-variation analyses**, not raw covariance-matrix estimation.

### Diet-specific correlation matrices

`08_a_covariance analysis_corr.R` uses the limma-corrected PTM matrix and calculates separate Spearman correlation matrices for NC and HC using `Hmisc::rcorr()`.

A shared PTM order is obtained by complete-linkage clustering of the all-sample correlation matrix using `1 − rho` as distance. For each diet, the script exports:

- Spearman correlation matrix;
- nominal p-value matrix;
- pairwise sample-size matrix;
- correlation plot.

The displayed matrix plots retain only correlations with nominal `p < 0.05`. This display threshold is not BH-adjusted and should be interpreted as exploratory.

### Differential correlation and variability

`08_b_covariance_supp_analysis_Cutler_style.R` uses the NC and HC correlation matrices from script 08a and the limma-corrected PTM matrix.

It performs:

1. an IQR-based comparison of PTM variability between diets, using a Wilcoxon rank-sum test across PTM-specific IQR values;

2. a 10,000-permutation Mantel comparison of the complete NC and HC correlation structures;

3. an exploratory distributional summary of correlation shifts using:

   ```text
   log2[(rho_HC + 1) / (rho_NC + 1)]
   ```

4. a pair-specific differential-correlation permutation test with 10,000 random reallocations of samples while preserving the observed NC and HC group sizes;

5. Benjamini–Hochberg adjustment across tested PTM pairs.

The differential-correlation matrix displays the log2 fold change of the shifted Spearman correlation,

```text
log2[(rho_HC + 1) / (rho_NC + 1)]
```

rather than raw Spearman rho values. Positive values indicate correlations higher in HC than NC, whereas negative values indicate correlations lower in HC. Only PTM pairs passing the selected BH-FDR threshold are displayed.

The pair-level output reports `rho_NC`, `rho_HC`, `delta_rho`, the transformed differential-correlation value, permutation p-value, BH FDR, and whether the correlation sign changes between diets.

### Relationship to Cutler et al. (2025)

The differential co-variation framework was adapted conceptually from the hPTM biological-variability analyses in Cutler et al. (2025), *Mass spectrometry-based profiling of single-cell histone post-translational modifications to dissect chromatin heterogeneity*, *Nature Communications* 16, 11100. The associated reproducibility repository is available at [cutleraging/single-cell-histone-ptm](https://github.com/cutleraging/single-cell-histone-ptm).

Specifically, the implementation was informed by the treatment-versus-control analysis in [`7-biological-variability/auto+nabut/biological-variability-auto+nabut.Rmd`](https://github.com/cutleraging/single-cell-histone-ptm/blob/main/7-biological-variability/auto%2Bnabut/biological-variability-auto%2Bnabut.Rmd). This includes the construction of group-specific PTM correlation matrices, comparison of overall correlation structure by a 10,000-permutation Mantel test, calculation of differential correlations as:

```text
log2[(rho_HC + 1) / (rho_NC + 1)]
```

and assessment of pair-specific differential correlations by 10,000 random reallocations of samples followed by Benjamini–Hochberg correction. The differential-correlation matrix displays this transformed log2 fold change rather than raw correlation coefficients.

The present implementation differs substantially in the analytical input and biological design. hPTM features were generated from bulk midgut or liver samples using MSqRobPTM processing, including normalisation, log2 transformation, and robust protein-level summarisation to adjust PTM abundance for the corresponding protein abundance. Correlations were then calculated on protein-adjusted, limma-corrected PTM-level features using Spearman correlation. In contrast, the source workflow analyses single-cell peptidoform-level MS1 ratio/relative-abundance data and calculates Pearson correlations.

The IQR-based PTM-variability comparison used here is an additional adaptation: it compares PTM-specific IQR values between diets using a Wilcoxon rank-sum test, whereas the source workflow evaluates CV differences using a t-test. Thus, the source repository provides the conceptual basis for the differential co-variation framework, but this workflow is not an identical reproduction of the original analysis pipeline.


## Main outputs

- `01_output_files/TroutGut_peptidoform_H3H4.csv` — protein-adjusted peptidoform matrix.
- `01_output_files/TroutGut_hPTM_H3H4.csv` — robustly summarised H3/H4 PTM matrix.
- `02_process_matrix/H3_H4_5PTM_midgut.rds` — processed matrix used for Bayesian modelling and batch correction.
- `03_batch_effect_pca_exploration/` — PCA panels and batch-assessment report.
- `04_brms_regression/` — candidate-family diagnostics, selected-family tables, per-feature models, posterior contrasts, and dot plot.
- `05_batch_correction_limma/` — limma-corrected matrix and before/after QC.
- `06_samples_features_heatmap_zscored/` — corrected-matrix z-score heatmap and annotation tables.
- `07_PCA_PERMANOVA/` — diet, phase, diet-phase, and day multivariate outputs.
- `08_a_covariance analysis_corr/` — NC/HC correlation, p-value, and pairwise-`n` matrices.
- `08_b_differential_covariance_analysis/` — IQR, Mantel, differential-correlation, and permutation-analysis outputs.

## Main R packages

Core packages are `QFeatures`, `msqrob2`, `MsCoreUtils`, `SummarizedExperiment`, `brms`, `loo`, `emmeans`, `tidybayes`, `bayesplot`, `limma`, `mixOmics`, `vegan`, `Hmisc`, `corrplot`, `ComplexHeatmap`, `circlize`, `ggtree`, `ggplot2`, `dplyr`, `tidyr`, `purrr`, `readr`, and `data.table`.
