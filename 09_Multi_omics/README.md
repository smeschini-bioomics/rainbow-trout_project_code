# Multi-omics integration

## Purpose

This module integrates six molecular datasets from independent rainbow trout sampled after exposure to the NC or HC diet at days 1, 2, 3, 4, 10, 15, and 22.

The integration method is **MOFA2** (Multi-Omics Factor Analysis). MOFA learns latent factors: continuous sample-level axes that capture coordinated variation across multiple omics layers. Each factor is interpreted jointly from:

- its sample scores;
- variance explained in each omics view;
- feature weights within each contributing view;
- associations with experimental variables and host traits.

The module also contains a retrospective **MultiGroupPower** analysis based on the six processed omics matrices.

Fish sampled at different time points are independent individuals. The integration does not treat the experiment as longitudinal.

## Integrated omics views

The fitted MOFA model contains six views:

## Input data preparation for MOFA2

Each molecular layer is processed and quality-controlled in its corresponding analysis module before multi-omics integration. The MOFA2 preparation script then aligns samples across views, applies the proteome variability filter, and feature-wise z-scores all retained matrices before model training.

| View                 | Samples | Features | Normalisation / transformation                         | Feature filtering                                 | Batch correction   | Additional MOFA2 preparation                           |
| -------------------- | ---------------: | -------: | ------------------------------------------------------ | ------------------------------------------------- | ------------------ | ------------------------------------------------------ |
| `Midgut hPTMs`       |          153/168 |       20 | Log2-transformed, median-centred, and protein-adjusted | Retained PTMs with ≤30% missing values            | Yes (`limma`)      | Feature-wise z-scoring                                 |
| `Liver hPTMs`        |          168/168 |       28 | Log2-transformed, median-centred, and protein-adjusted | Retained PTMs with ≤30% missing values            | Yes (`limma`)      | Feature-wise z-scoring                                 |
| `Midgut Proteome`    |          151/168 |    3,412 | MS-DAP variance-stabilising normalisation (VSN)        | Retained proteins with ≤30% missing values        | Yes (`limma`)      | Proteins with SD < 0.3 removed; feature-wise z-scoring |
| `Liver Proteome`     |          157/168 |    2,841 | MS-DAP variance-stabilising normalisation (VSN)        | Retained proteins with ≤30% missing values        | Yes (`limma`)      | Proteins with SD < 0.3 removed; feature-wise z-scoring |
| `Digesta Microbiota` | 165/168 | 140 | CLR transformation | ASV-level count matrix filtered at 20% prevalence | No; batch structure assessed using `PLSDAbatch`, with no correction applied | Feature-wise z-scoring |
| `Mucus Microbiota` | 163/168 | 47 | CLR transformation | ASV-level count matrix filtered at 20% prevalence | No; batch structure assessed using `PLSDAbatch`, with no correction applied | Feature-wise z-scoring |

Host physiological traits and global liver cytosine-composition variables are **not MOFA views**. They are retained as sample metadata and used to help interpret factor scores.

## Directory layout

- `input_files/` — processed omics matrices, experiment metadata, and host-trait source tables required to reproduce the integration.
- `metadata_for_omics_data.rds` — assembled design metadata, host traits, and liver cytosine-composition variables.
- `views.rds` — final feature-wise z-scored MOFA input matrices (`features × samples`).
- `views_not_zscored_for_MultiPower.rds` — aligned, filtered matrices before feature-wise z-scoring; used only by `07_MultiPower.R`.
- `covariates.rds` — aligned sample metadata used by MOFA and MultiGroupPower.
- `QC_views/` — independent input quality-control tables and plots.
- `mofa_models/` — trained MOFA HDF5 model files and session information.
- `elbo_mofa_models/` — ELBO and factor-reproducibility comparisons across trained models.
- `mofa_results/` — MOFA interpretation, overview, and enrichment outputs.
- `MultiGroupPower/` — equal-size and unequal-size retrospective multi-omics power outputs.

Run scripts from `09_Multi_omics/` so their relative paths resolve correctly.

## Required inputs

### Experiment metadata and host-trait source tables

`01_metadata_for_omics_preprocess_dataset.R` reads:

```text
input_files/STPN2309_metadata_all.tsv
input_files/plasma_glucose.tsv
input_files/plasma_dlactate.tsv
input_files/plasma_llactate.tsv
input_files/body_weight.tsv
input_files/hepatosomatic_index.tsv
input_files/methylome_for_dirichlet.tsv
```

The metadata export retains:

```text
day_categ
diet
group
phase
diet_phase
```

It also assembles the following sample-level covariates:

```text
plasma_glucose
HSI
body_weight
plasma_DLactate
plasma_LLactate
dC_rel
mdC_rel
hmdC_rel
```

For plasma glucose, replicate columns beginning with `replicate_` are averaged within sample.

### Processed omics matrices

`02_data_preparation_for_mofa.R` and `03_QC_views.R` expect:

```text
input_files/H3_H4_5PTM_midgut_batch_corrected_limma.rds
input_files/H3_H4_8PTM_liver_batch_corrected_limma.rds
input_files/midgut_proteome_filtered_batch_corrected_limma.rds
input_files/liver_proteome_filtered_batch_corrected_limma.rds
input_files/digesta_ASVbest_sample_by_feature.rds
input_files/mucus_ASVbest_sample_by_feature.rds
```

The hPTM and proteome RDS files contain named list objects. The scripts extract:

```text
ptm_batch_corrected_limma
proteome_batch_corrected_limma
```

respectively.

### STRING annotation resources

MOFA protein-weight interpretation and functional enrichment require two external STRING v12 rainbow trout files:

```text
STRG0A55HWH.protein.info.v12.0.txt
110079946.protein.enrichment.terms.v12.0.txt
```

They are not committed because they are large external annotation resources. Download them from the [STRING v12 rainbow trout organism page](https://version-12-0.string-db.org/organism/STRG0A55HWH) and place them in `09_Multi_omics/`.

## Script order

1. `01_metadata_for_omics_preprocess_dataset.R` — assembles design metadata, host traits, and liver cytosine-composition variables.
2. `02_data_preparation_for_mofa.R` — creates aligned MOFA views, the associated covariate table, and non-z-scored matrices for MultiGroupPower.
3. `03_QC_views.R` — independently evaluates original and feature-wise z-scored input distributions, missingness, feature dispersion, and variance concentration.
4. `04_training_models.R` — training template for one MOFA model per chosen factor number and random seed.
5. `05_compare_models.R` — compares trained models by ELBO and factor reproducibility.
6. `06_a_MOFA_downstream_selected_factors.R` — detailed factor interpretation, factor weights, factor–covariate associations, and GO Biological Process enrichment.
7. `06_b_mofa_overview.R` — generates the combined MOFA overview figure and reusable source tables.
8. `07_MultiPower.R` — retrospective MultiGroupPower analysis for the seven day-specific HC-versus-NC comparisons.

Scripts 02 and 03 both require the processed input matrices, but script 03 is an independent QC workflow: it does not overwrite `views.rds` or alter the model inputs produced by script 02.

## 1. Assembly of metadata and non-view covariates

`01_metadata_for_omics_preprocess_dataset.R` creates:

```text
metadata_for_omics_data.rds
```

This RDS file contains three elements:

```text
metadata
host_traits
Liver_DNA_methylation
```

`metadata` contains the categorical experimental design variables. `host_traits` contains glucose, HSI, body weight, and L-/D-lactate variables. `Liver_DNA_methylation` contains the relative global DNA fractions `dC_rel`, `mdC_rel`, and `hmdC_rel`.

These host and DNA variables are attached to the MOFA sample metadata during data preparation. They are later used for factor–covariate correlation plots, but they are not themselves included as latent-factor views.

## 2. Preparation of the six MOFA views

`02_data_preparation_for_mofa.R` performs the final data assembly for model training.

### Sample alignment

Each source matrix is converted from `samples × features` to MOFA's required orientation:

```text
features × samples
```

The script aligns every view to the same metadata sample universe.

A sample that is absent from a given molecular assay is represented by an all-`NA` column in that view. This allows MOFA2 to use all available measurements without requiring a complete six-assay profile for every fish.

### Feature handling before MOFA

The integration-specific filtering is intentionally limited:

- **Midgut and liver hPTMs:** retained without an additional MOFA-specific feature filter.
- **Midgut and liver proteomes:** retained only when the feature standard deviation across available samples is at least `0.30`.
- **Digesta and mucus microbiota:** prevalence is calculated for diagnostic purposes only; no additional feature removal is performed in this script. The input microbiota table is the same prevalence-filtered table used for the differential-abundance and correlation analyses.

All six final MOFA views are feature-wise z-scored across samples. Features with zero or non-finite standard deviation are excluded by the z-scoring helper.

Feature names are prefixed with their source view, for example:

```text
Midgut_
Liver_
Digesta_
Mucus_
```

This prevents collisions between identical IDs originating from different assays or tissues.

### Outputs

The script writes:

```text
views.rds
covariates.rds
views_not_zscored_for_MultiPower.rds
```

`views.rds` is the training input. It is a named list of the six z-scored matrices, each arranged as `features × samples`.

`covariates.rds` contains the aligned sample metadata, host traits, and liver cytosine-composition variables.

`views_not_zscored_for_MultiPower.rds` contains the same filtered and aligned six matrices before feature-wise z-scoring. It is reserved for `07_MultiPower.R`.

## 3. Input QC

`03_QC_views.R` is a standalone quality-control script. It reloads the processed molecular matrices and independently creates original and feature-wise z-scored representations for QC only.

For each view, it reports:

- sample and feature dimensions;
- feature-level missingness and zero fractions;
- distribution of original and z-scored values;
- feature mean, median, variance, standard deviation, range, and skewness;
- the number of features retained after QC z-scoring;
- feature-variance ranks and cumulative variance distributions.

It exports:

```text
QC_views/tables/omics_view_dimensions_original.csv
QC_views/tables/zscore_feature_filtering_report.csv
QC_views/tables/feature_level_QC_original_vs_zscore.csv
QC_views/tables/variance_rank_original_vs_zscore.csv
```

PDF figures compare original versus z-scored distributions, feature standard deviations, and cumulative variance concentration.

This QC script does not determine the final feature filtering applied in `02_data_preparation_for_mofa.R`; it documents the scale, missingness, and distributional consequences of the input processing.

## 4. MOFA2 model training

`04_training_models.R` is a model-training template.

Before running it, update the local Python executable used through `reticulate`:

```r
py <- "C:/mofa_venv/Scripts/python.exe"
```

The environment must contain the Python backend required by MOFA2.

The configured MOFA2 settings are:

```r
data_opts$scale_views <- TRUE

model_opts$num_factors <- 20
model_opts$spikeslab_factors <- FALSE
model_opts$spikeslab_weights <- TRUE

train_opts$convergence_mode <- "slow"
train_opts$maxiter <- 5000
train_opts$weight_views <- FALSE
```

The random seed and HDF5 output name are intentionally placeholders in the template and must be set for each fit.

The final downstream scripts currently reload:

```text
mofa_models/K20_setseed_42.hdf5
```

Thus, all interpretation figures and enrichment analyses must use the same selected model file unless the input model path is changed deliberately.

Although `views.rds` has already been feature-wise z-scored, `scale_views = TRUE` is also retained in the MOFA2 configuration so that the fitted model applies its configured view scaling.

## 5. Model comparison

`05_compare_models.R` loads all `.hdf5` files present in `mofa_models/`.

It generates:

- an evidence lower bound (ELBO) comparison table and plot;
- a clustered factor-comparison plot across models.

Outputs are written to:

```text
elbo_mofa_models/ELBO_comparison.pdf
elbo_mofa_models/ELBO_comparison.tiff
elbo_mofa_models/ELBO_comparison_data.csv
elbo_mofa_models/factor_comparison.pdf
elbo_mofa_models/factor_comparison.tiff
```

ELBO comparison is meaningful only among models fitted to the same input views and preprocessing. The factor-comparison heatmap helps assess whether corresponding latent patterns are reproducible across model choices or random seeds.

## 6. Factor interpretation and biological annotation

`06_a_MOFA_downstream_selected_factors.R` reloads the selected HDF5 model, standardises display labels for metadata, and writes all results under:

```text
mofa_results/
```

### Model structure and variance explained

The script exports:

- data availability overview;
- factor-correlation heatmap;
- per-factor, per-view variance explained;
- total variance explained per view;
- cumulative variance-explained curves;
- contribution-score and per-sample variance-explained objects.

These outputs identify which molecular views contribute most to each factor and how much of each view is represented by the model.

### Factor–covariate associations and visualisation

The script correlates all factors with:

```text
Plasma Glucose
Plasma L-Lactate
Plasma D-Lactate
HSI
Unmodified C
5-mC
5-hmC
```

and displays the associations as a log-p-value heatmap.

It also creates targeted visualisations for:

- Factor 3 scores by day and diet;
- Factor 2 versus Factor 3, coloured by diet and shaped by temporal phase;
- UMAP calculated from all model factors, coloured by diet and shaped by temporal phase;
- factor-specific heatmaps for Factor 3 across the six views;
- top positive and negative feature weights for Factor 3 in each view.

Protein labels used in interpretation plots are added only after model fitting. The original MOFA feature names are not changed.

The sign of a MOFA factor is arbitrary: an equivalent model may reverse both factor scores and feature weights. Interpret a factor from its score pattern and positive/negative weights together, not from the sign alone.


## 7. MOFA overview figure

`06_b_mofa_overview.R` produces the combined overview figure from the selected model:

- **Panel A:** sample availability across the six views;
- **Panel B:** total variance explained per view;
- **Panel C:** factor-wise variance explained across views.

It writes PDF and 1200 dpi TIFF versions of the tagged and untagged figure, together with reusable source tables:

```text
mofa_results/tables/MOFA_view_N_and_D_summary.csv
mofa_results/tables/MOFA_total_R2_per_view.csv
mofa_results/tables/MOFA_factor_R2_per_view.csv
```

These tables are the preferred source for reporting the final number of available samples (`N`) and retained features (`D`) per view, because they are calculated directly from the selected model.

## MOFA2 documentation and references

Model training, variance-explained assessment, factor interpretation, feature-weight extraction, and functional enrichment analyses were implemented following the [MOFA2 tutorials and documentation](https://biofam.github.io/MOFA2/tutorials.html).

## Reference

Argelaguet, R., Arnol, D., Bredikhin, D., et al. (2020). MOFA+: a statistical framework for comprehensive integration of multi-modal single-cell data. *Genome Biology*, 21, 111. https://doi.org/10.1186/s13059-020-02015-1


## 9. Retrospective multi-omics power analysis

`07_MultiPower.R` uses the following input objects:

```text
views_not_zscored_for_MultiPower.rds
covariates.rds
```

The analysis is performed directly on the six assay-specific data matrices and does not use MOFA factors or factor scores.

For each omics view, samples are matched to the metadata using their sample identifiers. Only samples assigned to the NC or HC diets and to one of the seven sampling days are retained. The script does not apply any additional feature filtering, transformation, imputation, or manual estimation of pilot parameters. The matrices are supplied to `MultiGroupPower` as normalized continuous data:

```r
type = 2
```

The following seven within-day comparisons are evaluated:

```text
NC_day_01 versus HC_day_01
NC_day_02 versus HC_day_02
NC_day_03 versus HC_day_03
NC_day_04 versus HC_day_04
NC_day_10 versus HC_day_10
NC_day_15 versus HC_day_15
NC_day_22 versus HC_day_22
```

The analysis includes the following omics views:

```text
Midgut hPTMs
Liver hPTMs
Midgut Proteome
Liver Proteome
Digesta Microbiota
Mucus Microbiota
```

`MultiGroupPower()` is run using an equal-size design:

```r
equalSize = TRUE
```

This configuration assumes balanced replication, with the same number of samples in the HC and NC groups for each comparison and across the included omics views.

The specified design targets are:

```text
minimum power per omics view: 0.60
average power across views:   0.80
feature-level FDR:            0.05
maximum search size:          2000 samples per group
relative omics cost:          1 for every view
```

The maximum sample size is a computational search limit and should not be interpreted as a recommended experimental sample size.

Following the `MultiGroupPower()` analysis, `postMultiPower()` is applied separately to each day-specific comparison. This retrospective step constrains the evaluation to a maximum of 12 biological replicates per diet and sampling day, corresponding to the replication of the original experimental design:

```r
max.size = 12
```

The resulting tables report, for each comparison and omics view, the estimated power, sample-size information, and number of detectable features under this replication constraint.

Outputs are written to:

```text
MultiGroupPower/
├── figures/
├── tables/
├── rds/
└── sessionInfo.txt
```

The exported files include:

* package-generated `MultiGroupPower` and `postMultiPower` plots;
* the matrix of requested day-specific comparisons;
* the `MultiGroupPower` global summary;
* day-specific `postMultiPower` power, sample-size, and detectable-feature tables;
* serialized input and result objects;
* R session information.

These analyses provide retrospective, pilot-data-based power and sample-size evaluations for balanced, within-day HC-versus-NC comparisons. They do not test the overall diet-by-day interaction and do not replace the statistical differential analyses performed separately for each omics dataset.

### MultiPower documentation and reference

The power and sample-size analyses were performed using the `MultiPower` R package. Script structure, input preparation, and function use were adapted from the example workflows and User’s Guide provided by the developers in the [MultiPower repository](https://github.com/ConesaLab/MultiPower).

### Reference

Tarazona, S., Balzano-Nogueira, L., Gómez-Cabrero, D., et al. (2020). Harmonization of quality metrics and power calculation in multi-omic studies. *Nature Communications*, 11, 3092. https://doi.org/10.1038/s41467-020-16937-8

## Main outputs

- `metadata_for_omics_data.rds` — integrated metadata, host traits, and liver cytosine-composition variables.
- `views.rds` — final six z-scored MOFA views.
- `covariates.rds` — aligned covariates and design metadata.
- `views_not_zscored_for_MultiPower.rds` — non-z-scored six-view input for MultiGroupPower.
- `QC_views/` — input-QC tables and figures.
- `mofa_models/K20_setseed_42.hdf5` — selected model used by the downstream scripts.
- `elbo_mofa_models/` — model ELBO and reproducibility comparisons.
- `mofa_results/` — factor interpretation, variance decomposition, weights, enrichment, and overview figure outputs.
- `MultiGroupPower/` — equal-size and unequal-size prospective sample-size outputs.

## Main R packages

Core packages are `MOFA2`, `reticulate`, `data.table`, `dplyr`, `tidyr`, `tibble`, `matrixStats`, `ggplot2`, `patchwork`, `ggpubr`, `ComplexHeatmap`, `circlize`, `readr`, `stringr`, `forcats`, `purrr`, `cowplot`, `scales`, `grid`, `MultiPower`, `FDRsampsize`, and `lpSolve`.
