# Host physiological traits

## Purpose

This module analyses body weight, plasma glucose, hepatosomatic index (HSI), plasma L-lactate, and plasma D-lactate in rainbow trout across dietary treatment and sampling day. It also provides pooled cross-trait Spearman correlations, descriptive day-specific coefficient-of-variation (CV) comparisons, and observed sample-size tables.

## Directory layout

- `input_files/` — tab-separated input tables for each physiological trait.
- `results/` — trait-specific modelling outputs and secondary analyses.
- `RDS/` — saved final `ggplot` objects, including the processed data used for the published trait panels.

Run all scripts from `01_Host_Physiological_Traits/`.

## Inputs

Place the following tab-separated files in `input_files/`:

- `body_weight.tsv`
- `plasma_glucose.tsv`
- `hepatosomatic_index.tsv`
- `plasma_llactate.tsv`
- `plasma_dlactate.tsv`

All trait tables require, at minimum, a diet variable (`diet`) and sampling-day variable (`day`).

Expected trait columns are:

| File | Modelled variable | Additional preprocessing |
|---|---|---|
| `body_weight.tsv` | `body_weight` | None |
| `plasma_glucose.tsv` | `plasma_glucose` | Calculated as the row mean of `replicate_1`, `replicate_2`, and `replicate_3` |
| `hepatosomatic_index.tsv` | `HSI` | None |
| `plasma_llactate.tsv` | `Plasma_L_Lactate` | None |
| `plasma_dlactate.tsv` | `Plasma_D_Lactate` | Values are multiplied by 1,000 before trait-specific modelling and plotting, yielding µM |

Sampling days are ordered as: day 1, 2, 3, 4, 10, 15, and 22. `NC` is the reference diet level in all trait models.

## Script order

Run the five trait-specific scripts first, as they generate the saved plot objects required by the CV and sample-size scripts.

1. `01_body_weight.R` — body-weight modelling and reporting.
2. `02_plasma_glucose.R` — glucose technical-replicate processing, modelling, and reporting.
3. `03_hepatosomatic_index.R` — HSI modelling and reporting.
4. `04_plasma_llactate.R` — plasma L-lactate modelling and reporting.
5. `05_plasma_dlactate.R` — plasma D-lactate unit conversion, modelling, and reporting.
6. `06_correlation_analysis.R` — pooled complete-case correlations among glucose, HSI, L-lactate, and D-lactate.
7. `07_CV_analysis.R` — descriptive day-specific CV analysis using data embedded in the saved trait plots.
8. `08_sample_size.R` — observed sample-size extraction from the same saved plot objects.

## Trait-specific modelling

### Fixed-effects design

All five trait-specific scripts fit the same mean model:

```r
trait ~ day * diet
```

### Candidate distributions and final model

Each script fits the following candidate models with `glmmTMB`:

- Gaussian model on the original response scale;
- Gaussian model on the log-transformed response scale;
- Student-t model on the original response scale;
- Gamma model with a log link.

For each distribution, a homoscedastic model and a model with condition-specific dispersion are fitted. The latter includes:

```r
dispformula = ~ day * diet
```

Candidate-fit summaries report family, link, mean and dispersion formulas, AIC, BIC, log-likelihood, convergence status, and Hessian status. Likelihood-ratio tests compare homoscedastic and variable-dispersion fits within the same family and response scale.

The final model used for all five traits is:

```r
trait ~ day * diet
family = Gamma(link = "log")
dispformula = ~ day * diet
```

This specification models positive, right-skewed trait values on a multiplicative scale while allowing residual dispersion to differ among diet-by-day conditions. The log-Gaussian and Gamma alternatives require strictly positive response values; scripts with explicit validation stop if non-positive values are detected.

### Diagnostics

Residual diagnostics are generated with `DHARMa`, using 1,000 simulated residual sets per model. The diagnostics include:

- residual uniformity;
- dispersion;
- outliers;
- graphical residual checks.

The final-model outputs also report `performance::check_convergence()`, `insight::is_converged()`, optimiser convergence information, and the positive-definite Hessian check.

### Diet contrasts and figures

Estimated marginal means are calculated with `emmeans` for diet within each day:

```r
emmeans(final_model, ~ diet | day)
```

HC-versus-NC contrasts are obtained separately for the seven sampling days. Raw contrast p-values are adjusted across those day-specific contrasts using the Benjamini–Hochberg procedure. Because the final model has a log link, the contrast output includes:

```r
fold_change = exp(estimate)
```

which represents the estimated HC/NC ratio.

Final trait figures show individual fish as beeswarm points, faceted by sampling day. Only day-specific contrasts with BH-adjusted p-values below 0.05 are annotated.

## Cross-trait correlation analysis

`06_correlation_analysis.R` harmonises the plasma L-Lactate, plasma D-Lactate, Hepato-somatic Index (HSI), and plasma glucose datasets by sample identifier, diet, day, and rearing tank.

The script retains complete cases and calculates the pairwise Spearman rank correlations among:

- plasma glucose;
- HSI;
- plasma L-lactate;
- plasma D-lactate.

Pairwise p-values are adjusted using the Benjamini–Hochberg procedure across the tests.

The script additionally generates a plasma-glucose-versus-HSI scatterplot, coloured by diet. Its fitted linear trend is for visualisation only; the displayed correlation annotation reports the unadjusted Spearman result.

## Coefficient-of-variation analysis

`07_CV_analysis.R` loads the five saved trait plot objects from `RDS/` and extracts their processed data. For each trait, diet, and sampling day, it calculates:

```text
CV (%) = 100 × SD / mean
```

This produces seven day-specific CV values per diet for each trait. These CV values are displayed as faceted NC-versus-HC boxplots with individual day-level points. Wilcoxon tests compare the distributions of day-specific CV values between diets.

## Observed sample-size extraction

`08_sample_size.R` scans the `.rds` plot objects in `RDS/`, detects each plot's response variable, and counts non-missing observations by trait, diet, and sampling day.

It writes:

- a long-format sample-size table;
- a day-by-diet wide summary used for reporting;
- a plain-text publication sentence with sample sizes shown as HC/NC for each time point.

## Outputs

### Trait-specific outputs

Each directory under `results/<trait>/` contains:

- final faceted beeswarm plot in PDF and TIFF formats;
- final DHARMa residual diagnostic plot in PDF and TIFF formats;
- candidate-model comparison table (`*_model_comparison.csv`);
- HC-versus-NC contrast table (`*_contrasts.csv`);
- final-model summary and diagnostic test output (`*_model_summary.txt`);
- `sessionInfo.txt`.

The associated `RDS/<trait>_plot.rds` file stores the final figure object and its processed underlying data.

### Secondary outputs

- `results/correlation_analysis/` — complete-case input table, pairwise Spearman results, rho and p-value matrices, glucose-versus-HSI statistics, and the corresponding figure.
- `results/CV_analysis/` — CV comparison figure in PDF and TIFF formats.
- `results/sample_size/` — observed sample-size CSV and publication-ready text.

## Main R packages

Core packages are `glmmTMB`, `emmeans`, `DHARMa`, `performance`, `insight`, `ggplot2`, `ggbeeswarm`, `ggsignif`, `dplyr`, `tidyr`, `purrr`, `readr`, `ggpubr`, and `rlang`.

