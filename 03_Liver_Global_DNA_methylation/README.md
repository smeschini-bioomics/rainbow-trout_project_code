# Global liver DNA cytosine composition

## Purpose

This module analyses global liver DNA cytosine composition measured by HPLC–UV. The three measured cytosine forms—unmodified cytosine (`dC`), 5-methylcytosine (`mdC`), and 5-hydroxymethylcytosine (`hmdC`) are modelled jointly as a composition rather than as three independent outcomes.

The workflow produces:

- relative-composition data for modelling;
- Bayesian Dirichlet regression results for `diet * day`;
- day-specific HC-minus-NC posterior contrasts;
- a contrast heatmap;
- raw-data beeswarm figures and descriptive summaries;
- day-specific coefficient-of-variation summaries.


## Script order

1. `01_preparing_data.R` — calculates relative cytosine composition and writes `methylome_for_dirichlet.tsv`.
2. `02_dirichlet_regression.R` — validates the composition, fits the Bayesian Dirichlet regression, calculates posterior predictions and HC-minus-NC contrasts, and saves model outputs.
3. `03_dirichlet_regression_heatmap_results.R` — converts posterior contrast results into the manuscript heatmap.
4. `04_data_beeswarm_plots.R` — generates raw-data beeswarm plots in percentage units and descriptive sample-level summaries.
5. `05_CV_analysis.R` — calculates day-specific CV values and compares their distributions between diets.

## Bayesian Dirichlet regression

### Model structure

`02_dirichlet_regression.R` fits a multivariate Dirichlet regression with `brms`:

```r
cbind(dC_rel, mdC_rel, hmdC_rel) ~ diet * day
```

The Dirichlet likelihood models the three fractions jointly and respects their unit-sum constraint. NC is set as the diet reference level. The `diet * day` interaction allows the HC-versus-NC composition difference to vary across sampling days.

### Input checks

Before model fitting, the script:

1. checks that `diet`, `day`, `dC_rel`, `mdC_rel`, and `hmdC_rel` are available;
2. removes rows with missing values in the required fields;
3. confirms that all composition values are finite, strictly positive, below one, and sum to one.

The script stops if these requirements are not met.

### Prior specification

The model uses weakly informative priors with empirical centring for baseline composition.

For the two non-reference Dirichlet intercepts, prior centres are calculated from NC at day 1 using log ratios of the mean component fractions relative to `dC`:

```r
mu_md_vs_dC  = log(mean(mdC_rel)  / mean(dC_rel))
mu_hmd_vs_dC = log(mean(hmdC_rel) / mean(dC_rel))
```

The corresponding intercept priors are:

```r
student_t(3, mu_md_vs_dC, 1)
student_t(3, mu_hmd_vs_dC, 1)
```

Diet, day, and diet-by-day regression coefficients use:

```r
normal(0, 0.5)
```

The Dirichlet precision parameter uses:

```r
exponential(0.25)
```

### MCMC settings and saved model

The fitted model uses:

```r
chains = 4
iter = 3000
warmup = 1000
seed = 123
init = 0
adapt_delta = 0.9995
max_treedepth = 15
```

All model parameters are retained for downstream diagnostics and model comparison. The fitted object is saved as:

```text
results/dirichlet_regression/fit_dirichlet_diet_day.rds
```

## Posterior contrasts and interpretation

### HC-minus-NC contrasts

For each day and each cytosine component, the script obtains posterior expected compositions for HC and NC with `posterior_epred()`. The day-specific contrast is:

```text
HC expected fraction − NC expected fraction
```

For each posterior contrast distribution, the output includes:

- posterior mean (`estimate`);
- central 95% credible interval (`lower95`, `upper95`);
- posterior probability that the contrast is positive (`prob_pos`);
- direction (`increase` or `decrease`);
- whether the 95% credible interval excludes zero (`sig95`).

The script also converts effects to percentage points:

```r
estimate_pp = 100 * estimate
```

and calculates posterior directional certainty:

```r
prob_certainty = max(P(contrast > 0), P(contrast < 0))
```

This is equivalent to a probability-of-direction measure. A value above 0.95 indicates strong posterior support for the estimated direction.

The complete contrast table is saved as:

```text
results/dirichlet_regression/contrasts_HC_minus_NC_by_day_response.csv
```

### Posterior prediction summaries

The script computes posterior expected compositions for every observed diet-by-day combination. It saves a figure combining individual fractions with posterior mean estimates and 95% credible intervals:

```text
results/dirichlet_regression/jitter_with_posterior_CI_dirichlet.pdf
```

The script also compares posterior fitted means with observed fractions in:

```text
results/dirichlet_regression/fitted_vs_observed_dirichlet.pdf
```

These are model-fit visualisations based on posterior expected values; they are not separate frequentist tests.

### Model summary and LOO output

The workflow writes:

- `diagnostics_summary.txt` — `brms` model summary;
- `brms_fixed_effects_dirichlet.csv` — fixed-effect posterior summaries;
- `loo_dirichlet.rds` and `loo_dirichlet.txt` — leave-one-out cross-validation results;
- `sessionInfo.txt`.

LOO is retained for Bayesian model diagnostic and comparison purposes. This workflow fits one prespecified `diet * day` model; it does not use LOO to select among a set of alternative final models.

## Contrast heatmap

`03_dirichlet_regression_heatmap_results.R` reads the HC-minus-NC posterior contrast table and creates a heatmap with:

- rows: `C`, `5-mC`, and `5-hmC`;
- columns: sampling days;
- point colour: direction of the HC-minus-NC effect;
- point size: absolute effect magnitude in percentage points;
- point outline: whether the 95% credible interval excludes zero;
- text label: posterior directional certainty, displayed only when above 0.95.

The heatmap is exported to:

```text
results/dirichlet_regression/heatmap_centered_probabilities.pdf
results/dirichlet_regression/heatmap_centered_probabilities.tiff
```

## Raw-data beeswarm plots and descriptive summaries

`04_data_beeswarm_plots.R` reads the derived composition table, removes `Fasted` samples, and converts fractions to percentages:

```r
dC_pct   = 100 * dC_rel
mdC_pct  = 100 * mdC_rel
hmdC_pct = 100 * hmdC_rel
```

It produces separate day-faceted beeswarm plots for:

- unmodified cytosine (`C`);
- 5-methylcytosine (`5-mC`);
- 5-hydroxymethylcytosine (`5-hmC`).

A vertically combined three-panel figure with a shared diet legend is also generated. These plots are descriptive and do not include statistical annotations.

The script exports:

```text
results/beeswarm_plot/raw_beeswarm_dC.pdf
results/beeswarm_plot/raw_beeswarm_mdC.pdf
results/beeswarm_plot/raw_beeswarm_hmdC.pdf
results/beeswarm_plot/raw_beeswarm_methylome_combined.pdf
```

It also writes:

- sample counts by component, diet, and day;
- a wide NC/HC sample-count table;
- mean, SD, median, first and third quartiles, minimum, and maximum by component, diet, and day;
- global descriptive summaries by cytosine form;
- session information.

## Coefficient-of-variation analysis

`05_CV_analysis.R` calculates the coefficient of variation for each component, diet, and day:

```text
CV (%) = 100 × SD / mean
```

The resulting dataset contains one CV value per `component × diet × day` group. The script compares the distributions of the seven day-specific CV values between NC and HC separately for `C`, `5-mC`, and `5-hmC` using Wilcoxon rank-sum tests.

Only unadjusted p-values at or below 0.05 are annotated in the CV figure. This analysis is descriptive: it compares day-level variability summaries, not individual-level CV values, and the current script does not apply multiple-testing correction.

Outputs include:

```text
results/CV_analysis/cv_percent_relative_all_data_by_day_diet.csv
results/CV_analysis/CV_boxjitter_relative_all_data.pdf
results/CV_analysis/CV_boxjitter_relative_all_data.tiff
results/CV_analysis/n_CV_by_diet_methylome_relative_all_data.csv
```

## Sample-size and descriptive-composition outputs

`02_dirichlet_regression.R` additionally writes a diet-by-day composition summary before fitting the model:

```text
results/sample_size/methylome_summary_by_diet_day.csv
```

This table contains `n`, mean, SD, and CV for each relative cytosine component by diet and day. It reports observed sample sizes and descriptive statistics only; it is not a prospective power or sample-size calculation.

## Main outputs

- `input_files/methylome_for_dirichlet.tsv` — relative composition table derived from the raw HPLC–UV measurements.
- `results/dirichlet_regression/` — fitted model, posterior contrasts, posterior expected-value plots, heatmap, fixed effects, LOO output, and session information.
- `results/beeswarm_plot/` — raw percentage-scale beeswarm figures, descriptive statistics, and sample-count tables.
- `results/CV_analysis/` — day-specific CV table and CV comparison figure.
- `results/sample_size/` — observed composition summary by diet and day.

## Main R packages

Core packages are `brms`, `posterior`, `loo`, `broom.mixed`, `readr`, `dplyr`, `tidyr`, `ggplot2`, `ggbeeswarm`, `ggnewscale`, `ggpubr`, `cowplot`, and `scales`.
