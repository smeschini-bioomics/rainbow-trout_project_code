# Targeted liver mRNA analysis

## Purpose

This module analyses targeted liver mRNA relative-concentration measurements across dietary treatment and sampling day. It produces gene-specific descriptive and inferential outputs, a supplementary multi-panel figure for selected glycolytic and gluconeogenic targets, observed sample-count tables, and a concise per-gene model-summary table.


## Script order

1. `01_mRNA_level_analysis.R` — reusable analysis template. Set the input path and selected final model, then run once for each target gene.
2. `02_sample_size.R` — derives observed sample counts from the retained gene-level plot objects.
3. `03_summary_all_genes.R` — creates a minimal all-gene table recording whether an inferential analysis was completed and the selected model family/link.

Run the gene-level analyses before the supplementary-figure and summary scripts.

## Single-gene workflow

### 1. Data preparation and eligibility for inferential analysis

For each target gene, `01_mRNA_level_analysis.R`:

1. imports the selected input table;
2. applies the defined day and diet factor ordering;
3. checks the number of non-missing observations in every planned `Diet × Day` group.
4. generates interactive exploratory histograms and density plots of `RC`;

Inferential analysis requires at least three non-missing `RC` observations in every NC-by-day and HC-by-day group. A minimum sample size alone is not sufficient: the response values must also be compatible with at least one valid candidate model, and the selected model must converge with acceptable diagnostics.

When this eligibility requirement is not met, or when no defensible model can be obtained, the inferential branch is skipped. No fitted model, contrasts, p-values, or significance annotations are generated. The descriptive beeswarm plot without statistical annotations is still retained.

Models requiring strictly positive response values cannot be fitted when `RC` contains zero or negative values:

- Gaussian on `log(RC)`;
- log-normal;
- Gamma with log link.

No arbitrary pseudocount is added automatically.

### 2. Candidate models and model selection

All candidate models use the same fixed-effects structure:

```r
RC ~ day * diet
```

The script evaluates five response-model specifications:

- Gaussian model with identity link;
- Gaussian model fitted to `log(RC)`;
- log-normal model with log link;
- Student-t model with identity link;
- Gamma model with log link.

For each candidate, DHARMa simulated residuals are generated using 1,000 simulations. Model choice is made separately for each gene, with primary emphasis on:

- numerical convergence and a positive-definite Hessian;
- DHARMa residual behaviour, including uniformity, dispersion, outliers, and graphical diagnostics.

AIC, BIC, and log-likelihood are used as supporting fit criteria, not as the sole selection rule. Where more than one family provides similarly adequate diagnostics and fit, a common model family is preferred across genes to maintain a parsimonious and consistent analytical framework.

The final family is set manually in the template for the gene being analysed, for example:

```r
final_model_id <- "student"
```

Valid model identifiers are:

```text
gaussian
gaussian_log
lognormal
student
gamma
```

The selected family is gene specific. It must not be assumed that all targets use the same response distribution.

### 3. Diet contrasts and effect-scale interpretation

Estimated marginal means are calculated for diet within each sampling day:

```r
emmeans(final_model, ~ diet | day, type = "link")
```

The script obtains reversed pairwise contrasts so that the estimate is reported in the direction `HC - NC`. Raw p-values are adjusted across the seven day-specific contrasts for each gene using the Benjamini–Hochberg procedure.

The scientific interpretation of the contrast depends on the final model scale:

| Final model | Meaning of the HC−NC estimate | `exp(estimate)` |
|---|---|---|
| Gaussian, identity link | HC−NC difference in relative-concentration units | Not reported |
| Student-t, identity link | HC−NC difference in relative-concentration units | Not reported |
| Gaussian on `log(RC)` | HC−NC contrast on the log scale | HC/NC ratio |
| Log-normal, log link | HC−NC contrast on the log scale | HC/NC ratio |
| Gamma, log link | HC−NC contrast on the log scale | HC/NC ratio |

For models fitted with a logarithmic response scale or log link, exponentiated contrasts are reported as HC/NC ratios. For Gaussian and Student-t models with an identity link, contrasts are reported as additive HC−NC differences in relative-concentration units.

### 4. Gene-level figures

For every analysed gene, the primary plot displays individual `RC` values as beeswarm points and is faceted by day.

Only HC-versus-NC comparisons with BH-adjusted p-values below 0.05 are annotated on the statistical plot.

A second descriptive beeswarm plot is also generated without brackets or p-value annotations. This no-stat plot is retained for genes where inferential analysis is intentionally omitted or cannot be justified.

## Model diagnostics and reporting

For each gene with a completed inferential analysis, the workflow saves:

- the selected final-model summary;
- a DHARMa residual-diagnostic figure based on 1,000 simulations;
- the day-specific HC−NC contrast table with raw and BH-adjusted p-values;
- effect-scale metadata describing whether the estimate is an additive difference or ratio-scale effect;
- session information.

For genes where statistical analysis is skipped, the workflow retains the descriptive no-stat figure and exports the reason for skipping together with the `Diet × Day` group-size table.


## Observed sample-size and all-gene model summary

`03_sample_size.R` reads the intended annotated plot objects in `RDS/` and counts non-missing observations by gene, day, and diet. It reports observed sample counts only; it does not perform power analysis or estimate future sample-size requirements.

`04_summary_all_genes.R` includes every input TSV in a concise summary table. For each gene, it records:

- whether an inferential analysis was completed;
- the selected model family;
- the selected link.

When no model summary is present because an analysis was skipped, the statistical-analysis and model fields are reported as `NA`. This applies, for example, to genes such as `gcka` and `gckb` when a defensible inferential model was not available.

## Outputs

### Gene-specific outputs

For a gene named `<gene>`, completed inferential analysis writes to `results/<gene>/`:

- `<gene>_facet_plot.pdf` and `.tiff` — day-faceted beeswarm plot with significant HC−NC annotations;
- `<gene>_facet_plot_no_stat.pdf` and `.tiff` — descriptive beeswarm plot without statistical annotations;
- `<gene>_DHARMa_residuals.pdf` and `.tiff` — selected-model residual diagnostics;
- `<gene>_contrasts.csv` — day-specific HC−NC contrasts;
- `<gene>_model_summary.txt` — selected-model summary and model-scale metadata;
- `sessionInfo.txt`.

The script also saves:

- `RDS/<gene>_plot.rds` — annotated plot object;
- `RDS/<gene>_plot_no_stat.rds` — unannotated plot object.

When inference is skipped, only outputs appropriate to the descriptive branch are produced.

### Supplementary and summary outputs

- `results/sample_size/` — observed sample-size table.
- `results/summary_all_genes/all_genes_model_summary.csv` — per-gene inferential-analysis status and selected family/link.

## Main R packages

Core packages are `glmmTMB`, `emmeans`, `DHARMa`, `ggplot2`, `ggbeeswarm`, `ggsignif`, `patchwork`, `dplyr`, `tidyr`, `purrr`, `readr`, `stringr`, and `rlang`.
