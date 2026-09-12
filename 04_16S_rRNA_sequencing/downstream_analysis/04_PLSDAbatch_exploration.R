# libraries
suppressPackageStartupMessages({
  library(microeco)
  library(mixOmics)
  library(PLSDAbatch)
  library(vegan)
  library(ggplot2)
  library(gridExtra)
  library(grid)
})

# select the dataset
input_file <- file.path("Rdata", "16s_dataset.RData")
# set the output directory
output_dir <- "04_PLSDAbatch_exploration"
# check
stopifnot(file.exists(input_file))

dir.create(file.path(output_dir, "Digesta"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "Mucus"), recursive = TRUE, showWarnings = FALSE)


# Load data and retain Digesta/Mucus samples from NC and HC fish
load(input_file)
# check
if (!exists("mt_raw")) {
  stop("Object 'mt_raw' was not found in: ", input_file)
}

# PLSDAbatch expects samples in rows and ASVs in columns.
otu_mat <- t(as.matrix(mt_raw$otu_table))
storage.mode(otu_mat) <- "numeric"

meta <- as.data.frame(mt_raw$sample_table)
meta <- meta[rownames(otu_mat), , drop = FALSE]

stopifnot(identical(rownames(meta), rownames(otu_mat)))

keep_samples <- meta$sample_type %in% c("Digesta", "Mucus") &
  meta$Diet %in% c("NC", "HC")

meta <- meta[keep_samples, , drop = FALSE]
otu_mat <- otu_mat[rownames(meta), , drop = FALSE]

meta$Diet <- factor(meta$Diet, levels = c("NC", "HC"))

stopifnot(identical(rownames(meta), rownames(otu_mat)))

# Remove unused factor levels after filtering
meta$sample_type <- droplevels(meta$sample_type)
meta$Diet <- droplevels(factor(meta$Diet, levels = c("NC", "HC")))

# Confirm that only the intended samples remain
stopifnot(all(as.character(meta$sample_type) %in% c("Digesta", "Mucus")))
stopifnot(all(as.character(meta$Diet) %in% c("NC", "HC")))

# Print sample numbers within each compartment
cat("\nSamples retained by compartment and diet:\n")

sample_counts <- as.data.frame.matrix(
  table(meta$sample_type, meta$Diet)
)

sample_counts$Total <- rowSums(sample_counts)
print(sample_counts)


# 2. Define biological and technical variables; split compartments

# Change this only if the sampling-day variable has another name.
day_var <- "Day"

if (!day_var %in% names(meta)) {
  stop(
    "Column '", day_var, "' was not found. Available metadata columns are:\n",
    paste(names(meta), collapse = ", ")
  )
}

batch_vars <- c(
  "DNA_extraction_batch",
  "grinding_method",
  "PCR1_batch",
  "PGTB_batch"
)

# Labels used in the combined pRDA plot and diagnostic-panel headings.
batch_labels <- c(
  DNA_extraction_batch = "DNA extraction date",
  grinding_method = "Mechanical lysis device",
  PCR1_batch = "PCR1 batch",
  PGTB_batch = "Illumina sequencing plate"
)

clean_batch_cols <- function(df, variables) {

  missing_vars <- setdiff(variables, names(df))

  if (length(missing_vars) > 0L) {
    stop("Missing technical variables: ", paste(missing_vars, collapse = ", "))
  }

  for (variable in variables) {
    df[[variable]] <- factor(trimws(as.character(df[[variable]])))
  }

  df
}

meta_digesta <- meta[meta$sample_type == "Digesta", , drop = FALSE]
meta_mucus <- meta[meta$sample_type == "Mucus", , drop = FALSE]

otu_digesta <- otu_mat[rownames(meta_digesta), , drop = FALSE]
otu_mucus <- otu_mat[rownames(meta_mucus), , drop = FALSE]

meta_digesta <- clean_batch_cols(meta_digesta, batch_vars)
meta_mucus <- clean_batch_cols(meta_mucus, batch_vars)


# Pre-filter raw ASV counts and CLR-transform each compartment

pf_digesta <- PreFL(data = otu_digesta, keep.spl = 10, keep.var = 0.001)
pf_mucus <- PreFL(data = otu_mucus, keep.spl = 10, keep.var = 0.001)

cat(
  "\nDigesta ASVs: ", ncol(otu_digesta), " before filtering; ",
  ncol(pf_digesta$data.filter), " after filtering.\n",
  sep = ""
)

cat(
  "Mucus ASVs: ", ncol(otu_mucus), " before filtering; ",
  ncol(pf_mucus$data.filter), " after filtering.\n",
  sep = ""
)

X_digesta <- logratio.transfo(
  X = pf_digesta$data.filter,
  logratio = "CLR",
  offset = 1
)

X_mucus <- logratio.transfo(
  X = pf_mucus$data.filter,
  logratio = "CLR",
  offset = 1
)

class(X_digesta) <- "matrix"
class(X_mucus) <- "matrix"


# plotting helpers
# Draw a list of ggplot/grid objects as one page.
draw_panel <- function(grobs, ncol = 2, top = NULL) {

  panel <- gridExtra::arrangeGrob(
    grobs = grobs,
    ncol = ncol,
    top = top
  )

  grid::grid.newpage()
  grid::grid.draw(panel)

  invisible(panel)
}

# Identify the PC1 ASV with the largest absolute loading.
get_top_pc1_asv <- function(pca_object) {

  top_asv <- unlist(
    mixOmics::selectVar(pca_object, comp = 1)$name,
    use.names = FALSE
  )[1]

  # Fallback for mixOmics versions returning a nested selectVar() object.
  if (
    length(top_asv) == 0L ||
    is.na(top_asv) ||
    !top_asv %in% rownames(pca_object$loadings$X)
  ) {
    top_asv <- names(
      sort(abs(pca_object$loadings$X[, 1]), decreasing = TRUE)
    )[1]
  }

  top_asv
}

# One compact batch-QC workflow per compartment

run_compartment_batch_qc <- function(
    X,
    meta_df,
    compartment,
    outdir,
    day_var,
    batch_vars,
    batch_labels
) {

  stopifnot(identical(rownames(X), rownames(meta_df)))

  # PCA as in the tutorial-style Scatter_Density panels

  pca_before <- mixOmics::pca(X, ncomp = 3, scale = TRUE)
  top_asv <- get_top_pc1_asv(pca_before)

  # One tutorial-style PCA density panel is produced for each technical factor.
  pca_density_panels <- list()
  boxplots <- list()
  densityplots <- list()
  pRDA_plot_rows <- list()
  report_lines <- c(
    paste0("Compartment: ", compartment),
    paste0("Top PC1 ASV: ", top_asv),
    "Biological covariate in linear models and pRDA:",
    "Diet-by-Day grouping factor (Diet × Day).",
    "",
    strrep("=", 78)
  )

  for (batch_var in batch_vars) {

    batch_label <- batch_labels[[batch_var]]
    batch_values <- trimws(as.character(meta_df[[batch_var]]))

    # Retain only samples with a valid level for the current technical factor.
    keep_batch <- !is.na(batch_values) & nzchar(batch_values)

    X_batch <- X[keep_batch, , drop = FALSE]
    meta_batch <- meta_df[keep_batch, , drop = FALSE]
    batch <- droplevels(factor(batch_values[keep_batch]))

    if (nlevels(batch) < 2L) {
      warning(
        compartment, " - ", batch_label,
        ": fewer than two levels; analysis skipped."
      )
      next
    }

    diet_batch <- droplevels(
      factor(meta_batch$Diet, levels = c("NC", "HC"))
    )

    day_batch <- droplevels(factor(meta_batch[[day_var]]))

    trt_model <- droplevels(
      interaction(diet_batch, day_batch, sep = "_", drop = TRUE)
    )

    names(batch) <- rownames(meta_batch)
    
    if (all(keep_batch)) {
      pca_for_panel <- pca_before
    } else {
      pca_for_panel <- mixOmics::pca(X_batch, ncomp = 3, scale = TRUE)
    }
    
    # PCA is unsupervised. This panel displays batch only:
    # point colours and marginal densities represent technical batch levels.
    # Diet is intentionally not displayed here, to avoid Batch × Diet density curves.
    pca_density_panels[[batch_var]] <- Scatter_Density(
      object = pca_for_panel,
      batch = batch,
      trt = NULL,
      title = batch_label
    )

    # Top selected-ASV boxplot and density plot

    top_asv_df <- data.frame(
      value = as.numeric(X_batch[, top_asv]),
      batch = batch
    )

    boxplots[[batch_var]] <- box_plot(
      df = top_asv_df,
      title = batch_label,
      x.angle = 30
    )

    densityplots[[batch_var]] <- density_plot(
      df = top_asv_df,
      title = batch_label
    )

    # Tutorial-style linear regression for the selected ASV
    # Use the largest batch as reference. This improves the stability and
    # readability of individual contrasts, especially when a batch is sparse.
    # It does not change fitted values, global batch tests, or pRDA results.
    batch_sizes <- table(batch)
    reference_batch <- names(batch_sizes)[which.max(batch_sizes)]
    batch_relevel <- stats::relevel(batch, ref = reference_batch)

    top_lm <- linear_regres(
      data = X_batch[, top_asv],
      trt = trt_model,
      batch.fix = batch_relevel,
      type = "linear model"
    )

    lm_fit <- top_lm$model$data
    lm_summary <- summary(lm_fit)

    batch_global_test <- drop1(lm_fit, test = "F")
    batch_global_test <- batch_global_test[
      rownames(batch_global_test) == "batch.fix",
      ,
      drop = FALSE
    ]

    # pRDA / variation partitioning
    # [a] = unique Diet × Day biology
    # [b] = unique technical batch effect after adjusting for biology
    # [c] = shared biological/technical variation
    # [d] = residual variation

    factors_df <- data.frame(
      trt = trt_model,
      batch = batch,
      row.names = rownames(X_batch)
    )

    rda_before <- vegan::varpart(
      X_batch,
      ~ trt,
      ~ batch,
      data = factors_df,
      scale = TRUE
    )

    indfract <- rda_before$part$indfract

    prop_raw <- data.frame(
      Treatment = indfract[1, "Adj.R.squared"],
      Intersection = indfract[3, "Adj.R.squared"],
      Batch = indfract[2, "Adj.R.squared"],
      Residuals = indfract[4, "Adj.R.squared"]
    )

    rownames(prop_raw) <- batch_label

    # partVar_plot() needs non-negative fractions summing to one. This
    # transformation is for display only; raw adjusted R² values are retained
    # in the text report and should be used for interpretation.
    prop_plot <- prop_raw
    prop_plot[prop_plot < 0] <- 0
    prop_plot <- prop_plot / sum(unlist(prop_plot))

    pRDA_plot_rows[[batch_var]] <- prop_plot

    # Collect all results in one text report
    report_lines <- c(
      report_lines,
      "",
      paste0("Technical factor: ", batch_label),
      paste0("Selected ASV: ", top_asv),
      paste0("Reference batch: ", reference_batch),
      "Samples per batch:",
      capture.output(table(batch)),
      "",
      "Linear model: CLR abundance ~ Diet-by-Day group + technical batch",
      capture.output(lm_summary),
      "",
      "Global batch test: partial F-test adjusted for Diet × Day",
      capture.output(batch_global_test),
      "",
      "Note: individual batch coefficients compare each level with the",
      "reference batch. Use the global batch test and pRDA fraction to",
      "evaluate the overall technical-factor effect.",
      "",
      "Raw pRDA adjusted R² fractions:",
      capture.output(prop_raw),
      strrep("-", 78)
    )
  }

  # Save one multipage diagnostic PDF

  diagnostics_pdf <- file.path(
    outdir,
    paste0(compartment, "_batch_diagnostics.pdf")
  )

  pdf(diagnostics_pdf, width = 12, height = 8)

  draw_panel(
    grobs = pca_density_panels,
    ncol = 2,
    top = "PCA with technical-factor density overlays"
  )

  draw_panel(
    grobs = boxplots,
    ncol = 2,
    top = paste0("Top PC1 ASV (", top_asv, "): batch-level boxplots")
  )

  draw_panel(
    grobs = densityplots,
    ncol = 2,
    top = paste0("Top PC1 ASV (", top_asv, "): batch-level density plots")
  )

  dev.off()

  # Display the same three pages interactively in RStudio.
  draw_panel(
    grobs = pca_density_panels,
    ncol = 2,
    top = NULL
  )
  draw_panel(
    grobs = boxplots,
    ncol = 2,
    top = paste0("Top PC1 ASV (", top_asv, ")")
  )
  draw_panel(
    grobs = densityplots,
    ncol = 2,
    top = paste0("Top PC1 ASV (", top_asv, ")")
  )

  # Save one pRDA summary plot containing all separate batch analyses
  pRDA_summary_plot <- do.call(rbind, pRDA_plot_rows)

  # Replace internal metadata-column names with descriptive technical labels.
  rownames(pRDA_summary_plot) <- unname(
    batch_labels[names(pRDA_plot_rows)]
  )

  pRDA_plot <- partVar_plot(prop.df = pRDA_summary_plot)

  pRDA_pdf <- file.path(
    outdir,
    paste0(compartment, "_batch_variance_summary.pdf")
  )

  print(pRDA_plot)

  pdf(pRDA_pdf, width = 8, height = 6)
  print(pRDA_plot)
  dev.off()

  # Save one text report with all regression and raw pRDA results

  report_file <- file.path(
    outdir,
    paste0(compartment, "_batch_QC_statistics.txt")
  )

  writeLines(report_lines, report_file)

  invisible(list(
    pca = pca_before,
    top_asv = top_asv,
    pRDA_relative = pRDA_summary_plot
  ))
}

# Run compact QC workflow for Digesta and Mucus

results_digesta <- run_compartment_batch_qc(
  X = X_digesta,
  meta_df = meta_digesta,
  compartment = "Digesta",
  outdir = file.path(output_dir, "Digesta"),
  day_var = day_var,
  batch_vars = batch_vars,
  batch_labels = batch_labels
)

results_mucus <- run_compartment_batch_qc(
  X = X_mucus,
  meta_df = meta_mucus,
  compartment = "Mucus",
  outdir = file.path(output_dir, "Mucus"),
  day_var = day_var,
  batch_vars = batch_vars,
  batch_labels = batch_labels
)

# Save session

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "sessionInfo.txt")
)

