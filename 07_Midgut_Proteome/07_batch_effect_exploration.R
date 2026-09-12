library(ggplot2)
library(dplyr)
library(patchwork)
library(data.table)
library(ggExtra)
library(vegan)

# Output directory
OUTDIR <- "07_batch_effect_exploration"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# INPUTS
INPUT_FILE    <- file.path("rds", "midgut_proteome_filtered.rds")
METADATA_FILE <- "STPN2309_metadata_all.tsv"

SAMPLE_ID_COL <- "sample_id"
DIET_COL      <- "diet"
DAY_COL       <- "day_categ"
TANK_COL      <- "rearing_tank"
BATCH_COL     <- "gut_proteome_batch"

DIETS_TO_KEEP <- c("NC", "HC")

# LOAD MATRIX AND METADATA
proteome_df <- as.data.frame(readRDS(INPUT_FILE))

meta_df <- fread(METADATA_FILE) %>%
  as.data.frame()

rownames(meta_df) <- meta_df[[1]]
meta_df <- meta_df[, -1, drop = FALSE]

if (!SAMPLE_ID_COL %in% colnames(meta_df)) {
  meta_df[[SAMPLE_ID_COL]] <- rownames(meta_df)
}

required_cols <- c(
  SAMPLE_ID_COL,
  DIET_COL,
  DAY_COL,
  TANK_COL,
  BATCH_COL
)

missing_cols <- setdiff(required_cols, colnames(meta_df))

if (length(missing_cols) > 0) {
  stop(
    "Missing metadata column(s): ",
    paste(missing_cols, collapse = ", ")
  )
}

# FILTER AND ALIGN SAMPLES
meta_df <- meta_df %>%
  mutate(
    across(
      all_of(required_cols),
      ~ trimws(as.character(.x))
    )
  ) %>%
  filter(.data[[DIET_COL]] %in% DIETS_TO_KEEP)

common_ids <- intersect(rownames(proteome_df), meta_df[[SAMPLE_ID_COL]])

if (length(common_ids) == 0) {
  stop("No shared sample IDs between proteome matrix and metadata.")
}

common_ids <- sort(common_ids)

proteome_df <- proteome_df[common_ids, , drop = FALSE]

meta_df <- meta_df[
  match(common_ids, meta_df[[SAMPLE_ID_COL]]),
  ,
  drop = FALSE
]

stopifnot(
  identical(rownames(proteome_df), meta_df[[SAMPLE_ID_COL]])
)

# PREPARE MATRIX FOR PCA
proteome_df[] <- lapply(proteome_df, as.numeric)

# Matrix was already filtered upstream at 30% missingness.
# These checks only remove proteins unsuitable for PCA.

# Remove proteins with <= 1 observed value
proteome_df <- proteome_df[
  ,
  colSums(!is.na(proteome_df)) > 1,
  drop = FALSE
]

# Remove constant proteins
is_constant <- vapply(
  proteome_df,
  function(x) {
    x <- x[!is.na(x)]
    length(x) < 2 || sd(x) < 1e-8
  },
  logical(1)
)

proteome_df <- proteome_df[, !is_constant, drop = FALSE]

# Median imputation, only for PCA visualization and distance analyses
for (j in seq_len(ncol(proteome_df))) {
  proteome_df[[j]][is.na(proteome_df[[j]])] <- median(
    proteome_df[[j]],
    na.rm = TRUE
  )
}

# PCA
pca_res <- prcomp(
  proteome_df,
  center = TRUE,
  scale. = TRUE
)

variance_explained <- 100 * (
  pca_res$sdev^2 / sum(pca_res$sdev^2)
)

scores_df <- as.data.frame(pca_res$x) %>%
  mutate(sample_id = rownames(.)) %>%
  left_join(
    meta_df,
    by = setNames(SAMPLE_ID_COL, "sample_id")
  ) %>%
  mutate(
    !!DIET_COL := factor(.data[[DIET_COL]], levels = DIETS_TO_KEEP),
    !!DAY_COL := factor(
      .data[[DAY_COL]],
      levels = c(
        "day_01", "day_02", "day_03", "day_04",
        "day_10", "day_15", "day_22"
      )
    ),
    !!TANK_COL := factor(.data[[TANK_COL]]),
    !!BATCH_COL := factor(.data[[BATCH_COL]])
  )

pc1_lab <- paste0(
  "PC1 (",
  sprintf("%.1f", variance_explained[1]),
  "%)"
)

pc2_lab <- paste0(
  "PC2 (",
  sprintf("%.1f", variance_explained[2]),
  "%)"
)

# PCA + DENSITY PANEL FUNCTION
make_pca_density_plot <- function(scores_df, group_var, group_label) {
  
  p <- ggplot(
    scores_df,
    aes(
      x = PC1,
      y = PC2,
      colour = .data[[group_var]]
    )
  ) +
    geom_point(size = 3, alpha = 0.5) +
    labs(
      title = paste("PCA coloured by", group_label),
      x = pc1_lab,
      y = pc2_lab,
      colour = group_label
    ) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold")
    )
  
  ggExtra::ggMarginal(
    p,
    type = "density",
    groupColour = TRUE,
    groupFill = TRUE,
    alpha = 0.10,
    size = 4
  )
}

p_diet <- make_pca_density_plot(
  scores_df,
  group_var = DIET_COL,
  group_label = "Diet"
)

p_day <- make_pca_density_plot(
  scores_df,
  group_var = DAY_COL,
  group_label = "Day"
)

p_tank <- make_pca_density_plot(
  scores_df,
  group_var = TANK_COL,
  group_label = "Rearing tank"
)

p_batch <- make_pca_density_plot(
  scores_df,
  group_var = BATCH_COL,
  group_label = "Gut proteome batch"
)

p_diet
p_day
p_tank
p_batch

ggsave(
  file.path(OUTDIR, "PCA_marginal_density_by_diet.pdf"),
  plot = p_diet,
  width = 8,
  height = 6
)

ggsave(
  file.path(OUTDIR, "PCA_marginal_density_by_day.pdf"),
  plot = p_day,
  width = 8,
  height = 6
)

ggsave(
  file.path(OUTDIR, "PCA_marginal_density_by_rearing_tank.pdf"),
  plot = p_tank,
  width = 8,
  height = 6
)

ggsave(
  file.path(OUTDIR, "PCA_marginal_density_by_gut_proteome_batch.pdf"),
  plot = p_batch,
  width = 8,
  height = 6
)

# BATCH-EFFECT ASSESSMENT
BATCH_REPORT <- file.path(
  OUTDIR,
  "batch_effect_assessment_midgut_proteome.txt"
)

# Standardize proteins before Euclidean-distance analyses
proteome_for_test <- scale(
  proteome_df,
  center = TRUE,
  scale = TRUE
)

dist_proteome <- dist(
  proteome_for_test,
  method = "euclidean"
)

# 1. Unadjusted batch association
set.seed(123)

permanova_batch <- vegan::adonis2(
  dist_proteome ~ gut_proteome_batch,
  data = meta_df,
  permutations = 9999
)

# 2. Homogeneity of multivariate dispersion
disp_batch <- vegan::betadisper(
  dist_proteome,
  group = meta_df$gut_proteome_batch
)

disp_batch_anova <- anova(disp_batch)

set.seed(123)

disp_batch_perm <- permutest(
  disp_batch,
  permutations = 9999
)

# 3. Independent diet, day, and batch effects
set.seed(123)

permanova_adjusted <- vegan::adonis2(
  dist_proteome ~ diet + day_categ + gut_proteome_batch,
  data = meta_df,
  permutations = 9999,
  by = "margin"
)

# 4. Batch × sampling-day allocation
tab_batch_day <- table(
  meta_df$gut_proteome_batch,
  meta_df$day_categ
)

set.seed(123)

batch_day_test <- chisq.test(
  tab_batch_day,
  simulate.p.value = TRUE,
  B = 9999
)

# 5. Batch × Diet-by-Day allocation
tab_batch_design <- table(
  meta_df$gut_proteome_batch,
  interaction(meta_df$diet, meta_df$day_categ)
)

set.seed(123)

batch_design_test <- chisq.test(
  tab_batch_design,
  simulate.p.value = TRUE,
  B = 9999
)

# 6. Check estimability of batch-adjusted model
design_check <- model.matrix(
  ~ gut_proteome_batch + diet * day_categ,
  data = meta_df
)

design_rank <- c(
  n_model_columns = ncol(design_check),
  model_rank = qr(design_check)$rank
)

# PRINT AND SAVE REPORT
report_lines <- capture.output({
  
  cat("MIDGUT PROTEOME BATCH-EFFECT ASSESSMENT\n")
  cat("========================================\n\n")
  
  cat("Input matrix:\n")
  cat("- Samples:", nrow(proteome_df), "\n")
  cat("- Proteins:", ncol(proteome_df), "\n")
  cat("- Distance: Euclidean on feature-wise standardized protein abundances\n")
  cat("- Permutations: 9,999\n\n")
  
  cat("1. Unadjusted batch association\n")
  cat("--------------------------------\n")
  print(permanova_batch)
  
  cat("\n2. Homogeneity of multivariate dispersion by batch\n")
  cat("---------------------------------------------------\n")
  print(disp_batch_anova)
  
  cat("\nPermutation test:\n")
  print(disp_batch_perm)
  
  cat("\n3. Marginal PERMANOVA: diet + day + batch\n")
  cat("------------------------------------------\n")
  print(permanova_adjusted)
  
  cat("\n4. Batch × sampling-day allocation\n")
  cat("-----------------------------------\n")
  print(tab_batch_day)
  
  cat("\nPermutation chi-square test:\n")
  print(batch_day_test)
  
  cat("\n5. Batch × Diet-by-Day allocation\n")
  cat("----------------------------------\n")
  print(tab_batch_design)
  
  cat("\nPermutation chi-square test:\n")
  print(batch_design_test)
  
  cat("\n6. Estimability of batch-adjusted biological model\n")
  cat("-------------------------------------------------\n")
  print(design_rank)
  
  if (design_rank["model_rank"] == design_rank["n_model_columns"]) {
    cat(
      "\nResult: full-rank model. ",
      "The model gut_proteome_batch + diet * day_categ is estimable ",
      "without exact confounding.\n",
      sep = ""
    )
  } else {
    cat(
      "\nResult: rank-deficient model. ",
      "Some batch and biological effects are exactly confounded.\n",
      sep = ""
    )
  }
  
  cat("\nInterpretation guide:\n")
  cat("- Significant batch PERMANOVA indicates different batch centroids.\n")
  cat("- Non-significant dispersion indicates the batch result is not due to unequal within-batch variability.\n")
  cat("- Marginal PERMANOVA quantifies independent diet, day, and batch contributions.\n")
  cat("- Batch × day testing evaluates whether batches differ in day composition.\n")
  cat("- Batch × Diet × Day testing evaluates imbalance across the full biological design.\n")
  cat("- A full-rank model supports including batch as a fixed covariate:\n")
  cat("  ~ gut_proteome_batch + diet * day_categ\n")
})

cat(paste(report_lines, collapse = "\n"), "\n")

writeLines(
  report_lines,
  con = BATCH_REPORT
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(OUTDIR, "sessionInfo.txt")
)

cat(
  "\nSaved batch-effect report to:\n",
  BATCH_REPORT,
  "\n",
  sep = ""
)
