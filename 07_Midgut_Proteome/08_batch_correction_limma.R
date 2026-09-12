# Libraries
library(limma)
library(mixOmics)
library(ggplot2)
library(dplyr)
library(patchwork)

# Settings
INPUT_MATRIX_RDS <- file.path(
  "rds",
  "midgut_proteome_filtered.rds"
)

METADATA_FILE <- "STPN2309_metadata_all.tsv"

OUTPUT_DIR <- "08_batch_correction_limma"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

SAMPLE_ID_COL <- "sampleID_STPN2309"
DIET_COL      <- "diet"
DAY_COL       <- "day_categ"
BATCH_COL     <- "gut_proteome_batch"

DIETS_TO_KEEP <- c("NC", "HC")

CENTER_DATA <- TRUE
SCALE_DATA  <- TRUE
BASE_SIZE   <- 14

# Helper functions
theme_clean <- function() {
  theme_bw(base_size = BASE_SIZE) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      plot.title.position = "plot"
    )
}

theme_density <- function() {
  theme_bw(base_size = BASE_SIZE) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      plot.title.position = "plot",
      legend.position = "right"
    )
}

clean_feature_names <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\\s+", "", x)
  make.unique(x)
}

print_matrix_summary <- function(mat, label = "Matrix") {
  
  vals <- as.numeric(as.matrix(mat))
  vals <- vals[!is.na(vals)]
  
  cat("\n=== ", label, " ===\n", sep = "")
  cat("Dimensions:", nrow(mat), "samples x", ncol(mat), "proteins\n")
  cat("Observed values:", length(vals), "\n")
  cat("Missing values:", sum(is.na(mat)), "\n")
  
  if (length(vals) > 0) {
    print(summary(vals))
    cat("Minimum:", min(vals), "\n")
    cat("Maximum:", max(vals), "\n")
  }
  
  cat("\n")
}

get_aov_p <- function(formula_obj, data_obj) {
  fit <- aov(formula_obj, data = data_obj)
  summary(fit)[[1]][["Pr(>F)"]][1]
}

make_pca_density_panel <- function(
    scores_df,
    pct_var,
    batch_col,
    diet_col,
    title_text
) {
  
  scatter_plot <- ggplot(
    scores_df,
    aes(
      x = PC1,
      y = PC2,
      color = .data[[batch_col]],
      shape = .data[[diet_col]]
    )
  ) +
    geom_point(size = 3, alpha = 0.9) +
    labs(
      title = title_text,
      x = paste0("PC1 (", sprintf("%.1f", pct_var[1]), "%)"),
      y = paste0("PC2 (", sprintf("%.1f", pct_var[2]), "%)"),
      color = "Batch",
      shape = "Diet"
    ) +
    theme_clean()
  
  pc1_density <- ggplot(
    scores_df,
    aes(
      x = PC1,
      fill = .data[[batch_col]],
      color = .data[[batch_col]]
    )
  ) +
    geom_density(alpha = 0.30, linewidth = 0.7) +
    labs(
      title = "PC1 density by batch",
      x = paste0("PC1 (", sprintf("%.1f", pct_var[1]), "%)"),
      y = "Density",
      fill = "Batch",
      color = "Batch"
    ) +
    theme_density()
  
  pc2_density <- ggplot(
    scores_df,
    aes(
      x = PC2,
      fill = .data[[batch_col]],
      color = .data[[batch_col]]
    )
  ) +
    geom_density(alpha = 0.30, linewidth = 0.7) +
    labs(
      title = "PC2 density by batch",
      x = paste0("PC2 (", sprintf("%.1f", pct_var[2]), "%)"),
      y = "Density",
      fill = "Batch",
      color = "Batch"
    ) +
    theme_density()
  
  (scatter_plot | pc1_density | pc2_density) +
    plot_layout(widths = c(1.2, 1, 1))
}

get_batch_p <- function(mat, meta) {
  
  sapply(colnames(mat), function(feature) {
    
    df <- data.frame(
      y = mat[[feature]],
      meta,
      check.names = FALSE
    )
    
    df <- df[
      complete.cases(df[, c("y", DIET_COL, DAY_COL, BATCH_COL)]),
      ,
      drop = FALSE
    ]
    
    if (nrow(df) < 5) {
      return(NA_real_)
    }
    
    fit <- lm(
      as.formula(
        paste(
          "y ~",
          DIET_COL,
          "*",
          DAY_COL,
          "+",
          BATCH_COL
        )
      ),
      data = df
    )
    
    an <- anova(fit)
    
    if (!BATCH_COL %in% rownames(an)) {
      return(NA_real_)
    }
    
    an[BATCH_COL, "Pr(>F)"]
  })
}

get_diet_delta <- function(mat, meta) {
  
  sapply(colnames(mat), function(feature) {
    
    mean(
      mat[meta[[DIET_COL]] == "HC", feature],
      na.rm = TRUE
    ) -
      mean(
        mat[meta[[DIET_COL]] == "NC", feature],
        na.rm = TRUE
      )
  })
}

# Check inputs
if (!file.exists(INPUT_MATRIX_RDS)) {
  stop("Proteome matrix RDS not found: ", INPUT_MATRIX_RDS)
}

if (!file.exists(METADATA_FILE)) {
  stop("Metadata TSV not found: ", METADATA_FILE)
}

# Load proteome matrix and metadata
proteome_df <- readRDS(INPUT_MATRIX_RDS)
proteome_df <- as.data.frame(proteome_df, stringsAsFactors = FALSE)

meta_df <- read.delim(
  METADATA_FILE,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (is.null(rownames(proteome_df))) {
  stop("The proteome matrix must have sample IDs stored as rownames.")
}

required_cols <- c(
  SAMPLE_ID_COL,
  DIET_COL,
  DAY_COL,
  BATCH_COL
)

missing_cols <- setdiff(required_cols, names(meta_df))

if (length(missing_cols) > 0) {
  stop(
    "Missing required metadata column(s): ",
    paste(missing_cols, collapse = ", ")
  )
}

# Clean metadata and IDs
rownames(proteome_df) <- trimws(rownames(proteome_df))

meta_df[[SAMPLE_ID_COL]] <- trimws(
  as.character(meta_df[[SAMPLE_ID_COL]])
)

meta_df[[DIET_COL]] <- trimws(
  as.character(meta_df[[DIET_COL]])
)

meta_df[[DAY_COL]] <- trimws(
  as.character(meta_df[[DAY_COL]])
)

meta_df[[BATCH_COL]] <- trimws(
  as.character(meta_df[[BATCH_COL]])
)

if (anyDuplicated(rownames(proteome_df)) > 0) {
  stop("Duplicated sample IDs found in proteome matrix rownames.")
}

if (anyDuplicated(meta_df[[SAMPLE_ID_COL]]) > 0) {
  stop(
    "Duplicated sample IDs found in metadata column: ",
    SAMPLE_ID_COL
  )
}

meta_df <- meta_df[
  meta_df[[DIET_COL]] %in% DIETS_TO_KEEP,
  ,
  drop = FALSE
]

# Align proteome matrix and metadata
common_ids <- intersect(
  rownames(proteome_df),
  meta_df[[SAMPLE_ID_COL]]
)

if (length(common_ids) == 0) {
  stop(
    "No common sample IDs found between proteome matrix rownames and metadata."
  )
}

common_ids <- rownames(proteome_df)[
  rownames(proteome_df) %in% common_ids
]

proteome_df <- proteome_df[
  common_ids,
  ,
  drop = FALSE
]

meta_df <- meta_df[
  match(common_ids, meta_df[[SAMPLE_ID_COL]]),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(proteome_df),
    meta_df[[SAMPLE_ID_COL]]
  )
)

cat("Aligned data:\n")
cat(
  "Proteome matrix:",
  nrow(proteome_df),
  "samples x",
  ncol(proteome_df),
  "proteins\n"
)
cat("Metadata:", nrow(meta_df), "samples\n")

# Clean proteome matrix
colnames(proteome_df) <- clean_feature_names(colnames(proteome_df))

proteome_df[] <- lapply(
  proteome_df,
  function(x) suppressWarnings(as.numeric(x))
)

# The input was already filtered at 30% missingness.
# These are only safety checks for unusable columns.

keep_nonmissing <- colSums(!is.na(proteome_df)) > 1

if (any(!keep_nonmissing)) {
  
  message(
    "Dropping ",
    sum(!keep_nonmissing),
    " protein(s) with <= 1 observed value."
  )
  
  proteome_df <- proteome_df[
    ,
    keep_nonmissing,
    drop = FALSE
  ]
}

is_constant <- vapply(
  proteome_df,
  function(x) {
    x <- x[!is.na(x)]
    length(x) < 2 || sd(x) < 1e-8
  },
  logical(1)
)

if (any(is_constant)) {
  
  message(
    "Dropping ",
    sum(is_constant),
    " constant protein(s)."
  )
  
  proteome_df <- proteome_df[
    ,
    !is_constant,
    drop = FALSE
  ]
}

if (ncol(proteome_df) < 2) {
  stop("Fewer than two proteins remain after filtering.")
}

print_matrix_summary(
  proteome_df,
  "Proteome matrix before batch correction"
)

# Prepare factors
meta_df[[DIET_COL]] <- factor(
  meta_df[[DIET_COL]],
  levels = DIETS_TO_KEEP
)

meta_df[[DAY_COL]] <- factor(
  meta_df[[DAY_COL]],
  levels = c(
    "day_01", "day_02", "day_03", "day_04",
    "day_10", "day_15", "day_22"
  )
)

meta_df[[BATCH_COL]] <- factor(meta_df[[BATCH_COL]])

# Batch balance diagnostics
tab_batch_diet <- table(
  meta_df[[BATCH_COL]],
  meta_df[[DIET_COL]]
)

tab_batch_day <- table(
  meta_df[[BATCH_COL]],
  meta_df[[DAY_COL]]
)

prop_batch_diet <- prop.table(
  tab_batch_diet,
  margin = 1
)

prop_batch_day <- prop.table(
  tab_batch_day,
  margin = 1
)

write.csv(
  as.data.frame.matrix(tab_batch_diet),
  file.path(OUTPUT_DIR, "batch_by_diet_counts.csv")
)

write.csv(
  as.data.frame.matrix(round(prop_batch_diet, 4)),
  file.path(OUTPUT_DIR, "batch_by_diet_row_proportions.csv")
)

write.csv(
  as.data.frame.matrix(tab_batch_day),
  file.path(OUTPUT_DIR, "batch_by_day_counts.csv")
)

write.csv(
  as.data.frame.matrix(round(prop_batch_day, 4)),
  file.path(OUTPUT_DIR, "batch_by_day_row_proportions.csv")
)

# limma batch correction
# Preserve Diet, Day, and Diet × Day effects
proteome_limma <- proteome_df
meta_limma <- meta_df

design <- model.matrix(
  as.formula(
    paste("~", DIET_COL, "*", DAY_COL)
  ),
  data = meta_limma
)

mat_for_limma <- t(as.matrix(proteome_limma))

mat_corrected <- limma::removeBatchEffect(
  x = mat_for_limma,
  batch = meta_limma[[BATCH_COL]],
  design = design
)

proteome_corrected <- as.data.frame(t(mat_corrected))

rownames(proteome_corrected) <- rownames(proteome_limma)
colnames(proteome_corrected) <- colnames(proteome_limma)

stopifnot(
  identical(rownames(proteome_limma), rownames(proteome_corrected)),
  identical(colnames(proteome_limma), colnames(proteome_corrected))
)

write.csv(
  cbind(
    sampleID_STPN2309 = rownames(proteome_corrected),
    proteome_corrected
  ),
  file.path(
    OUTPUT_DIR,
    "midgut_proteome_filtered_batch_corrected_limma.csv"
  ),
  row.names = FALSE
)

# PCA before and after correction
pca_before <- mixOmics::pca(
  X = proteome_limma,
  ncomp = 5,
  center = CENTER_DATA,
  scale = SCALE_DATA
)

pca_after <- mixOmics::pca(
  X = proteome_corrected,
  ncomp = 5,
  center = CENTER_DATA,
  scale = SCALE_DATA
)

scores_before <- as.data.frame(pca_before$variates$X)
scores_after <- as.data.frame(pca_after$variates$X)

pct_before <- 100 * pca_before$prop_expl_var$X
pct_after <- 100 * pca_after$prop_expl_var$X

scores_before[[SAMPLE_ID_COL]] <- rownames(scores_before)

scores_before <- left_join(
  scores_before,
  meta_limma,
  by = SAMPLE_ID_COL
)

scores_after[[SAMPLE_ID_COL]] <- rownames(scores_after)

scores_after <- left_join(
  scores_after,
  meta_limma,
  by = SAMPLE_ID_COL
)

write.csv(
  scores_before,
  file.path(
    OUTPUT_DIR,
    "PCA_scores_before_batch_correction.csv"
  ),
  row.names = FALSE
)

write.csv(
  scores_after,
  file.path(
    OUTPUT_DIR,
    "PCA_scores_after_batch_correction.csv"
  ),
  row.names = FALSE
)

# PC association tests
pc_assoc_before <- data.frame(
  variable = c(BATCH_COL, DIET_COL, DAY_COL),
  p_PC1 = c(
    get_aov_p(as.formula(paste("PC1 ~", BATCH_COL)), scores_before),
    get_aov_p(as.formula(paste("PC1 ~", DIET_COL)), scores_before),
    get_aov_p(as.formula(paste("PC1 ~", DAY_COL)), scores_before)
  ),
  p_PC2 = c(
    get_aov_p(as.formula(paste("PC2 ~", BATCH_COL)), scores_before),
    get_aov_p(as.formula(paste("PC2 ~", DIET_COL)), scores_before),
    get_aov_p(as.formula(paste("PC2 ~", DAY_COL)), scores_before)
  )
)

pc_assoc_after <- data.frame(
  variable = c(BATCH_COL, DIET_COL, DAY_COL),
  p_PC1 = c(
    get_aov_p(as.formula(paste("PC1 ~", BATCH_COL)), scores_after),
    get_aov_p(as.formula(paste("PC1 ~", DIET_COL)), scores_after),
    get_aov_p(as.formula(paste("PC1 ~", DAY_COL)), scores_after)
  ),
  p_PC2 = c(
    get_aov_p(as.formula(paste("PC2 ~", BATCH_COL)), scores_after),
    get_aov_p(as.formula(paste("PC2 ~", DIET_COL)), scores_after),
    get_aov_p(as.formula(paste("PC2 ~", DAY_COL)), scores_after)
  )
)

write.csv(
  pc_assoc_before,
  file.path(OUTPUT_DIR, "PCA_PC_association_tests_before.csv"),
  row.names = FALSE
)

write.csv(
  pc_assoc_after,
  file.path(OUTPUT_DIR, "PCA_PC_association_tests_after.csv"),
  row.names = FALSE
)

# PCA plots
p_before_density <- make_pca_density_panel(
  scores_df = scores_before,
  pct_var = pct_before,
  batch_col = BATCH_COL,
  diet_col = DIET_COL,
  title_text = "PCA before batch correction"
)

p_after_density <- make_pca_density_panel(
  scores_df = scores_after,
  pct_var = pct_after,
  batch_col = BATCH_COL,
  diet_col = DIET_COL,
  title_text = "PCA after batch correction"
)

ggsave(
  file.path(
    OUTPUT_DIR,
    "PCA_with_density_before_batch_correction.pdf"
  ),
  p_before_density,
  width = 15,
  height = 5
)

ggsave(
  file.path(
    OUTPUT_DIR,
    "PCA_with_density_before_batch_correction.tiff"
  ),
  p_before_density,
  width = 15,
  height = 5,
  dpi = 1200,
  compression = "lzw"
)

ggsave(
  file.path(
    OUTPUT_DIR,
    "PCA_with_density_after_batch_correction.pdf"
  ),
  p_after_density,
  width = 15,
  height = 5
)

ggsave(
  file.path(
    OUTPUT_DIR,
    "PCA_with_density_after_batch_correction.tiff"
  ),
  p_after_density,
  width = 15,
  height = 5,
  dpi = 1200,
  compression = "lzw"
)

# Feature-level batch-effect diagnostics
batch_p_before <- get_batch_p(
  proteome_limma,
  meta_limma
)

batch_p_after <- get_batch_p(
  proteome_corrected,
  meta_limma
)

batch_check <- data.frame(
  protein = colnames(proteome_limma),
  batch_p_before = batch_p_before,
  batch_p_after = batch_p_after,
  stringsAsFactors = FALSE
)

batch_check$batch_FDR_before <- p.adjust(
  batch_check$batch_p_before,
  method = "BH"
)

batch_check$batch_FDR_after <- p.adjust(
  batch_check$batch_p_after,
  method = "BH"
)

write.csv(
  batch_check,
  file.path(
    OUTPUT_DIR,
    "batch_pvalues_before_after_limma.csv"
  ),
  row.names = FALSE
)

cat(
  "\nBatch-associated proteins before correction, nominal p < 0.05:",
  sum(batch_check$batch_p_before < 0.05, na.rm = TRUE),
  "\n"
)

cat(
  "Batch-associated proteins after correction, nominal p < 0.05:",
  sum(batch_check$batch_p_after < 0.05, na.rm = TRUE),
  "\n"
)

cat(
  "Batch-associated proteins before correction, BH FDR < 0.05:",
  sum(batch_check$batch_FDR_before < 0.05, na.rm = TRUE),
  "\n"
)

cat(
  "Batch-associated proteins after correction, BH FDR < 0.05:",
  sum(batch_check$batch_FDR_after < 0.05, na.rm = TRUE),
  "\n"
)

# Check preservation of HC - NC differences
diet_delta_before <- get_diet_delta(
  proteome_limma,
  meta_limma
)

diet_delta_after <- get_diet_delta(
  proteome_corrected,
  meta_limma
)

diet_check <- data.frame(
  protein = colnames(proteome_limma),
  diet_delta_before = diet_delta_before,
  diet_delta_after = diet_delta_after,
  stringsAsFactors = FALSE
)

diet_effect_correlation <- cor(
  diet_check$diet_delta_before,
  diet_check$diet_delta_after,
  use = "complete.obs"
)

cat(
  "\nCorrelation of HC - NC protein differences before versus after correction:",
  round(diet_effect_correlation, 4),
  "\n"
)

write.csv(
  diet_check,
  file.path(
    OUTPUT_DIR,
    "diet_delta_before_after_limma.csv"
  ),
  row.names = FALSE
)

# Save objects
saveRDS(
  list(
    proteome_original = proteome_limma,
    proteome_batch_corrected_limma = proteome_corrected,
    metadata = meta_limma,
    pca_before = pca_before,
    pca_after = pca_after,
    scores_before = scores_before,
    scores_after = scores_after,
    batch_check = batch_check,
    diet_check = diet_check
  ),
  file = file.path(
    OUTPUT_DIR,
    "midgut_proteome_filtered_batch_corrected_limma.rds"
  )
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(OUTPUT_DIR, "sessionInfo.txt")
)

cat("\nFinal dimensions:\n")
cat("Original matrix:", dim(proteome_limma), "\n")
cat("Corrected matrix:", dim(proteome_corrected), "\n")

cat(
  "Missing values before correction:",
  sum(is.na(proteome_limma)),
  "\n"
)

cat(
  "Missing values after correction:",
  sum(is.na(proteome_corrected)),
  "\n"
)

cat(
  "\nSaved outputs to:",
  normalizePath(OUTPUT_DIR),
  "\n"
)

# quick check
sum(batch_check$batch_FDR_before < 0.05, na.rm = TRUE)
sum(batch_check$batch_FDR_after < 0.05, na.rm = TRUE)
diet_effect_correlation