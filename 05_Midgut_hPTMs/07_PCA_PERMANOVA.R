library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

## ---------------------------------------------------------
## SETTINGS
## ---------------------------------------------------------
INPUT_RDS <- "05_batch_correction_limma/H3_H4_5PTM_midgut_batch_corrected_limma.rds"
METADATA_FILE <- "STPN2309_metadata_all.tsv"

OUTPUT_DIR <- "07_PCA_PERMANOVA"
DIET_OUTPUT_DIR <- file.path(OUTPUT_DIR, "by_diet")

dir.create(DIET_OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

SAMPLE_ID_COL <- "sampleID_STPN2309"
DIET_COL      <- "diet"
GROUP_A <- "NC"
GROUP_B <- "HC"
DIETS_TO_KEEP <- c(GROUP_A, GROUP_B)

SEED   <- 123
N_PERM <- 9999

BASE_SIZE <- 14
TITLE_FAMILY <- "mono"
TITLE_SIZE   <- 16

cols2diet <- c(
  "NC" = "#44a7c4",
  "HC" = "#f99943"
)

APPLY_Z_CAP <- FALSE
Z_CAP <- 4

## ---------------------------------------------------------
## HELPER FUNCTIONS
## ---------------------------------------------------------
clean_feature_names <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("\\s+", "", x)
  make.unique(x)
}

theme_cutler <- function() {
  theme_bw(base_size = BASE_SIZE) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(
        family = TITLE_FAMILY,
        size = TITLE_SIZE,
        face = "bold"
      ),
      plot.title.position = "plot",
      plot.margin = margin(t = 10, r = 5, b = 5, l = 5)
    )
}

p_to_symbol <- function(p) {
  if (is.na(p)) return("NA")
  if (p <= 1e-4) return("****")
  if (p <= 1e-3) return("***")
  if (p <= 1e-2) return("**")
  if (p <= 0.05) return("*")
  "ns"
}

format_p_value_plot <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 1e-4) return("p < 1e-4")
  if (p < 0.001) return(paste0("p = ", formatC(p, format = "e", digits = 2)))
  paste0("p = ", signif(p, 3))
}

print_matrix_summary <- function(mat, label = "Matrix") {
  vals <- as.numeric(as.matrix(mat))
  vals <- vals[!is.na(vals)]
  
  cat("=== ", label, " ===\n", sep = "")
  cat("Dimensions:", nrow(mat), "rows x", ncol(mat), "cols\n")
  cat("Observed values:", length(vals), "\n")
  cat("Missing values :", sum(is.na(mat)), "\n")
  
  if (length(vals) > 0) {
    print(summary(vals))
    cat("Min:", min(vals), "\n")
    cat("Max:", max(vals), "\n")
  }
  
  cat("\n")
}

## ---------------------------------------------------------
## LOAD LIMMA-CORRECTED hPTM MATRIX AND METADATA
## ---------------------------------------------------------

cat("Loading limma-corrected hPTM matrix and metadata...\n")

if (!file.exists(INPUT_RDS)) {
  stop("Corrected hPTM RDS not found: ", INPUT_RDS)
}

if (!file.exists(METADATA_FILE)) {
  stop("Metadata TSV not found: ", METADATA_FILE)
}

rds_obj <- readRDS(INPUT_RDS)

if (!"ptm_batch_corrected_limma" %in% names(rds_obj)) {
  stop(
    "Object 'ptm_batch_corrected_limma' not found in: ",
    INPUT_RDS
  )
}

ptm_df <- as.data.frame(
  rds_obj$ptm_batch_corrected_limma,
  stringsAsFactors = FALSE
)

meta_df <- read.delim(
  METADATA_FILE,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (is.null(rownames(ptm_df))) {
  stop("Corrected hPTM matrix must have sample IDs as rownames.")
}

required_cols <- c(
  SAMPLE_ID_COL,
  DIET_COL
)

missing_cols <- setdiff(required_cols, names(meta_df))

if (length(missing_cols) > 0) {
  stop(
    "Missing required metadata column(s): ",
    paste(missing_cols, collapse = ", ")
  )
}

cat(
  "Corrected hPTM matrix:",
  nrow(ptm_df), "samples x",
  ncol(ptm_df), "features\n"
)

cat(
  "Metadata:",
  nrow(meta_df), "samples x",
  ncol(meta_df), "columns\n\n"
)

# Keep NC and HC, then align matrix and metadata by sample ID

rownames(ptm_df) <- trimws(rownames(ptm_df))

meta_df[[SAMPLE_ID_COL]] <- trimws(
  as.character(meta_df[[SAMPLE_ID_COL]])
)

meta_df[[DIET_COL]] <- trimws(
  as.character(meta_df[[DIET_COL]])
)

meta_df <- meta_df %>%
  filter(.data[[DIET_COL]] %in% DIETS_TO_KEEP) %>%
  droplevels()

common_ids <- rownames(ptm_df)[
  rownames(ptm_df) %in% meta_df[[SAMPLE_ID_COL]]
]

if (length(common_ids) == 0) {
  stop("No matching sample IDs between hPTM matrix and metadata.")
}

ptm_df <- ptm_df[
  common_ids,
  ,
  drop = FALSE
]

meta_df <- meta_df[
  match(common_ids, meta_df[[SAMPLE_ID_COL]]),
  ,
  drop = FALSE
]

rownames(meta_df) <- meta_df[[SAMPLE_ID_COL]]

meta_df[[DIET_COL]] <- factor(
  meta_df[[DIET_COL]],
  levels = c(GROUP_A, GROUP_B)
)

stopifnot(identical(rownames(ptm_df), rownames(meta_df)))

cat("After diet filtering:\n")
cat("hPTM matrix:", nrow(ptm_df), "samples x", ncol(ptm_df), "features\n")
cat(GROUP_A, "samples:", sum(meta_df[[DIET_COL]] == GROUP_A), "\n")
cat(GROUP_B, "samples:", sum(meta_df[[DIET_COL]] == GROUP_B), "\n\n")

if (sum(meta_df[[DIET_COL]] == GROUP_A) < 4 || sum(meta_df[[DIET_COL]] == GROUP_B) < 4) {
  stop("Need at least 4 samples per diet group.")
}

## ---------------------------------------------------------
## CLEAN hPTM MATRIX
## ---------------------------------------------------------
cat("Cleaning hPTM matrix...\n")

colnames(ptm_df) <- clean_feature_names(colnames(ptm_df))
ptm_df[] <- lapply(ptm_df, function(x) suppressWarnings(as.numeric(x)))

keep_nonmissing <- colSums(!is.na(ptm_df)) > 1

if (any(!keep_nonmissing)) {
  message("Dropping ", sum(!keep_nonmissing), " hPTM(s) with <=1 non-missing value globally.")
  ptm_df <- ptm_df[, keep_nonmissing, drop = FALSE]
}

is_const_global <- vapply(
  ptm_df,
  function(v) {
    v <- v[!is.na(v)]
    length(v) < 2 || sd(v) < 1e-8
  },
  logical(1)
)

if (any(is_const_global)) {
  message("Dropping ", sum(is_const_global), " constant hPTM(s) globally.")
  ptm_df <- ptm_df[, !is_const_global, drop = FALSE]
}

if (ncol(ptm_df) < 2) {
  stop("Not enough hPTMs remain after filtering.")
}

print_matrix_summary(
  ptm_df,
  "hPTM matrix before z-scoring"
)


## ---------------------------------------------------------
## Z-SCORE hPTMs ACROSS SAMPLES
## ---------------------------------------------------------
cat("Z-scoring hPTMs across samples...\n")

ptm_scaled <- scale(ptm_df)
ptm_scaled <- as.data.frame(ptm_scaled)

rownames(ptm_scaled) <- rownames(ptm_df)
colnames(ptm_scaled) <- colnames(ptm_df)

cat("Sanity check after z-scoring:\n")

col_means <- colMeans(ptm_scaled, na.rm = TRUE)
col_sds   <- apply(ptm_scaled, 2, sd, na.rm = TRUE)

cat("Summary of column means:\n")
print(summary(col_means))

cat("\nSummary of column SDs:\n")
print(summary(col_sds))

cat("\nLargest absolute mean deviation:", max(abs(col_means), na.rm = TRUE), "\n")
cat("Range of SD values:", paste(range(col_sds, na.rm = TRUE), collapse = " to "), "\n\n")

if (APPLY_Z_CAP) {
  ptm_scaled_mat <- as.matrix(ptm_scaled)
  
  ptm_scaled_mat[ptm_scaled_mat > Z_CAP] <- Z_CAP
  ptm_scaled_mat[ptm_scaled_mat < -Z_CAP] <- -Z_CAP
  
  ptm_scaled <- as.data.frame(ptm_scaled_mat)
  rownames(ptm_scaled) <- rownames(ptm_df)
  colnames(ptm_scaled) <- colnames(ptm_df)
  
  cat("Applied z-score cap at ±", Z_CAP, "\n\n", sep = "")
} else {
  cat("No z-score cap applied. PCA/PERMANOVA uses pure z-scores.\n\n")
}

print_matrix_summary(
  ptm_scaled,
  "Z-scored hPTM matrix"
)

## PCA and Euclidean distance cannot handle missing values.
## Since data are z-scored, replacing NA by 0 imputes the feature mean.
ptm_scaled_pca <- ptm_scaled
ptm_scaled_pca[is.na(ptm_scaled_pca)] <- 0


## ---------------------------------------------------------
## METADATA FOR PCA / PERMANOVA
## ---------------------------------------------------------
meta_test <- data.frame(
  sample_id = rownames(meta_df),
  diet = factor(meta_df[[DIET_COL]], levels = c(GROUP_A, GROUP_B))
)

rownames(meta_test) <- meta_test$sample_id

stopifnot(identical(rownames(ptm_scaled_pca), rownames(meta_test)))


## ---------------------------------------------------------
## PERMANOVA
## ---------------------------------------------------------
cat("Running PERMANOVA on z-scored hPTM profiles...\n")

dist_mat <- dist(ptm_scaled_pca, method = "euclidean")

set.seed(SEED)

permanova_res <- adonis2(
  dist_mat ~ diet,
  data = meta_test,
  permutations = N_PERM
)

permanova_tab <- as.data.frame(permanova_res)
permanova_tab$term <- rownames(permanova_tab)
rownames(permanova_tab) <- NULL

diet_row <- permanova_tab[
  !(permanova_tab$term %in% c("Residual", "Total")),
  ,
  drop = FALSE
]

if (nrow(diet_row) < 1) {
  stop("Could not extract the PERMANOVA diet term.")
}

diet_row <- diet_row[1, , drop = FALSE]
diet_row$term <- "diet"

permanova_R2 <- diet_row$R2
permanova_F  <- diet_row$F
permanova_p  <- diet_row$`Pr(>F)`



permanova_summary <- data.frame(
  test = "PERMANOVA",
  omics = "hPTM",
  factor = "diet",
  df = diet_row$Df,
  sum_of_squares = diet_row$SumOfSqs,
  R2 = permanova_R2,
  F_value = permanova_F,
  p_value = permanova_p,
  significant = permanova_p < 0.05
)

write.csv(
  permanova_summary,
  file.path(DIET_OUTPUT_DIR, "hPTM_PERMANOVA_summary.csv"),
  row.names = FALSE
)

print(permanova_res)

## ---------------------------------------------------------
## BETADISPER
## ---------------------------------------------------------
cat("Checking homogeneity of multivariate dispersion with betadisper...\n")

bd <- betadisper(dist_mat, meta_test$diet)

bd_anova <- anova(bd)

set.seed(SEED)
bd_perm <- permutest(bd, permutations = N_PERM)



print(bd_anova)
print(bd_perm)

plot(bd, main = "hPTM betadisper")
dev.off()

boxplot(
  bd,
  main = "Distance to group centroid",
  ylab = "Distance"
)
dev.off()



## ---------------------------------------------------------
## PCA DIRECTLY ON Z-SCORED hPTM MATRIX
## ---------------------------------------------------------
cat("Running PCA directly on z-scored hPTM matrix...\n")

pca_res <- prcomp(
  ptm_scaled_pca,
  center = FALSE,
  scale. = FALSE
)

pca_var <- summary(pca_res)$importance[2, ] * 100

pca_df <- as.data.frame(pca_res$x[, 1:5, drop = FALSE])
pca_df$sample_id <- rownames(ptm_scaled_pca)
pca_df$diet <- meta_test[pca_df$sample_id, "diet"]




permanova_label <- paste0(
  "PERMANOVA: R² = ",
  sprintf("%.4f", permanova_R2),
  ", F = ",
  sprintf("%.3f", permanova_F),
  ", ",
  format_p_value_plot(permanova_p)
)

## ---------------------------------------------------------
## PCA PLOT WITH PC1/PC2 DENSITIES
## Same layout/style as proteome script
## ---------------------------------------------------------
p_pca <- ggplot(
  pca_df,
  aes(x = PC1, y = PC2, color = diet, fill = diet)
) +
  geom_point(size = 2.7, alpha = 0.85) +
  stat_ellipse(
    geom = "polygon",
    alpha = 0.15,
    color = NA,
    show.legend = FALSE
  ) +
  stat_ellipse(linewidth = 0.8, show.legend = FALSE) +
  scale_color_manual(values = cols2diet) +
  scale_fill_manual(values = cols2diet) +
  labs(
    x = paste0("PC1 (", round(pca_var[1], 1), "%)"),
    y = paste0("PC2 (", round(pca_var[2], 1), "%)"),
    #title = "Midgut hPTM profiles",
    subtitle = permanova_label
  ) +
  theme_cutler() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(
      family = TITLE_FAMILY,
      size = TITLE_SIZE,
      face = "bold"
    )
  )

p_density_x <- ggplot(
  pca_df,
  aes(x = PC1, fill = diet, color = diet)
) +
  geom_density(alpha = 0.25, linewidth = 0.8) +
  scale_color_manual(values = cols2diet) +
  scale_fill_manual(values = cols2diet) +
  theme_void() +
  theme(legend.position = "none")

p_density_y <- ggplot(
  pca_df,
  aes(x = PC2, fill = diet, color = diet)
) +
  geom_density(alpha = 0.25, linewidth = 0.8) +
  scale_color_manual(values = cols2diet) +
  scale_fill_manual(values = cols2diet) +
  coord_flip() +
  theme_void() +
  theme(legend.position = "none")

p_empty <- ggplot() + theme_void()

p_pca_density <-
  (p_density_x + p_empty + plot_layout(widths = c(4, 1))) /
  (p_pca + p_density_y + plot_layout(widths = c(4, 1))) +
  plot_layout(heights = c(1, 4))

ggsave(
  file.path(DIET_OUTPUT_DIR, "hPTM_PCA_with_PC1_PC2_density.pdf"),
  p_pca_density,
  width = 7,
  height = 6
)

ggsave(
  file.path(DIET_OUTPUT_DIR, "hPTM_PCA_with_PC1_PC2_density.tiff"),
  p_pca_density,
  width = 7,
  height = 6,
  dpi = 1200,
  compression = "lzw"
)

## ---------------------------------------------------------
## PC1 / PC2 BOXPLOTS WITH WILCOXON TESTS
## Same layout/style as proteome script
## ---------------------------------------------------------
cat("Testing PC1 and PC2 scores by diet...\n")

pc_box_df <- pca_df %>%
  dplyr::select(sample_id, diet, PC1, PC2) %>%
  pivot_longer(
    cols = c(PC1, PC2),
    names_to = "PC",
    values_to = "score"
  )

pc_tests <- pc_box_df %>%
  group_by(PC) %>%
  summarise(
    p_wilcox = wilcox.test(score ~ diet)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_wilcox, method = "BH"),
    label = vapply(p_adj, p_to_symbol, character(1))
  )

write.csv(
  pc_tests,
  file.path(DIET_OUTPUT_DIR, "hPTM_PC1_PC2_wilcox_tests.csv"),
  row.names = FALSE
)

print(pc_tests)

label_df <- pc_box_df %>%
  group_by(PC) %>%
  summarise(
    y = max(score, na.rm = TRUE) * 1.08,
    .groups = "drop"
  ) %>%
  left_join(pc_tests, by = "PC")

p_pc_box <- ggplot(
  pc_box_df,
  aes(x = diet, y = score, fill = diet)
) +
  geom_violin(trim = FALSE, alpha = 0.5, color = "black") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  geom_jitter(width = 0.08, size = 2, alpha = 0.8) +
  geom_text(
    data = label_df,
    aes(x = 1.5, y = y, label = label),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 6
  ) +
  facet_wrap(~ PC, scales = "free_y") +
  scale_fill_manual(values = cols2diet) +
  labs(
    x = "",
    y = "PC score"
    #title = "Diet differences along hPTM PCA axes"
  ) +
  theme_cutler() +
  theme(legend.position = "none")

ggsave(
  file.path(DIET_OUTPUT_DIR, "hPTM_PC1_PC2_boxplots.pdf"),
  p_pc_box,
  width = 6,
  height = 4
)

ggsave(
  file.path(DIET_OUTPUT_DIR, "hPTM_PC1_PC2_boxplots.tiff"),
  p_pc_box,
  width = 6,
  height = 4,
  dpi = 1200,
  compression = "lzw"
)

## ---------------------------------------------------------
## COMBINED SUMMARY TABLE
## ---------------------------------------------------------
cat("Saving combined summary tables...\n")

bd_perm_tab <- as.data.frame(bd_perm$tab)
bd_perm_tab$term <- rownames(bd_perm_tab)
rownames(bd_perm_tab) <- NULL

bd_group_row <- bd_perm_tab[
  bd_perm_tab$term == "Groups",
  ,
  drop = FALSE
]

betadisper_summary <- data.frame(
  test = "betadisper_permutation",
  omics = "hPTM",
  factor = "diet",
  F_value = bd_group_row$F,
  p_value = bd_group_row$`Pr(>F)`,
  significant = bd_group_row$`Pr(>F)` < 0.05
)

write.csv(
  betadisper_summary,
  file.path(DIET_OUTPUT_DIR, "hPTM_betadisper_summary.csv"),
  row.names = FALSE
)



## ---------------------------------------------------------
## SAVE SESSION INFO
## ---------------------------------------------------------
capture.output(
  sessionInfo(),
  file = file.path(DIET_OUTPUT_DIR, "sessionInfo.txt")
)

cat("\nDone.\n")
cat("Outputs saved in:", OUTPUT_DIR, "\n")

############################################################
## ADDITIONAL METADATA GROUPINGS
############################################################

ADDITIONAL_GROUPS <- list(
  list(column = "phase",       display = "Phase",                   folder = "by_phase"),
  list(column = "diet_phase",  display = "Diet x Temporal Phase",  folder = "by_diet_phase"),
  list(column = "day_categ",   display = "Day",                     folder = "by_day_categ")
)

run_additional_grouping <- function(group_column, group_display, folder_name) {
  cat("\n=========================================================\n")
  cat("Analysing grouping variable:", group_column,
      " (displayed as: ", group_display, ")\n", sep = "")

  if (!group_column %in% names(meta_df)) {
    warning("Metadata column not found; skipping: ", group_column)
    return(invisible(NULL))
  }

  group_values <- trimws(as.character(meta_df[[group_column]]))
  group_values[group_values %in% c("", "NA", "NaN")] <- NA_character_
  keep_group <- !is.na(group_values)

  if (sum(keep_group) < 3) {
    warning("Fewer than 3 samples have non-missing values for ", group_column,
            "; skipping this grouping.")
    return(invisible(NULL))
  }

  ## Keep the order in which levels occur in the metadata.
  group_levels <- unique(group_values[keep_group])

  if (length(group_levels) < 2) {
    warning("Fewer than 2 groups are present for ", group_column,
            "; skipping this grouping.")
    return(invisible(NULL))
  }

  group_output_dir <- file.path(OUTPUT_DIR, folder_name)
  dir.create(group_output_dir, showWarnings = FALSE, recursive = TRUE)

  ## Same matrix and same sample order as the original analysis.
  meta_group <- data.frame(
    sample_id = rownames(meta_df)[keep_group],
    group = factor(group_values[keep_group], levels = group_levels),
    stringsAsFactors = FALSE
  )
  rownames(meta_group) <- meta_group$sample_id

  stopifnot(all(meta_group$sample_id %in% rownames(ptm_scaled_pca)))

  ## Subsetting the original Euclidean distance matrix is equivalent to
  ## recalculating Euclidean distances on the same already-imputed z-scores.
  dist_group <- as.dist(
    as.matrix(dist_mat)[meta_group$sample_id, meta_group$sample_id, drop = FALSE]
  )

  pca_group_df <- pca_df[
    match(meta_group$sample_id, pca_df$sample_id),
    ,
    drop = FALSE
  ]
  pca_group_df$group <- meta_group[pca_group_df$sample_id, "group"]

  stopifnot(identical(pca_group_df$sample_id, meta_group$sample_id))



  ## -------------------------------------------------------
  ## PERMANOVA
  ## Exact same extraction logic as the original DIET block:
  ## take the first non-Residual/non-Total row.
  ## -------------------------------------------------------
  cat("Running PERMANOVA...\n")

  set.seed(SEED)
  permanova_res_group <- adonis2(
    dist_group ~ group,
    data = meta_group,
    permutations = N_PERM
  )

  permanova_tab_group <- as.data.frame(permanova_res_group)
  permanova_tab_group$term <- rownames(permanova_tab_group)
  rownames(permanova_tab_group) <- NULL

  tested_row_group <- permanova_tab_group[
    !(permanova_tab_group$term %in% c("Residual", "Total")),
    ,
    drop = FALSE
  ]

  if (nrow(tested_row_group) < 1) {
    stop(
      "Could not extract the PERMANOVA term for: ", group_column,
      ". Returned rows: ",
      paste(permanova_tab_group$term, collapse = ", ")
    )
  }

  tested_row_group <- tested_row_group[1, , drop = FALSE]
  tested_row_group$term <- group_column

  permanova_R2_group <- tested_row_group$R2
  permanova_F_group  <- tested_row_group$F
  permanova_p_group  <- tested_row_group$`Pr(>F)`



  permanova_summary_group <- data.frame(
    test = "PERMANOVA",
    omics = "hPTM",
    factor = group_column,
    factor_display = group_display,
    df = tested_row_group$Df,
    sum_of_squares = tested_row_group$SumOfSqs,
    R2 = permanova_R2_group,
    F_value = permanova_F_group,
    p_value = permanova_p_group,
    significant = permanova_p_group < 0.05
  )

  write.csv(
    permanova_summary_group,
    file.path(group_output_dir, "hPTM_PERMANOVA_summary.csv"),
    row.names = FALSE
  )

  print(permanova_res_group)

  ## -------------------------------------------------------
  ## BETADISPER
  ## Same procedure as the original DIET block.
  ## -------------------------------------------------------
  cat("Checking homogeneity of multivariate dispersion with betadisper...\n")

  bd_group <- betadisper(dist_group, meta_group$group)
  bd_anova_group <- anova(bd_group)

  set.seed(SEED)
  bd_perm_group <- permutest(bd_group, permutations = N_PERM)



  print(bd_anova_group)
  print(bd_perm_group)

  plot(bd_group, main = paste0("hPTM betadisper: ", group_display))
  dev.off()

  boxplot(
    bd_group,
    main = paste0("Distance to group centroid: ", group_display),
    ylab = "Distance"
  )
  dev.off()



  bd_perm_tab_group <- as.data.frame(bd_perm_group$tab)
  bd_perm_tab_group$term <- rownames(bd_perm_tab_group)
  rownames(bd_perm_tab_group) <- NULL

  bd_group_row <- bd_perm_tab_group[
    bd_perm_tab_group$term == "Groups",
    ,
    drop = FALSE
  ]

  betadisper_summary_group <- data.frame(
    test = "betadisper_permutation",
    omics = "hPTM",
    factor = group_column,
    factor_display = group_display,
    F_value = bd_group_row$F,
    p_value = bd_group_row$`Pr(>F)`,
    significant = bd_group_row$`Pr(>F)` < 0.05
  )

  write.csv(
    betadisper_summary_group,
    file.path(group_output_dir, "hPTM_betadisper_summary.csv"),
    row.names = FALSE
  )

  ## -------------------------------------------------------
  ## PCA PLOT WITH PC1/PC2 DENSITIES
  ## PCA itself is the original PCA above; this only changes
  ## the grouping used for colour, fill, legend and annotation.
  ## -------------------------------------------------------
  group_cols <- stats::setNames(
    grDevices::hcl.colors(nlevels(meta_group$group), palette = "Dark 3"),
    levels(meta_group$group)
  )

  permanova_label_group <- paste0(
    "PERMANOVA: R² = ",
    sprintf("%.4f", permanova_R2_group),
    ", F = ",
    sprintf("%.3f", permanova_F_group),
    ", ",
    format_p_value_plot(permanova_p_group)
  )

  p_pca_group <- ggplot(
    pca_group_df,
    aes(x = PC1, y = PC2, color = group, fill = group)
  ) +
    geom_point(size = 2.7, alpha = 0.85) +
    stat_ellipse(
      geom = "polygon",
      alpha = 0.15,
      color = NA,
      show.legend = FALSE
    ) +
    stat_ellipse(linewidth = 0.8, show.legend = FALSE) +
    scale_color_manual(values = group_cols, name = group_display) +
    scale_fill_manual(values = group_cols, name = group_display) +
    labs(
      x = paste0("PC1 (", round(pca_var[1], 1), "%)"),
      y = paste0("PC2 (", round(pca_var[2], 1), "%)"),
      subtitle = permanova_label_group
    ) +
    theme_cutler() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(
        family = TITLE_FAMILY,
        size = TITLE_SIZE,
        face = "bold"
      )
    )

  p_density_x_group <- ggplot(
    pca_group_df,
    aes(x = PC1, fill = group, color = group)
  ) +
    geom_density(alpha = 0.25, linewidth = 0.8) +
    scale_color_manual(values = group_cols) +
    scale_fill_manual(values = group_cols) +
    theme_void() +
    theme(legend.position = "none")

  p_density_y_group <- ggplot(
    pca_group_df,
    aes(x = PC2, fill = group, color = group)
  ) +
    geom_density(alpha = 0.25, linewidth = 0.8) +
    scale_color_manual(values = group_cols) +
    scale_fill_manual(values = group_cols) +
    coord_flip() +
    theme_void() +
    theme(legend.position = "none")

  p_empty_group <- ggplot() + theme_void()

  p_pca_density_group <-
    (p_density_x_group + p_empty_group + plot_layout(widths = c(4, 1))) /
    (p_pca_group + p_density_y_group + plot_layout(widths = c(4, 1))) +
    plot_layout(heights = c(1, 4))

  ggsave(
    file.path(group_output_dir, "hPTM_PCA_with_PC1_PC2_density.pdf"),
    p_pca_density_group,
    width = 7,
    height = 6
  )

  ggsave(
    file.path(group_output_dir, "hPTM_PCA_with_PC1_PC2_density.tiff"),
    p_pca_density_group,
    width = 7,
    height = 6,
    dpi = 1200,
    compression = "lzw"
  )

  ## -------------------------------------------------------
  ## PC1 / PC2 BOXPLOTS
  ## For 2 groups, this is the original Wilcoxon test.
  ## For 3+ groups, one global Kruskal-Wallis test is required.
  ## -------------------------------------------------------
  cat("Testing PC1 and PC2 scores by ", group_column, "...\n", sep = "")

  pc_box_df_group <- pca_group_df %>%
    dplyr::select(sample_id, group, PC1, PC2) %>%
    pivot_longer(
      cols = c(PC1, PC2),
      names_to = "PC",
      values_to = "score"
    )

  pc_tests_group <- pc_box_df_group %>%
    group_by(PC) %>%
    group_modify(~ {
      n_groups <- nlevels(droplevels(.x$group))

      if (n_groups == 2) {
        data.frame(
          test = "Wilcoxon rank-sum",
          p_raw = wilcox.test(score ~ group, data = .x)$p.value
        )
      } else {
        data.frame(
          test = "Kruskal-Wallis",
          p_raw = kruskal.test(score ~ group, data = .x)$p.value
        )
      }
    }) %>%
    ungroup() %>%
    mutate(
      p_adj = p.adjust(p_raw, method = "BH"),
      label = vapply(p_adj, p_to_symbol, character(1))
    )

  write.csv(
    pc_tests_group,
    file.path(group_output_dir, "hPTM_PC1_PC2_group_tests.csv"),
    row.names = FALSE
  )

  print(pc_tests_group)

  label_df_group <- pc_box_df_group %>%
    group_by(PC) %>%
    summarise(
      score_min = min(score, na.rm = TRUE),
      score_max = max(score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      y = score_max + pmax(0.08 * (score_max - score_min), 0.05)
    ) %>%
    left_join(pc_tests_group, by = "PC")

  p_pc_box_group <- ggplot(
    pc_box_df_group,
    aes(x = group, y = score, fill = group)
  ) +
    geom_violin(trim = FALSE, alpha = 0.5, color = "black") +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
    geom_jitter(width = 0.08, size = 2, alpha = 0.8) +
    geom_text(
      data = label_df_group,
      aes(x = (nlevels(meta_group$group) + 1) / 2, y = y, label = label),
      inherit.aes = FALSE,
      fontface = "bold",
      size = 6
    ) +
    facet_wrap(~ PC, scales = "free_y") +
    scale_fill_manual(values = group_cols) +
    labs(
      x = group_display,
      y = "PC score"
      #title = paste0(group_display, " differences along hPTM PCA axes")
    ) +
    theme_cutler() +
    theme(legend.position = "none")

  ggsave(
    file.path(group_output_dir, "hPTM_PC1_PC2_boxplots.pdf"),
    p_pc_box_group,
    width = 6,
    height = 4
  )

  ggsave(
    file.path(group_output_dir, "hPTM_PC1_PC2_boxplots.tiff"),
    p_pc_box_group,
    width = 6,
    height = 4,
    dpi = 1200,
    compression = "lzw"
  )



  capture.output(
    sessionInfo(),
    file = file.path(group_output_dir, "sessionInfo.txt")
  )

  cat("Additional outputs saved in:", group_output_dir, "\n")

  invisible(
    list(
      permanova_summary = permanova_summary_group,
      betadisper_summary = betadisper_summary_group,
      pc_tests = pc_tests_group
    )
  )
}

## Run only the extra metadata analyses. The original diet section above is not changed.
additional_results <- lapply(
  ADDITIONAL_GROUPS,
  function(x) {
    run_additional_grouping(
      group_column = x$column,
      group_display = x$display,
      folder_name = x$folder
    )
  }
)
names(additional_results) <- vapply(ADDITIONAL_GROUPS, `[[`, character(1), "column")


cat("\nAll additional metadata analyses completed.\n")
