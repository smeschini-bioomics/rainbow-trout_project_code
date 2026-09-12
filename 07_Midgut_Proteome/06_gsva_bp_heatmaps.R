
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(limma)
  library(sva)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(grid)
  library(circlize)
})

# Settings
TISSUE <- "Liver"   # change to "Midgut" in the midgut module

INPUT_RDS <- "05_gsva/sample_geneterms_matrix_fullnames.rds"
OUTPUT_DIR <- "06_gsva_bp_heatmaps"
STAT_DIR <- file.path(OUTPUT_DIR, "statistics")

SEED <- 123
FDR_MAIN <- 0.05
FDR_SUPPLEMENTARY <- if (TISSUE == "Liver") 1e-4 else 0.05
TOP_N_UP <- 3
TOP_N_DOWN <- 3

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(STAT_DIR, recursive = TRUE, showWarnings = FALSE)

# Load GSVA scores and metadata
gsva_obj <- readRDS(INPUT_RDS)

if (!all(c("es_bp", "metadata") %in% names(gsva_obj))) {
  stop("The GSVA RDS must contain objects named 'es_bp' and 'metadata'.")
}

es_bp <- gsva_obj$es_bp
meta_filt <- gsva_obj$metadata

required_meta <- c("sample", "diet", "day_categ")
missing_meta <- setdiff(required_meta, names(meta_filt))

if (length(missing_meta) > 0) {
  stop(
    "Missing metadata column(s): ",
    paste(missing_meta, collapse = ", ")
  )
}

# Prepare and align metadata
meta_filt <- meta_filt %>%
  mutate(
    diet = factor(diet, levels = c("NC", "HC")),
    day_categ = factor(
      day_categ,
      levels = c(
        "day_01", "day_02", "day_03", "day_04",
        "day_10", "day_15", "day_22"
      )
    )
  )

keep_samples <- intersect(colnames(es_bp), meta_filt$sample)

es_bp_filt <- es_bp[, keep_samples, drop = FALSE]

meta_filt <- meta_filt %>%
  filter(sample %in% keep_samples) %>%
  arrange(match(sample, keep_samples)) %>%
  mutate(
    diet = droplevels(diet),
    day_categ = droplevels(day_categ)
  )

stopifnot(identical(colnames(es_bp_filt), meta_filt$sample))

# SVA: preserve diet and day main effects
mod_sva <- model.matrix(
  ~ diet + day_categ,
  data = meta_filt
)

mod0_sva <- model.matrix(
  ~ 1,
  data = meta_filt
)

if (qr(mod_sva)$rank != ncol(mod_sva)) {
  stop("The SVA design matrix is not full rank.")
}

set.seed(SEED)

sv <- sva(
  es_bp_filt,
  mod_sva,
  mod0_sva,
  vfilter = min(2500, nrow(es_bp_filt))
)

# Full limma interaction model with surrogate variables
mod_full_clean <- model.matrix(
  ~ diet * day_categ,
  data = meta_filt
)

colnames(mod_full_clean) <- make.names(
  colnames(mod_full_clean)
)

sv_mat <- sv$sv

if (is.null(sv_mat)) {
  sv_mat <- matrix(
    numeric(0),
    nrow = nrow(meta_filt),
    ncol = 0
  )
} else {
  colnames(sv_mat) <- paste0(
    "SV",
    seq_len(ncol(sv_mat))
  )
}

mod_sv <- cbind(
  mod_full_clean,
  sv_mat
)

fit <- lmFit(
  es_bp_filt,
  mod_sv
)

coef_names <- colnames(coef(fit))

# HC versus NC contrast at each sampling day
day_levels <- levels(meta_filt$day_categ)

contrast_strings <- vapply(
  day_levels,
  function(day_level) {
    if (day_level == day_levels[1]) {
      return("dietHC")
    }

    interaction_name <- paste0(
      "dietHC.day_categ",
      day_level
    )

    if (!interaction_name %in% coef_names) {
      stop(
        "Interaction coefficient not found: ",
        interaction_name
      )
    }

    paste0(
      "dietHC + ",
      interaction_name
    )
  },
  character(1)
)

names(contrast_strings) <- paste0(
  "HC_vs_NC_",
  day_levels
)

cont_mat <- makeContrasts(
  contrasts = contrast_strings,
  levels = coef_names
)

# Standard limma sequence: fit -> contrasts -> empirical Bayes
fit_contr <- contrasts.fit(
  fit,
  cont_mat
)

fit_contr <- eBayes(
  fit_contr,
  robust = TRUE
)

contrast_names <- colnames(cont_mat)

# Complete pathway-level results
all_results <- bind_rows(
  lapply(
    contrast_names,
    function(contrast_name) {
      topTable(
        fit_contr,
        coef = contrast_name,
        number = Inf,
        sort.by = "none"
      ) %>%
        rownames_to_column("pathway") %>%
        mutate(
          tissue = TISSUE,
          contrast = contrast_name,
          day = as.integer(
            sub(
              "^HC_vs_NC_day_0?",
              "",
              contrast_name
            )
          ),
          direction = case_when(
            logFC > 0 ~ "Higher in HC",
            logFC < 0 ~ "Lower in HC",
            TRUE ~ "No difference"
          )
        ) %>%
        select(
          tissue,
          contrast,
          day,
          pathway,
          logFC,
          AveExpr,
          t,
          P.Value,
          adj.P.Val,
          B,
          direction
        )
    }
  )
)

# Sample sizes for each day-specific contrast
sample_counts <- meta_filt %>%
  mutate(
    day = as.integer(
      sub(
        "^day_0?",
        "",
        as.character(day_categ)
      )
    )
  ) %>%
  count(
    day,
    diet,
    name = "n"
  ) %>%
  complete(
    day = c(1, 2, 3, 4, 10, 15, 22),
    diet = c("NC", "HC"),
    fill = list(n = 0)
  ) %>%
  pivot_wider(
    names_from = diet,
    values_from = n
  ) %>%
  mutate(
    n_total = HC + NC
  )

# Numbers tested and significant for each day
summary_by_day <- all_results %>%
  group_by(
    tissue,
    contrast,
    day
  ) %>%
  summarise(
    n_GO_BP_pathways_tested = n(),
    n_significant_FDR_main = sum(
      adj.P.Val < FDR_MAIN,
      na.rm = TRUE
    ),
    n_higher_HC_FDR_main = sum(
      adj.P.Val < FDR_MAIN & logFC > 0,
      na.rm = TRUE
    ),
    n_lower_HC_FDR_main = sum(
      adj.P.Val < FDR_MAIN & logFC < 0,
      na.rm = TRUE
    ),
    n_significant_FDR_supplementary = sum(
      adj.P.Val < FDR_SUPPLEMENTARY,
      na.rm = TRUE
    ),
    n_higher_HC_FDR_supplementary = sum(
      adj.P.Val < FDR_SUPPLEMENTARY & logFC > 0,
      na.rm = TRUE
    ),
    n_lower_HC_FDR_supplementary = sum(
      adj.P.Val < FDR_SUPPLEMENTARY & logFC < 0,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  left_join(
    sample_counts,
    by = "day"
  ) %>%
  mutate(
    FDR_main = FDR_MAIN,
    FDR_supplementary = FDR_SUPPLEMENTARY
  ) %>%
  select(
    tissue,
    contrast,
    day,
    HC,
    NC,
    n_total,
    n_GO_BP_pathways_tested,
    FDR_main,
    n_significant_FDR_main,
    n_higher_HC_FDR_main,
    n_lower_HC_FDR_main,
    FDR_supplementary,
    n_significant_FDR_supplementary,
    n_higher_HC_FDR_supplementary,
    n_lower_HC_FDR_supplementary
  ) %>%
  arrange(day)

# Top pathways for the main contrast heatmap
selected_per_day_df <- all_results %>%
  filter(
    adj.P.Val < FDR_MAIN,
    logFC != 0
  ) %>%
  group_by(
    contrast,
    day,
    direction
  ) %>%
  arrange(
    adj.P.Val,
    desc(abs(logFC)),
    .by_group = TRUE
  ) %>%
  mutate(rank_in_group = row_number()) %>%
  filter(
    rank_in_group <= if_else(direction == "Higher in HC", TOP_N_UP, TOP_N_DOWN)
  ) %>%
  select(-rank_in_group) %>%
  ungroup()

DE_bp <- unique(selected_per_day_df$pathway)

if (length(DE_bp) == 0) {
  stop("No pathways met the main-figure threshold.")
}

main_figure_summary <- data.frame(
  tissue = TISSUE,
  FDR = FDR_MAIN,
  n_unique_pathways_displayed = length(DE_bp)
)

# Export statistics from this authoritative run
write.csv(
  summary_by_day,
  file.path(
    STAT_DIR,
    "GSVA_sample_size_tested_significant_by_day.csv"
  ),
  row.names = FALSE
)

write.csv(
  all_results,
  file.path(
    STAT_DIR,
    "GSVA_all_pathways_all_day_contrasts.csv"
  ),
  row.names = FALSE
)

write.csv(
  all_results %>%
    filter(adj.P.Val < FDR_MAIN) %>%
    arrange(
      day,
      adj.P.Val,
      desc(abs(logFC))
    ),
  file.path(
    STAT_DIR,
    "GSVA_significant_pathways_FDR_main.csv"
  ),
  row.names = FALSE
)

write.csv(
  all_results %>%
    filter(adj.P.Val < FDR_SUPPLEMENTARY) %>%
    arrange(
      day,
      adj.P.Val,
      desc(abs(logFC))
    ),
  file.path(
    STAT_DIR,
    "GSVA_significant_pathways_FDR_supplementary.csv"
  ),
  row.names = FALSE
)

write.csv(
  selected_per_day_df,
  file.path(
    STAT_DIR,
    "GSVA_main_figure_top3_up_down.csv"
  ),
  row.names = FALSE
)

write.csv(
  main_figure_summary,
  file.path(
    STAT_DIR,
    "GSVA_main_figure_unique_pathway_count.csv"
  ),
  row.names = FALSE
)

saveRDS(
  list(
    tissue = TISSUE,
    seed = SEED,
    FDR_main = FDR_MAIN,
    FDR_supplementary = FDR_SUPPLEMENTARY,
    es_bp_filt = es_bp_filt,
    metadata = meta_filt,
    mod_sva = mod_sva,
    mod0_sva = mod0_sva,
    sv = sv,
    sv_mat = sv_mat,
    mod_full_clean = mod_full_clean,
    mod_sv = mod_sv,
    cont_mat = cont_mat,
    fit_contr = fit_contr,
    all_results = all_results,
    sample_counts = sample_counts,
    summary_by_day = summary_by_day,
    selected_per_day_df = selected_per_day_df
  ),
  file.path(
    STAT_DIR,
    "GSVA_SVA_limma_complete_results.rds"
  )
)

# Main heatmap matrices
adjP_mat <- sapply(
  contrast_names,
  function(contrast_name) {
    tt <- topTable(
      fit_contr,
      coef = contrast_name,
      number = Inf,
      sort.by = "none"
    )

    tt$adj.P.Val
  }
)

rownames(adjP_mat) <- rownames(fit_contr)
colnames(adjP_mat) <- contrast_names

logFC_mat <- coef(fit_contr)[
  DE_bp,
  contrast_names,
  drop = FALSE
]

adjP_DE <- adjP_mat[
  DE_bp,
  contrast_names,
  drop = FALSE
]

mat_to_plot <- logFC_mat

max_abs_main <- max(
  abs(mat_to_plot),
  na.rm = TRUE
)

pal <- brewer.pal(
  11,
  "RdBu"
)

col_fun_main <- colorRamp2(
  c(
    -max_abs_main,
    0,
    max_abs_main
  ),
  pal[c(11, 6, 1)]
)

pval_symbols <- matrix(
  "",
  nrow = nrow(adjP_DE),
  ncol = ncol(adjP_DE),
  dimnames = dimnames(adjP_DE)
)

pval_symbols[adjP_DE < 0.05] <- "*"
pval_symbols[adjP_DE < 0.01] <- "**"
pval_symbols[adjP_DE < 0.001] <- "***"
pval_symbols[adjP_DE < 1e-4] <- "****"

new_colnames <- paste0(
  "Day ",
  sub(
    "^HC_vs_NC_day_0?",
    "",
    colnames(mat_to_plot)
  )
)

colnames(mat_to_plot) <- new_colnames
colnames(pval_symbols) <- new_colnames

heatmap_top_bp <- Heatmap(
  mat_to_plot,
  name = "GSVA score (HC - NC)",
  col = col_fun_main,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_column_names = TRUE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 24),
  column_names_gp = gpar(fontsize = 20),
  column_names_rot = 90,
  row_names_max_width = unit(15, "cm"),
  width = unit(10, "cm"),
  heatmap_legend_param = list(
    title = "GSVA score (HC - NC)",
    direction = "horizontal",
    title_position = "topcenter",
    title_gp = gpar(
      fontsize = 20,
      fontface = "bold"
    ),
    labels_gp = gpar(fontsize = 18),
    legend_width = unit(10, "cm"),
    grid_width = unit(0.9, "cm"),
    grid_height = unit(0.9, "cm")
  ),
  cell_fun = function(j, i, x, y, w, h, fill) {
    symbol <- pval_symbols[i, j]

    if (nzchar(symbol)) {
      grid.text(
        symbol,
        x = x,
        y = y,
        gp = gpar(
          fontsize = 30,
          fontface = "bold",
          col = "black"
        )
      )
    }
  }
)

sig_legend <- Legend(
  title = "Significance",
  labels = c(
    "*   FDR < 0.05",
    "**  FDR < 0.01",
    "*** FDR < 0.001",
    "**** FDR < 1e-4"
  ),
  type = "text",
  title_gp = gpar(
    fontsize = 20,
    fontface = "bold"
  ),
  labels_gp = gpar(fontsize = 18)
)

pdf(
  file.path(
    OUTPUT_DIR,
    "heatmap_top_bp.pdf"
  ),
  width = 13.5,
  height = 12,
  useDingbats = FALSE
)

draw(
  heatmap_top_bp,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "right",
  heatmap_legend_list = list(sig_legend),
  padding = unit(
    c(2, 2, 2, 40),
    "mm"
  )
)

dev.off()

tiff(
  file.path(
    OUTPUT_DIR,
    "heatmap_top_bp.tiff"
  ),
  width = 13.5,
  height = 12,
  units = "in",
  res = 1200,
  compression = "lzw"
)

draw(
  heatmap_top_bp,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "right",
  heatmap_legend_list = list(sig_legend),
  padding = unit(
    c(2, 2, 2, 40),
    "mm"
  )
)

dev.off()

# Supplementary sample-level heatmap
DE_union <- unique(
  all_results$pathway[
    all_results$adj.P.Val < FDR_SUPPLEMENTARY
  ]
)

if (length(DE_union) == 0) {
  stop(
    "No pathways met the supplementary threshold: FDR < ",
    FDR_SUPPLEMENTARY
  )
}

logFC_all_DE <- coef(fit_contr)[
  DE_union,
  contrast_names,
  drop = FALSE
]

max_abs_effect <- apply(
  abs(logFC_all_DE),
  1,
  max,
  na.rm = TRUE
)

DE_union <- DE_union[
  order(
    max_abs_effect,
    decreasing = TRUE
  )
]

DE_es <- es_bp_filt[
  DE_union,
  ,
  drop = FALSE
]

if (ncol(sv_mat) > 0) {
  DE_es_clean <- removeBatchEffect(
    DE_es,
    covariates = sv_mat,
    design = mod_full_clean
  )
} else {
  DE_es_clean <- DE_es
}

do_zscore <- FALSE

if (do_zscore) {
  mat_plot <- t(
    scale(
      t(DE_es_clean)
    )
  )
  mat_plot[is.na(mat_plot)] <- 0
  hm_name <- "GSVA (z)"
} else {
  mat_plot <- DE_es_clean
  hm_name <- "GSVA score"
}

max_abs_supp <- max(
  abs(mat_plot),
  na.rm = TRUE
)

col_fun_supp <- colorRamp2(
  c(
    -max_abs_supp,
    0,
    max_abs_supp
  ),
  rev(
    colorRampPalette(
      brewer.pal(11, "RdBu")
    )(256)
  )[c(1, 128, 256)]
)

legend_breaks <- round(
  seq(
    -max_abs_supp,
    max_abs_supp,
    length.out = 5
  ),
  2
)

sample_hc <- hclust(
  as.dist(
    1 - cor(
      mat_plot,
      method = "spearman",
      use = "pairwise.complete.obs"
    )
  ),
  method = "complete"
)

path_hc <- hclust(
  as.dist(
    1 - cor(
      t(mat_plot),
      method = "pearson",
      use = "pairwise.complete.obs"
    )
  ),
  method = "complete"
)

day_num <- as.integer(
  sub(
    "^day_0?",
    "",
    as.character(meta_filt$day_categ)
  )
)

meta_filt <- meta_filt %>%
  mutate(
    phase = factor(
      if_else(
        day_num %in% c(1, 2, 3, 4),
        "Early",
        "Late"
      ),
      levels = c("Early", "Late")
    ),
    diet_phase = factor(
      paste(
        diet,
        phase,
        sep = " - "
      ),
      levels = c(
        "NC - Early",
        "NC - Late",
        "HC - Early",
        "HC - Late"
      )
    )
  )

day_lab <- factor(
  paste0(
    "Day ",
    day_num
  ),
  levels = c(
    "Day 1",
    "Day 2",
    "Day 3",
    "Day 4",
    "Day 10",
    "Day 15",
    "Day 22"
  )
)

diet_cols <- c(
  "NC" = "#44a7c4",
  "HC" = "#f99943"
)

day_cols <- c(
  "Day 1" = "#1b9e77",
  "Day 2" = "#d95f02",
  "Day 3" = "#7570b3",
  "Day 4" = "#e7298a",
  "Day 10" = "#66a61e",
  "Day 15" = "#e6ab02",
  "Day 22" = "#a6761d"
)

phase_cols <- c(
  "Early" = "#7b6fd6",
  "Late" = "#4daf4a"
)

dtg_cols <- c(
  "NC - Early" = "#6baed6",
  "NC - Late" = "#2171b5",
  "HC - Early" = "#fdae6b",
  "HC - Late" = "#d94801"
)

ha_all <- HeatmapAnnotation(
  diet = meta_filt$diet,
  day = day_lab,
  phase = meta_filt$phase,
  diet_phase = meta_filt$diet_phase,
  col = list(
    diet = diet_cols,
    day = day_cols,
    phase = phase_cols,
    diet_phase = dtg_cols
  ),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    diet = list(
      title = "Diet",
      title_gp = gpar(
        fontsize = 18,
        fontface = "bold"
      ),
      labels_gp = gpar(fontsize = 16)
    ),
    day = list(
      title = "Day",
      title_gp = gpar(
        fontsize = 18,
        fontface = "bold"
      ),
      labels_gp = gpar(fontsize = 16),
      ncol = 2
    ),
    phase = list(
      title = "Temporal Phase",
      title_gp = gpar(
        fontsize = 18,
        fontface = "bold"
      ),
      labels_gp = gpar(fontsize = 16)
    ),
    diet_phase = list(
      title = "Diet x Temporal Phase",
      title_gp = gpar(
        fontsize = 18,
        fontface = "bold"
      ),
      labels_gp = gpar(fontsize = 16),
      ncol = 2
    )
  )
)

ht_bp_samples <- Heatmap(
  mat_plot,
  name = hm_name,
  col = col_fun_supp,
  top_annotation = ha_all,
  height = unit(14, "cm"),
  use_raster = TRUE,
  raster_device = "png",
  raster_quality = 5,
  cluster_columns = as.dendrogram(sample_hc),
  cluster_rows = as.dendrogram(path_hc),
  show_column_names = FALSE,
  show_row_names = FALSE,
  column_title = NULL,
  row_title = paste0(
    "DE pathways FDR < ",
    format(
      FDR_SUPPLEMENTARY,
      scientific = FALSE
    )
  ),
  row_title_gp = gpar(fontsize = 18),
  heatmap_legend_param = list(
    title = "GSVA score",
    at = legend_breaks,
    labels = legend_breaks,
    title_gp = gpar(
      fontsize = 20,
      fontface = "bold"
    ),
    labels_gp = gpar(fontsize = 16),
    direction = "vertical",
    legend_height = unit(14, "cm"),
    title_position = "leftcenter-rot"
  )
)

pdf(
  file.path(
    OUTPUT_DIR,
    "BP_sample_pathway_heatmap.pdf"
  ),
  width = 12,
  height = 8,
  useDingbats = FALSE
)

draw(
  ht_bp_samples,
  heatmap_legend_side = "right",
  annotation_legend_side = "bottom",
  align_heatmap_legend = "heatmap_center",
  merge_legends = FALSE
)

dev.off()

tiff(
  file.path(
    OUTPUT_DIR,
    "BP_sample_pathway_heatmap.tiff"
  ),
  width = 10,
  height = 10,
  units = "in",
  res = 1200,
  compression = "lzw"
)

draw(
  ht_bp_samples,
  heatmap_legend_side = "right",
  annotation_legend_side = "bottom",
  align_heatmap_legend = "heatmap_center",
  merge_legends = FALSE
)

dev.off()

writeLines(
  capture.output(sessionInfo()),
  file.path(
    OUTPUT_DIR,
    "sessionInfo.txt"
  )
)

print(summary_by_day)
print(main_figure_summary)

message("Analysis and exports completed: ", OUTPUT_DIR)
