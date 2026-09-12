suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# ---------- CONFIG ------------------------------------------------------------

infile <- "export_for_R/spearman_corr_digesta/digesta_spearman_CLR_results_LONG.tsv"

physio_cols <- c(
  "Body_Weight",
  "Plasma_Glucose",
  "Plasma_L-Lactate",
  "Plasma_D-Lactate",
  "HSI"
)

physio_labels <- c(
  "Body Weight",
  "Plasma Glucose",
  "Plasma L-Lactate",
  "Plasma D-Lactate",
  "HSI"
)
names(physio_labels) <- physio_cols

genus_pattern <- "g__"
fdr_keep_asv <- 0.001

plot_title <- "Spearman correlations"

# Output directory: same directory as input file
out_dir <- dirname(infile)

pdf_file <- file.path(out_dir, "digesta_spearman_CLR_heatmap.pdf")
tiff_file <- file.path(out_dir, "digesta_spearman_CLR_heatmap.tiff")

# ---------- LOAD PYTHON RESULTS -----------------------------------------------

df <- read_tsv(infile, show_col_types = FALSE)

required_cols <- c(
  "asv_id",
  "asv_label",
  "tax_label",
  "physiological_trait",
  "spearman_rho",
  "BH_adjusted_p",
  "significance"
)

missing_cols <- setdiff(required_cols, names(df))

if (length(missing_cols) > 0) {
  stop(
    "Missing columns in input TSV: ",
    paste(missing_cols, collapse = ", ")
  )
}

# Python exported physiological-trait display labels.
python_to_original_physio <- c(
  "Body Weight" = "Body_Weight",
  "Plasma Glucose" = "Plasma_Glucose",
  "Plasma L-Lactate" = "Plasma_L-Lactate",
  "Plasma D-Lactate" = "Plasma_D-Lactate",
  "Hepato-Somatic Index (HSI)" = "HSI"
)

corr_df <- df %>%
  mutate(
    physio = unname(python_to_original_physio[physiological_trait]),
    rho = as.numeric(spearman_rho),
    p_adj = as.numeric(BH_adjusted_p),
    stars = replace_na(significance, "")
  ) %>%
  filter(!is.na(physio))

if (nrow(corr_df) == 0) {
  stop("No valid physiological-trait names were recognised in the input file.")
}

# ---------- FILTER GENUS-LEVEL ASVs AND FDR ----------------------------------

corr_filt <- corr_df %>%
  filter(str_detect(tax_label, fixed(genus_pattern))) %>%
  group_by(asv_id, asv_label, tax_label) %>%
  filter(any(p_adj < fdr_keep_asv, na.rm = TRUE)) %>%
  ungroup()

if (nrow(corr_filt) == 0) {
  stop(
    "No ASVs remained after genus-level and FDR filtering ",
    "(FDR < ", fdr_keep_asv, ")."
  )
}

# Shorten only the very long genus name for heatmap display.
# The original taxonomy remains unchanged in the input/output tables.
corr_filt <- corr_filt %>%
  mutate(
    tax_label_plot = str_replace(
      tax_label,
      "g__Burkholderia-Caballeronia-Paraburkholderia",
      "g__B-C-P"
    ),
    plot_label = if_else(
      duplicated(tax_label_plot) | duplicated(tax_label_plot, fromLast = TRUE),
      paste0(asv_id, " | ", tax_label_plot),
      tax_label_plot
    )
  )

# ---------- BUILD HEATMAP MATRICES --------------------------------------------

mat_rho <- corr_filt %>%
  select(plot_label, physio, rho) %>%
  pivot_wider(
    names_from = physio,
    values_from = rho
  ) %>%
  column_to_rownames("plot_label") %>%
  as.matrix()

mat_star <- corr_filt %>%
  select(plot_label, physio, stars) %>%
  pivot_wider(
    names_from = physio,
    values_from = stars
  ) %>%
  column_to_rownames("plot_label") %>%
  as.matrix()

# Add missing physiological traits as empty columns if needed.
missing_physio <- setdiff(physio_cols, colnames(mat_rho))

if (length(missing_physio) > 0) {
  for (trait in missing_physio) {
    mat_rho <- cbind(
      mat_rho,
      setNames(rep(NA_real_, nrow(mat_rho)), trait)
    )
    
    mat_star <- cbind(
      mat_star,
      setNames(rep("", nrow(mat_star)), trait)
    )
  }
}

# Apply the required physiological-trait order.
mat_rho <- mat_rho[, physio_cols, drop = FALSE]
mat_star <- mat_star[, physio_cols, drop = FALSE]

# Neutral display for missing correlations; no star for missing tests.
mat_rho[is.na(mat_rho)] <- 0
mat_star[is.na(mat_star)] <- ""

# Rename columns only for display.
colnames(mat_rho) <- unname(physio_labels[colnames(mat_rho)])
colnames(mat_star) <- unname(physio_labels[colnames(mat_star)])

# ---------- HEATMAP SETTINGS --------------------------------------------------

col_fun <- colorRamp2(
  breaks = c(-1, -0.5, 0, 0.5, 1),
  colors = c("#2166AC", "#67A9CF", "white", "#EF8A62", "#B2182B")
)

n_asv <- nrow(mat_rho)

pdf_height <- max(7, 0.35 * n_asv + 3)
tiff_height <- 5

ht <- Heatmap(
  mat_rho,
  name = "Spearman \u03c1",
  col = col_fun,
  
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  
  show_row_dend = TRUE,
  row_dend_side = "left",
  row_dend_width = unit(18, "mm"),
  show_column_dend = FALSE,
  
  row_names_side = "left",
  row_names_max_width = unit(95, "mm"),
  row_names_gp = gpar(fontsize = 11),
  column_names_gp = gpar(fontsize = 13, fontface = "bold"),
  
  #column_title = plot_title,
  column_title_side = "top",
  column_title_gp = gpar(
    fontsize = 16,
    fontface = "bold",
    fontfamily = "mono"
  ),
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    star <- mat_star[i, j]
    
    if (!is.na(star) && nzchar(star)) {
      grid.text(
        star,
        x = x,
        y = y,
        gp = gpar(
          fontsize = 14,
          fontface = "bold",
          col = "black"
        )
      )
    }
  },
  
  heatmap_legend_param = list(
    direction = "horizontal",
    title_gp = gpar(fontsize = 12, fontface = "bold"),
    labels_gp = gpar(fontsize = 11),
    legend_width = unit(9, "cm"),
    grid_height = unit(5, "mm")
  )
)

# ---------- DISPLAY ------------------------------------------------------------

draw(
  ht,
  heatmap_legend_side = "bottom",
  padding = unit(c(4, 20, 4, 20), "mm")
)

# ---------- EXPORT -------------------------------------------------------------

pdf(
  pdf_file,
  width = 8.5,
  height = pdf_height,
  useDingbats = FALSE
)

draw(
  ht,
  heatmap_legend_side = "bottom",
  padding = unit(c(4, 20, 4, 20), "mm")
)

dev.off()

tiff(
  tiff_file,
  width = 8.5,
  height = tiff_height,
  units = "in",
  res = 1200,
  compression = "lzw"
)

draw(
  ht,
  heatmap_legend_side = "bottom",
  padding = unit(c(4, 20, 4, 20), "mm")
)

dev.off()

message("Saved:")
message("  ", pdf_file)
message("  ", tiff_file)

