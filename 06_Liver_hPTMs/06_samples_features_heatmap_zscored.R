# Libraries
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(grid)

# Settings
INPUT_RDS <- "05_batch_correction_limma/H3_H4_8PTM_liver_batch_corrected_limma.rds"
METADATA_FILE <- "STPN2309_metadata_all.tsv"

OUTPUT_DIR <- "06_samples_features_heatmap_zscored"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

SAMPLE_ID_COL <- "sampleID_STPN2309"
DIET_COL <- "diet"
DAY_COL <- "day_categ"

DIETS_TO_KEEP <- c("NC", "HC")

HEATMAP_LEGEND_TITLE <- "Protein-adjusted Histone PTM intensities (z-score)"

DO_ZSCORE <- TRUE

APPLY_CAP <- TRUE
CAP_QUANTILE <- 0.99

CLUSTER_COLUMNS <- TRUE
CLUSTER_ROWS <- TRUE
SAMPLE_COR_METHOD <- "spearman"
PTM_COR_METHOD <- "pearson"
HCLUST_METHOD <- "complete"

PDF_WIDTH <- 12
PDF_HEIGHT <- 8
TIFF_WIDTH <- 12
TIFF_HEIGHT <- 8
TIFF_RES <- 1200

NA_HEATMAP_COL <- "grey85"

# Helper functions
row_zscore_na <- function(mat) {
  
  z <- t(apply(mat, 1, function(x) {
    
    observed <- !is.na(x)
    
    if (sum(observed) < 2) {
      return(rep(NA_real_, length(x)))
    }
    
    mu <- mean(x[observed])
    s <- sd(x[observed])
    
    if (is.na(s) || s < 1e-8) {
      return(rep(NA_real_, length(x)))
    }
    
    out <- rep(NA_real_, length(x))
    out[observed] <- (x[observed] - mu) / s
    
    out
  }))
  
  rownames(z) <- rownames(mat)
  colnames(z) <- colnames(mat)
  
  z
}

# Check input files
if (!file.exists(INPUT_RDS)) {
  stop("Corrected hPTM RDS not found: ", INPUT_RDS)
}

if (!file.exists(METADATA_FILE)) {
  stop("Metadata TSV not found: ", METADATA_FILE)
}

# Load corrected hPTM matrix and metadata
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
  stop("Corrected hPTM matrix must contain sample IDs as rownames.")
}

required_cols <- c(SAMPLE_ID_COL, DIET_COL, DAY_COL)

missing_cols <- setdiff(required_cols, names(meta_df))

if (length(missing_cols) > 0) {
  stop(
    "Missing required metadata columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

# Align matrix and metadata
rownames(ptm_df) <- trimws(rownames(ptm_df))

meta_df[[SAMPLE_ID_COL]] <- trimws(
  as.character(meta_df[[SAMPLE_ID_COL]])
)

meta_df[[DIET_COL]] <- trimws(
  as.character(meta_df[[DIET_COL]])
)

meta_df[[DAY_COL]] <- trimws(
  as.character(meta_df[[DAY_COL]])
)

if (anyDuplicated(rownames(ptm_df)) > 0) {
  stop("Duplicated sample IDs found in corrected hPTM matrix.")
}

if (anyDuplicated(meta_df[[SAMPLE_ID_COL]]) > 0) {
  stop("Duplicated sample IDs found in metadata.")
}

meta_df <- meta_df %>%
  filter(.data[[DIET_COL]] %in% DIETS_TO_KEEP) %>%
  droplevels()

common_ids <- rownames(ptm_df)[
  rownames(ptm_df) %in% meta_df[[SAMPLE_ID_COL]]
]

if (length(common_ids) == 0) {
  stop("No matching sample IDs found between corrected matrix and metadata.")
}

ptm_df <- ptm_df[common_ids, , drop = FALSE]

meta_df <- meta_df[
  match(common_ids, meta_df[[SAMPLE_ID_COL]]),
  ,
  drop = FALSE
]

rownames(meta_df) <- meta_df[[SAMPLE_ID_COL]]

stopifnot(
  identical(
    rownames(ptm_df),
    rownames(meta_df)
  )
)

meta_df[[DIET_COL]] <- factor(
  meta_df[[DIET_COL]],
  levels = DIETS_TO_KEEP
)

cat("Aligned hPTM matrix:", nrow(ptm_df), "samples x", ncol(ptm_df), "PTMs\n")
print(table(meta_df[[DIET_COL]]))

# Prepare PTM matrix
ptm_df[] <- lapply(
  ptm_df,
  function(x) suppressWarnings(as.numeric(x))
)

mat <- t(as.matrix(ptm_df))

stopifnot(
  identical(
    colnames(mat),
    meta_df[[SAMPLE_ID_COL]]
  )
)

# Z-score PTMs across samples
if (DO_ZSCORE) {
  mat_z <- row_zscore_na(mat)
} else {
  mat_z <- mat
}

# Compute display cap
z_vals <- as.numeric(mat_z)
z_vals <- z_vals[!is.na(z_vals)]

z_cap <- as.numeric(
  quantile(
    abs(z_vals),
    probs = CAP_QUANTILE,
    na.rm = TRUE
  )
)

if (APPLY_CAP) {
  
  mat_plot <- mat_z
  
  mat_plot[mat_plot > z_cap] <- z_cap
  mat_plot[mat_plot < -z_cap] <- -z_cap
  
} else {
  
  mat_plot <- mat_z
  z_cap <- max(abs(mat_plot), na.rm = TRUE)
}

# Save matrix used for the heatmap
write.table(
  mat_plot,
  file = file.path(
    OUTPUT_DIR,
    "histone_ptm_matrix_for_heatmap.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

write.table(
  data.frame(
    cap_quantile = CAP_QUANTILE,
    chosen_z_cap = z_cap,
    cells_above_cap = sum(mat_z > z_cap, na.rm = TRUE),
    cells_below_cap = sum(mat_z < -z_cap, na.rm = TRUE)
  ),
  file = file.path(
    OUTPUT_DIR,
    "cap_analysis.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Create metadata annotations
day_num <- as.numeric(
  gsub(
    "^day_0?",
    "",
    as.character(meta_df[[DAY_COL]])
  )
)

meta_df$temporal_phase <- ifelse(
  day_num %in% c(1, 2, 3, 4),
  "Early",
  "Late"
)

meta_df$temporal_phase <- factor(
  meta_df$temporal_phase,
  levels = c("Early", "Late")
)

meta_df$diet_time_group <- paste(
  meta_df[[DIET_COL]],
  meta_df$temporal_phase,
  sep = " - "
)

meta_df$diet_time_group <- factor(
  meta_df$diet_time_group,
  levels = c(
    "NC - Early",
    "NC - Late",
    "HC - Early",
    "HC - Late"
  )
)

day_lab <- paste0("Day ", day_num)

day_lab <- factor(
  day_lab,
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

annotation_df <- data.frame(
  sample_id = meta_df[[SAMPLE_ID_COL]],
  diet = meta_df[[DIET_COL]],
  day = day_lab,
  temporal_phase = meta_df$temporal_phase,
  diet_time_group = meta_df$diet_time_group,
  stringsAsFactors = FALSE
)

write.table(
  annotation_df,
  file = file.path(
    OUTPUT_DIR,
    "sample_annotations_for_heatmap.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Annotation colors
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

diet_time_cols <- c(
  "NC - Early" = "#6baed6",
  "NC - Late" = "#2171b5",
  "HC - Early" = "#fdae6b",
  "HC - Late" = "#d94801"
)

# Heatmap colors
col_fun <- colorRamp2(
  c(-z_cap, 0, z_cap),
  rev(
    colorRampPalette(
      brewer.pal(11, "RdBu")
    )(256)
  )[c(1, 128, 256)]
)

legend_breaks <- round(c(-z_cap, 0, z_cap), 2)

# Cluster samples
sample_cor <- cor(
  mat_z,
  method = SAMPLE_COR_METHOD,
  use = "pairwise.complete.obs"
)

sample_hc <- hclust(
  as.dist(1 - sample_cor),
  method = HCLUST_METHOD
)

# Heatmap annotation
ha_all <- HeatmapAnnotation(
  diet = meta_df[[DIET_COL]],
  day = day_lab,
  temporal_phase = meta_df$temporal_phase,
  diet_time_group = meta_df$diet_time_group,
  col = list(
    diet = diet_cols,
    day = day_cols,
    temporal_phase = phase_cols,
    diet_time_group = diet_time_cols
  ),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    diet = list(
      title = "Diet",
      title_gp = gpar(fontsize = 16),
      labels_gp = gpar(fontsize = 14),
      grid_height = unit(0.5, "cm"),
      grid_width = unit(0.5, "cm")
    ),
    day = list(
      title = "Day",
      title_gp = gpar(fontsize = 16),
      labels_gp = gpar(fontsize = 14),
      grid_height = unit(0.5, "cm"),
      grid_width = unit(0.5, "cm")
    ),
    temporal_phase = list(
      title = "Temporal Phase",
      title_gp = gpar(fontsize = 16),
      labels_gp = gpar(fontsize = 14),
      grid_height = unit(0.5, "cm"),
      grid_width = unit(0.5, "cm")
    ),
    diet_time_group = list(
      title = "Diet × Temporal Phase",
      title_gp = gpar(fontsize = 16),
      labels_gp = gpar(fontsize = 14),
      grid_height = unit(0.5, "cm"),
      grid_width = unit(0.5, "cm")
    )
  )
)

# Build heatmap
ht_ptm_samples <- Heatmap(
  mat_plot,
  name = HEATMAP_LEGEND_TITLE,
  col = col_fun,
  na_col = NA_HEATMAP_COL,
  top_annotation = ha_all,
  
  row_split = NULL,
  cluster_row_slices = FALSE,
  
  height = unit(14, "cm"),
  
  use_raster = TRUE,
  raster_device = "png",
  raster_quality = 5,
  
  cluster_columns = if (CLUSTER_COLUMNS) {
    as.dendrogram(sample_hc)
  } else {
    FALSE
  },
  
  cluster_rows = CLUSTER_ROWS,
  
  clustering_distance_rows = function(x) {
    as.dist(
      1 - cor(
        t(x),
        method = PTM_COR_METHOD,
        use = "pairwise.complete.obs"
      )
    )
  },
  
  clustering_method_rows = HCLUST_METHOD,
  
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 10),
  row_names_max_width = unit(8, "cm"),
  
  column_title = NULL,
  row_title = "Histone PTMs",
  
  heatmap_legend_param = list(
    title = HEATMAP_LEGEND_TITLE,
    at = legend_breaks,
    labels = legend_breaks,
    title_gp = gpar(fontsize = 18, fontface = "plain"),
    labels_gp = gpar(fontsize = 15),
    direction = "horizontal",
    legend_width = unit(10, "cm"),
    grid_width = unit(0.7, "cm"),
    grid_height = unit(0.7, "cm"),
    title_position = "topcenter"
  )
)

# Save heatmap
pdf(
  file = file.path(
    OUTPUT_DIR,
    "histone_ptm_sample_heatmap_zscore.pdf"
  ),
  width = PDF_WIDTH,
  height = PDF_HEIGHT,
  useDingbats = FALSE
)

draw(
  ht_ptm_samples,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "right",
  align_heatmap_legend = "heatmap_center",
  merge_legends = FALSE
)

dev.off()

tiff(
  filename = file.path(
    OUTPUT_DIR,
    "histone_ptm_sample_heatmap_zscore.tiff"
  ),
  width = TIFF_WIDTH,
  height = TIFF_HEIGHT,
  units = "in",
  res = TIFF_RES,
  compression = "lzw"
)

draw(
  ht_ptm_samples,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "right",
  align_heatmap_legend = "heatmap_center",
  merge_legends = FALSE
)

dev.off()

writeLines(
  capture.output(sessionInfo()),
  con = file.path(OUTPUT_DIR, "sessionInfo.txt")
)

cat("Done. Outputs saved in:", OUTPUT_DIR, "\n")
