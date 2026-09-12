# libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(matrixStats)
  library(ggplot2)
  library(patchwork)
})

set.seed(42)
options(stringsAsFactors = FALSE)

# 1. Input and output paths



input_dir <- "input_files"

metadata_file <- "metadata_for_omics_data.rds"

liver_hptm_file <- file.path(
  input_dir,
  "H3_H4_8PTM_liver_batch_corrected_limma.rds"
)

midgut_hptm_file <- file.path(
  input_dir,
  "H3_H4_5PTM_midgut_batch_corrected_limma.rds"
)

digesta_microbiota_file <- file.path(
  input_dir,
  "digesta_ASVbest_sample_by_feature.rds"
)

mucus_microbiota_file <- file.path(
  input_dir,
  "mucus_ASVbest_sample_by_feature.rds"
)

midgut_proteome_file <- file.path(
  input_dir,
  "midgut_proteome_filtered_batch_corrected_limma.rds"
)

liver_proteome_file <- file.path(
  input_dir,
  "liver_proteome_filtered_batch_corrected_limma.rds"
)

out_dir <- "QC_views"
fig_dir <- file.path(out_dir, "figures")
pdf_dir <- file.path(fig_dir, "PDF")
tiff_dir <- file.path(fig_dir, "TIFF_1200dpi")
tab_dir <- file.path(out_dir, "tables")
rds_dir <- file.path(out_dir, "rds")

dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tiff_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(
  metadata_file,
  liver_hptm_file,
  midgut_hptm_file,
  digesta_microbiota_file,
  mucus_microbiota_file,
  midgut_proteome_file,
  liver_proteome_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing input file(s):\n",
    paste(missing_files, collapse = "\n")
  )
}



# 2. Helpers



as_numeric_matrix <- function(x) {
  mat <- as.matrix(x)
  storage.mode(mat) <- "numeric"
  mat[!is.finite(mat)] <- NA_real_
  mat
}

# Input matrices are expected to be samples x features.

# MOFA-style matrices are returned as features x samples.

to_feature_sample_matrix <- function(x) {
  x <- as_numeric_matrix(x)
  
  if (is.null(rownames(x)) || is.null(colnames(x))) {
    stop("Each omics matrix must have row names and column names.")
  }
  
  t(x)
}

# Keep only samples present in both metadata and the current omics view.

# This is performed independently for every view.

keep_nonfasted_samples <- function(mat, metadata_df) {
  metadata_samples <- rownames(metadata_df)
  
  retained_samples <- metadata_samples[
    metadata_samples %in% colnames(mat)
  ]
  
  mat[, retained_samples, drop = FALSE]
}

zscore_features <- function(mat, min_observed = 3) {
  mat <- as_numeric_matrix(mat)
  
  n_observed <- rowSums(!is.na(mat))
  feature_mean <- rowMeans(mat, na.rm = TRUE)
  feature_sd <- matrixStats::rowSds(mat, na.rm = TRUE)
  
  keep <- n_observed >= min_observed &
    is.finite(feature_sd) &
    feature_sd > 0
  
  dropped_features <- data.frame(
    feature = rownames(mat)[!keep],
    n_observed = n_observed[!keep],
    sd = feature_sd[!keep]
  )
  
  mat <- mat[keep, , drop = FALSE]
  feature_mean <- feature_mean[keep]
  feature_sd <- feature_sd[keep]
  
  z_mat <- sweep(mat, 1, feature_mean, "-")
  z_mat <- sweep(z_mat, 1, feature_sd, "/")
  z_mat[!is.finite(z_mat)] <- NA_real_
  
  attr(z_mat, "dropped_features") <- dropped_features
  
  z_mat
}

safe_min <- function(x) {
  x <- x[is.finite(x)]
  
  if (length(x) == 0) {
    return(NA_real_)
  }
  
  min(x)
}

safe_max <- function(x) {
  x <- x[is.finite(x)]
  
  if (length(x) == 0) {
    return(NA_real_)
  }
  
  max(x)
}

safe_quantile <- function(x, probability) {
  x <- x[is.finite(x)]
  
  if (length(x) == 0) {
    return(NA_real_)
  }
  
  as.numeric(
    stats::quantile(
      x,
      probs = probability,
      names = FALSE,
      type = 7
    )
  )
}

row_skewness <- function(x) {
  x <- x[is.finite(x)]
  
  if (length(x) < 3) {
    return(NA_real_)
  }
  
  x_sd <- stats::sd(x)
  
  if (!is.finite(x_sd) || x_sd == 0) {
    return(NA_real_)
  }
  
  mean((x - mean(x))^3) / x_sd^3
}

theme_publication <- function(base_size = 13) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      strip.background = element_rect(fill = "grey95", colour = "grey70"),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.title = element_text(face = "bold")
    )
}

save_publication_plot <- function(plot, filename, width, height, dpi = 1200) {
  pdf_file <- file.path(pdf_dir, paste0(filename, ".pdf"))
  tif_file <- file.path(tiff_dir, paste0(filename, ".tiff"))
  
  ggsave(
    filename = pdf_file,
    plot = plot,
    width = width,
    height = height,
    units = "in"
  )
  
  ggsave(
    filename = tif_file,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    compression = "lzw"
  )
  
  message("Saved PDF: ", pdf_file)
  message("Saved TIFF: ", tif_file)
}



# 3. Load metadata and non-omics covariates



metadata_data <- readRDS(metadata_file)

metadata <- metadata_data$metadata
host_traits_all <- as.data.frame(metadata_data$host_traits)
methylome <- as.data.frame(metadata_data$Liver_DNA_methylation)

required_metadata_cols <- c(
  "day_categ",
  "diet",
  "group",
  "diet_phase"
)

if (!all(required_metadata_cols %in% colnames(metadata))) {
  stop(
    "Missing metadata column(s): ",
    paste(
      setdiff(required_metadata_cols, colnames(metadata)),
      collapse = ", "
    )
  )
}

metadata <- metadata[, required_metadata_cols, drop = FALSE] |>
  mutate(across(everything(), as.factor))

methylome <- methylome |>
  select(dC_rel, mdC_rel, hmdC_rel)

if (is.null(rownames(metadata))) {
  stop("Metadata must have sample IDs as row names.")
}

if (is.null(rownames(host_traits_all))) {
  stop("Host traits must have sample IDs as row names.")
}

if (is.null(rownames(methylome))) {
  stop("Methylome data must have sample IDs as row names.")
}

# This keeps samples present in metadata and host traits.

# Methylome values may remain NA for samples not measured in that dataset.

metadata_nonfasted <- metadata[
  rownames(metadata) %in% rownames(host_traits_all),
  ,
  drop = FALSE
]

host_traits_nonfasted <- host_traits_all[
  rownames(metadata_nonfasted),
  ,
  drop = FALSE
]

methylome_nonfasted <- methylome[
  rownames(metadata_nonfasted),
  ,
  drop = FALSE
]

metadata_nonfasted <- cbind(
  metadata_nonfasted,
  host_traits_nonfasted,
  methylome_nonfasted
)

write.csv(
  metadata_nonfasted |>
    rownames_to_column("sample_id"),
  file.path(tab_dir, "metadata_nonfasted_for_MOFA_QC.csv"),
  row.names = FALSE
)



# 4. Load omics data
extract_rds_element <- function(file, element_name) {
  rds_object <- readRDS(file)
  
  if (!is.list(rds_object)) {
    stop(
      "Expected a list in: ", file,
      "\nBut readRDS() returned: ", class(rds_object)[1]
    )
  }
  
  if (!element_name %in% names(rds_object)) {
    stop(
      "Element '", element_name, "' was not found in:\n", file,
      "\nAvailable elements: ",
      paste(names(rds_object), collapse = ", ")
    )
  }
  
  extracted_object <- rds_object[[element_name]]
  
  if (!is.matrix(extracted_object) && !is.data.frame(extracted_object)) {
    stop(
      "Element '", element_name, "' in ", file,
      " is not a matrix or data frame."
    )
  }
  
  extracted_object
}

# hPTM matrices: RDS files are lists
liver_hPTM <- extract_rds_element(
  file = liver_hptm_file,
  element_name = "ptm_batch_corrected_limma"
)
midgut_hPTM <- extract_rds_element(
  file = midgut_hptm_file,
  element_name = "ptm_batch_corrected_limma"
)

# Proteome matrices: RDS files are lists
midgut_proteome <- extract_rds_element(
  file = midgut_proteome_file,
  element_name = "proteome_batch_corrected_limma"
)
liver_proteome <- extract_rds_element(
  file = liver_proteome_file,
  element_name = "proteome_batch_corrected_limma"
)

# Microbiota files are expected to contain a matrix/data frame directly
digesta_microbiota <- readRDS(digesta_microbiota_file)
mucus_microbiota <- readRDS(mucus_microbiota_file)

if (!is.matrix(digesta_microbiota) && !is.data.frame(digesta_microbiota)) {
  stop("digesta_microbiota RDS must contain a matrix or data frame directly.")
}

if (!is.matrix(mucus_microbiota) && !is.data.frame(mucus_microbiota)) {
  stop("mucus_microbiota RDS must contain a matrix or data frame directly.")
}

message("Loaded omics datasets:")
message("  Liver hPTM: ", nrow(liver_hPTM), " samples × ", ncol(liver_hPTM), " features")
message("  Midgut hPTM: ", nrow(midgut_hPTM), " samples × ", ncol(midgut_hPTM), " features")
message("  Midgut proteome: ", nrow(midgut_proteome), " samples × ", ncol(midgut_proteome), " features")
message("  Liver proteome: ", nrow(liver_proteome), " samples × ", ncol(liver_proteome), " features")
message("  Digesta microbiota: ", nrow(digesta_microbiota), " samples × ", ncol(digesta_microbiota), " features")
message("  Mucus microbiota: ", nrow(mucus_microbiota), " samples × ", ncol(mucus_microbiota), " features")




# 5. Build independent omics views
# Each dataset retains only its own valid non-fasted samples.
# No global intersection across omics datasets is performed here.



view_list_original <- list(
  "Midgut hPTMs" = to_feature_sample_matrix(midgut_hPTM),
  "Liver hPTMs" = to_feature_sample_matrix(liver_hPTM),
  "Midgut proteome" = to_feature_sample_matrix(midgut_proteome),
  "Liver proteome" = to_feature_sample_matrix(liver_proteome),
  "Digesta microbiota" = to_feature_sample_matrix(digesta_microbiota),
  "Mucus microbiota" = to_feature_sample_matrix(mucus_microbiota)
)

view_list_original <- lapply(
  view_list_original,
  keep_nonfasted_samples,
  metadata_df = metadata_nonfasted
)

empty_views <- names(view_list_original)[
  vapply(view_list_original, ncol, integer(1)) == 0
]

if (length(empty_views) > 0) {
  stop(
    "No metadata-matched samples found for: ",
    paste(empty_views, collapse = ", ")
  )
}

view_dimensions <- data.frame(
  view = names(view_list_original),
  n_features = vapply(view_list_original, nrow, integer(1)),
  n_samples = vapply(view_list_original, ncol, integer(1))
)

write.csv(
  view_dimensions,
  file.path(tab_dir, "omics_view_dimensions_original.csv"),
  row.names = FALSE
)



# 6. Feature-wise z-scoring
view_list_zscore <- lapply(
  view_list_original,
  zscore_features
)

zscore_report <- bind_rows(
  lapply(names(view_list_original), function(view_name) {
    dropped <- attr(view_list_zscore[[view_name]], "dropped_features")
    
    data.frame(
      view = view_name,
      n_features_original = nrow(view_list_original[[view_name]]),
      n_features_zscore = nrow(view_list_zscore[[view_name]]),
      n_features_dropped = nrow(dropped)
    )
    
  })
)

write.csv(
  zscore_report,
  file.path(tab_dir, "zscore_feature_filtering_report.csv"),
  row.names = FALSE
)

saveRDS(
  view_list_original,
  file.path(rds_dir, "view_list_original_nonfasted.rds")
)

saveRDS(
  view_list_zscore,
  file.path(rds_dir, "view_list_zscore_nonfasted.rds")
)



# 7. QC tables
make_value_df <- function(view_list,
                          transform_label,
                          max_values_per_view = 200000) {
  bind_rows(
    lapply(names(view_list), function(view_name) {
      values <- as.numeric(view_list[[view_name]])
      values <- values[is.finite(values)]
      
      
      if (length(values) > max_values_per_view) {
        values <- sample(values, max_values_per_view)
      }
      
      data.frame(
        transform = transform_label,
        view = view_name,
        value = values
      )
    })
    
    
  )
}

make_feature_stats <- function(mat, view_name, transform_label) {
  mat <- as_numeric_matrix(mat)
  
  data.frame(
    transform = transform_label,
    view = view_name,
    feature = rownames(mat),
    n_samples = ncol(mat),
    n_observed = rowSums(!is.na(mat)),
    missing_fraction = rowMeans(is.na(mat)),
    zero_fraction = rowMeans(mat == 0, na.rm = TRUE),
    mean = rowMeans(mat, na.rm = TRUE),
    median = matrixStats::rowMedians(mat, na.rm = TRUE),
    variance = matrixStats::rowVars(mat, na.rm = TRUE),
    sd = matrixStats::rowSds(mat, na.rm = TRUE),
    min = apply(mat, 1, safe_min),
    q05 = apply(mat, 1, safe_quantile, probability = 0.05),
    q95 = apply(mat, 1, safe_quantile, probability = 0.95),
    max = apply(mat, 1, safe_max),
    skewness = apply(mat, 1, row_skewness)
  )
}

make_variance_rank_df <- function(feature_stats) {
  feature_stats |>
    filter(is.finite(variance), variance > 0) |>
    group_by(transform, view) |>
    arrange(desc(variance), .by_group = TRUE) |>
    mutate(
      feature_rank = row_number(),
      feature_fraction = feature_rank / n(),
      cumulative_variance_fraction = cumsum(variance) / sum(variance)
    ) |>
    ungroup()
}

feature_stats_original <- bind_rows(
  lapply(names(view_list_original), function(view_name) {
    make_feature_stats(
      mat = view_list_original[[view_name]],
      view_name = view_name,
      transform_label = "Original"
    )
  })
)

feature_stats_zscore <- bind_rows(
  lapply(names(view_list_zscore), function(view_name) {
    make_feature_stats(
      mat = view_list_zscore[[view_name]],
      view_name = view_name,
      transform_label = "Z-score"
    )
  })
)

qc_values <- bind_rows(
  make_value_df(view_list_original, "Original"),
  make_value_df(view_list_zscore, "Z-score")
)

qc_feature_stats <- bind_rows(
  feature_stats_original,
  feature_stats_zscore
)

qc_variance_rank <- make_variance_rank_df(qc_feature_stats)

write.csv(
  qc_feature_stats,
  file.path(tab_dir, "feature_level_QC_original_vs_zscore.csv"),
  row.names = FALSE
)

write.csv(
  qc_variance_rank,
  file.path(tab_dir, "variance_rank_original_vs_zscore.csv"),
  row.names = FALSE
)



# 8. QC figures ---------------------------------------------------------------

theme_qc <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "grey92", colour = NA),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      plot.margin = margin(6, 8, 6, 6)
    )
}

p01_value_density <- ggplot(qc_values, aes(value)) +
  geom_density(linewidth = 0.8, na.rm = TRUE) +
  facet_grid(transform ~ view, scales = "free") +
  theme_qc() +
  labs(x = "Feature value", y = "Density")

p02_feature_sd <- qc_feature_stats |>
  filter(is.finite(sd), sd > 0) |>
  ggplot(aes(sd)) +
  geom_density(linewidth = 0.8) +
  scale_x_log10() +
  facet_grid(transform ~ view, scales = "free_y") +
  theme_qc() +
  labs(x = "Feature SD across samples", y = "Density")

p03_cumulative_variance <- ggplot(
  qc_variance_rank,
  aes(feature_fraction, cumulative_variance_fraction)
) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_line(linewidth = 0.9) +
  facet_grid(transform ~ view) +
  theme_qc() +
  labs(
    x = "Fraction of features ranked by variance",
    y = "Cumulative fraction of total variance"
  )

p04_zscore_density <- qc_values |>
  filter(transform == "Z-score") |>
  ggplot(aes(value, colour = view)) +
  geom_density(linewidth = 0.9, na.rm = TRUE) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
  theme_qc() +
  labs(
    x = "Feature-wise z-score",
    y = "Density",
    colour = "Omic view"
  )

# Combined QC figure: no title, no panel letters/tags
p_QC_combined <- (
  (p01_value_density | p02_feature_sd) /
    (p03_cumulative_variance | p04_zscore_density)
) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

print(p_QC_combined)

plots_to_save <- list(
  "01_value_density_original_vs_zscore" = list(p01_value_density, 10, 4.5),
  "02_feature_SD_original_vs_zscore" = list(p02_feature_sd, 10, 4.5),
  "03_cumulative_variance_original_vs_zscore" = list(p03_cumulative_variance, 12, 4.5),
  "04_zscore_density_overlap_across_omics" = list(p04_zscore_density, 7, 4.5),
  "05_combined_omics_QC" = list(p_QC_combined, 18, 11)
)

invisible(lapply(names(plots_to_save), function(filename) {
  x <- plots_to_save[[filename]]
  save_publication_plot(
    plot = x[[1]],
    filename = filename,
    width = x[[2]],
    height = x[[3]]
  )
}))


message("QC output directory: ", out_dir)
