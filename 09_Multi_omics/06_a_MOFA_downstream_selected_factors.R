suppressPackageStartupMessages({
  library(MOFA2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(data.table)
  library(purrr)
  library(ggpubr)
  library(readr)
  library(tibble)
  library(ggrepel)
  library(forcats)
  library(cowplot)
  library(scales)
  library(grid)
})

set.seed(42)

# Input files
model_file <- "mofa_models/K20_setseed_123.hdf5"
protein_info_file <- "STRG0A55HWH.protein.info.v12.0.txt"
string_terms_file <- "110079946.protein.enrichment.terms.v12.0.txt"

# Output directories
out_dir <- "mofa_results"
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
rds_dir <- file.path(out_dir, "rds")
txt_dir <- file.path(out_dir, "text")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(txt_dir, recursive = TRUE, showWarnings = FALSE)

# General settings
cols2diet <- c(
  "NC" = "#44a7c4",
  "HC" = "#f99943"
)

mono_bold_title <- theme(
  plot.title = element_text(
    family = "mono",
    face = "bold"
  )
)

theme_mofa <- theme_bw(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    panel.grid.minor = element_blank()
  )


safe_output_file <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  if (!file.exists(path)) {
    return(path)
  }

  removed <- suppressWarnings(file.remove(path))

  if (isTRUE(removed) && !file.exists(path)) {
    return(path)
  }

  extension <- tools::file_ext(path)
  stem <- tools::file_path_sans_ext(basename(path))
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

  fallback <- file.path(
    dirname(path),
    paste0(stem, "_", timestamp, ".", extension)
  )

  warning(
    "The requested output file is open or locked. Saving as: ",
    fallback
  )

  fallback
}

save_plot_to <- function(plot,
                         filename,
                         width,
                         height,
                         output_dir,
                         dpi = 1200) {

  if (inherits(plot, "ggplot")) {
    plot <- plot +
      theme(
        plot.title = element_text(
          family = "mono",
          face = "bold"
        ),
        plot.background = element_rect(
          fill = "transparent",
          colour = NA
        ),
        panel.background = element_rect(
          fill = "transparent",
          colour = NA
        ),
        legend.background = element_rect(
          fill = "transparent",
          colour = NA
        ),
        legend.box.background = element_rect(
          fill = "transparent",
          colour = NA
        )
      )
  }

  pdf_file <- safe_output_file(
    file.path(output_dir, paste0(filename, ".pdf"))
  )

  tiff_file <- safe_output_file(
    file.path(output_dir, paste0(filename, ".tiff"))
  )

  ggsave(
    filename = pdf_file,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = "pdf",
    bg = "transparent"
  )

  ggsave(
    filename = tiff_file,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    compression = "lzw",
    bg = "transparent"
  )

  message("Saved: ", pdf_file)
  message("Saved: ", tiff_file)
}

save_plot <- function(plot, filename, width, height, dpi = 1200) {
  save_plot_to(
    plot = plot,
    filename = filename,
    width = width,
    height = height,
    output_dir = fig_dir,
    dpi = dpi
  )
}

save_gtable_plot_to <- function(plot_object,
                                filename,
                                width,
                                height,
                                output_dir,
                                dpi = 1200) {

  pdf_file <- safe_output_file(
    file.path(output_dir, paste0(filename, ".pdf"))
  )

  tiff_file <- safe_output_file(
    file.path(output_dir, paste0(filename, ".tiff"))
  )

  pdf_device <- NULL

  tryCatch(
    {
      grDevices::pdf(
        pdf_file,
        width = width,
        height = height
      )
      pdf_device <- grDevices::dev.cur()
      grid::grid.newpage()
      grid::grid.draw(plot_object$gtable)
    },
    finally = {
      open_devices <- grDevices::dev.list()
      if (
        !is.null(pdf_device) &&
        !is.null(open_devices) &&
        pdf_device %in% open_devices
      ) {
        grDevices::dev.off(pdf_device)
      }
    }
  )

  tiff_device <- NULL

  tryCatch(
    {
      grDevices::tiff(
        filename = tiff_file,
        width = width,
        height = height,
        units = "in",
        res = dpi,
        compression = "lzw"
      )
      tiff_device <- grDevices::dev.cur()
      grid::grid.newpage()
      grid::grid.draw(plot_object$gtable)
    },
    finally = {
      open_devices <- grDevices::dev.list()
      if (
        !is.null(tiff_device) &&
        !is.null(open_devices) &&
        tiff_device %in% open_devices
      ) {
        grDevices::dev.off(tiff_device)
      }
    }
  )

  message("Saved: ", pdf_file)
  message("Saved: ", tiff_file)
}

save_gtable_plot <- function(plot_object,
                             filename,
                             width,
                             height,
                             dpi = 1200) {
  save_gtable_plot_to(
    plot_object = plot_object,
    filename = filename,
    width = width,
    height = height,
    output_dir = fig_dir,
    dpi = dpi
  )
}

save_heatmap_pdf_to <- function(expr,
                                filename,
                                width,
                                height,
                                output_dir) {

  plot_object <- force(expr)

  pdf_file <- safe_output_file(
    file.path(output_dir, paste0(filename, ".pdf"))
  )

  pdf_device <- NULL

  tryCatch(
    {
      grDevices::pdf(
        file = pdf_file,
        width = width,
        height = height,
        onefile = TRUE
      )

      pdf_device <- grDevices::dev.cur()

      grid::grid.newpage()

      if (inherits(plot_object, "pheatmap")) {
        grid::grid.draw(plot_object$gtable)

      } else if (
        inherits(plot_object, "gtable") ||
        inherits(plot_object, "grob") ||
        inherits(plot_object, "gTree")
      ) {
        grid::grid.draw(plot_object)

      } else if (inherits(plot_object, "ggplot")) {
        print(plot_object)

      } else if (
        inherits(plot_object, "Heatmap") ||
        inherits(plot_object, "HeatmapList")
      ) {
        if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
          stop(
            "The ComplexHeatmap package is required to draw this object."
          )
        }

        ComplexHeatmap::draw(plot_object)

      } else if (
        is.list(plot_object) &&
        !is.null(plot_object$gtable)
      ) {
        grid::grid.draw(plot_object$gtable)

      } else {
        print(plot_object)
      }
    },
    finally = {
      open_devices <- grDevices::dev.list()

      if (
        !is.null(pdf_device) &&
        !is.null(open_devices) &&
        pdf_device %in% open_devices
      ) {
        grDevices::dev.off(pdf_device)
      }
    }
  )

  if (
    !file.exists(pdf_file) ||
    is.na(file.info(pdf_file)$size) ||
    file.info(pdf_file)$size == 0
  ) {
    stop("The heatmap PDF was not written correctly: ", pdf_file)
  }

  message("Saved: ", pdf_file)

  invisible(pdf_file)
}

save_heatmap_pdf <- function(expr, filename, width, height) {
  save_heatmap_pdf_to(
    expr = expr,
    filename = filename,
    width = width,
    height = height,
    output_dir = fig_dir
  )
}

clean_protein_id <- function(x) {
  x |>
    stringr::str_remove("^Midgut_") |>
    stringr::str_remove("^MidGut_") |>
    stringr::str_remove("^Liver_") |>
    stringr::str_remove("^STRG0A55HWH\\.") |>
    stringr::str_remove(";.*$")
}

# Load MOFA model
stopifnot(file.exists(model_file))

model <- load_model(
  model_file,
  load_interpol_Z = FALSE,
  sort_factors = TRUE,
  verbose = TRUE
)

plot_ascii_data(model, nonzero = FALSE)


# Normalize hPTM feature naming inside the MOFA model:
# keep H3/H4 and K/R + residue number as-is,
# lowercase only the chemical mark (Me, Ac, La, Me2, Me3, ...)
# ------------------------------------------------------------
normalize_ptm <- function(x) {
  vapply(x, function(entry) {
    parts <- strsplit(entry, "\\|")[[1]]
    norm_parts <- vapply(parts, function(p) {
      g <- regmatches(
        p,
        regexec("^(Liver_|Midgut_|MidGut_)?(H[34])(K|R)([0-9]{1,3})([A-Za-z0-9]+)$", p)
      )[[1]]
      if (length(g) < 6 || is.na(g[1])) return(NA_character_)
      prefix <- g[2]  # "Liver_", "Midgut_", "MidGut_", or "" if absent
      paste0(prefix, g[3], g[4], g[5], tolower(g[6]))
    }, character(1))
    if (any(is.na(norm_parts))) return(NA_character_)
    paste(norm_parts, collapse = "|")
  }, character(1), USE.NAMES = FALSE)
}

hptm_views <- c("Liver hPTMs", "Midgut hPTMs")
fn <- features_names(model)

for (v in hptm_views) {
  if (!v %in% names(fn)) {
    stop("Expected view not found in model: ", v)
  }
  original   <- fn[[v]]
  normalized <- normalize_ptm(original)
  
  if (any(is.na(normalized))) {
    cat("Unparsed hPTM feature names in view:", v, "\n")
    print(original[is.na(normalized)])
    stop("Some hPTM feature names failed to parse — fix regex or naming before continuing.")
  }
  
  if (anyDuplicated(normalized) > 0) {
    stop("Normalization created duplicate feature names in view: ", v)
  }
  
  fn[[v]] <- normalized
}

features_names(model) <- fn

cat("Liver hPTMs after normalization:\n")
print(features_names(model)[["Liver hPTMs"]])
cat("Midgut hPTMs after normalization:\n")
print(features_names(model)[["Midgut hPTMs"]])


write.csv(
  data.frame(view = views_names(model)),
  file.path(tab_dir, "mofa_view_names.csv"),
  row.names = FALSE
)

# Clean sample metadata
md <- samples_metadata(model)

md <- md |>
  dplyr::rename(
    `5-hmC` = X5.hmC,
    `5-mC` = X5.mC,
    `Unmodified C` = C,
    `Temporal Phase` = phase,
    `Diet x Temporal Phase` = diet_phase,
    `Body Weight` = body_weight,
    `Plasma Glucose` = plasma_glucose,
    `Plasma L-Lactate` = plasma_LLactate,
    `Plasma D-Lactate` = plasma_DLactate,
    Diet = diet,
    Day = day_categ
  ) |>
  dplyr::mutate(
    Day = paste(
      "Day",
      as.numeric(gsub("day_", "", Day))
    ),
    Day = factor(
      Day,
      levels = c(
        "Day 1",
        "Day 2",
        "Day 3",
        "Day 4",
        "Day 10",
        "Day 15",
        "Day 22"
      )
    ),
    
    Diet = factor(
      Diet,
      levels = c("NC", "HC")
    ),
    
    `Temporal Phase` = stringr::str_to_title(
      as.character(`Temporal Phase`)
    ),
    `Temporal Phase` = factor(
      `Temporal Phase`,
      levels = c("Early", "Late")
    ),
    
    `Diet x Temporal Phase` = paste(
      Diet,
      `Temporal Phase`,
      sep = " - "
    ),
    `Diet x Temporal Phase` = factor(
      `Diet x Temporal Phase`,
      levels = c(
        "NC - Early",
        "NC - Late",
        "HC - Early",
        "HC - Late"
      )
    )
  )

samples_metadata(model) <- md

write.csv(
  samples_metadata(model),
  file.path(tab_dir, "mofa_sample_metadata_cleaned.csv"),
  row.names = TRUE
)

saveRDS(
  samples_metadata(model),
  file.path(rds_dir, "mofa_sample_metadata_cleaned.rds")
)

# Data overview
data_plot <- plot_data_overview(
  model,
  show_covariate = FALSE,
  show_dimensions = TRUE
) +
  theme(
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size = 25, face = "bold"),
    strip.text.x = element_text(size = 30, face = "bold"),
    plot.subtitle = element_text(size = 25),
    legend.title = element_text(size = 13, face = "bold")
  )

print(data_plot)

save_plot(
  data_plot,
  filename = "data_overview",
  width = 7.5,
  height = 5
)

# Factor correlation
p_factor_corr <- plot_factor_cor(model)
p_factor_corr

save_heatmap_pdf(
  p_factor_corr,
  filename = "factor_correlation",
  width = 7,
  height = 6
)

# Variance explained
r2 <- model@cache$variance_explained$r2_per_factor[[1]]
r2_total <- model@cache$variance_explained$r2_total[[1]]

r2_dt <- as.data.table(r2)

r2_dt[, factor := as.factor(seq_len(model@dimensions$K))]

r2_dt <- melt(
  r2_dt,
  id.vars = "factor",
  variable.name = "view",
  value.name = "r2"
)

r2_dt[, cum_r2 := cumsum(r2), by = "view"]

write.csv(
  as.data.frame(r2),
  file.path(tab_dir, "variance_explained_per_factor_by_view.csv"),
  row.names = TRUE
)

write.csv(
  as.data.frame(r2_total),
  file.path(tab_dir, "variance_explained_total_by_view.csv"),
  row.names = TRUE
)

saveRDS(
  list(
    r2_per_factor = r2,
    r2_total = r2_total,
    r2_long = as.data.frame(r2_dt)
  ),
  file.path(rds_dir, "variance_explained_results.rds")
)

p_cum_r2 <- ggpubr::ggline(
  r2_dt,
  x = "factor",
  y = "cum_r2",
  color = "view"
) +
  labs(
    x = "Factor number",
    y = "Cumulative variance explained (%)"
  ) +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    axis.text = element_text(size = rel(0.8))
  )

print(p_cum_r2)

save_plot(
  p_cum_r2,
  filename = "cumulative_variance_explained_by_view",
  width = 5,
  height = 4.5
)

r2_total_factor <- rowSums(r2, na.rm = TRUE) |>
  enframe(name = "factor", value = "total_r2") |>
  mutate(
    factor = factor(factor, levels = factor)
  )

write.csv(
  r2_total_factor,
  file.path(tab_dir, "variance_explained_total_per_factor.csv"),
  row.names = FALSE
)

p_r2_total_factor <- ggplot(
  r2_total_factor,
  aes(x = factor, y = total_r2)
) +
  geom_col(
    fill = "#00627d",
    color = "black",
    linewidth = 0.4,
    width = 0.9
  ) +
  theme_classic(base_size = 18) +
  labs(
    x = NULL,
    y = "Variance explained (%)"
  ) +
  theme(
    axis.text.x = element_text(
      size = 20,
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 24, face = "bold"),
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6),
    plot.margin = margin(10, 10, 10, 20)
  )

print(p_r2_total_factor)

save_plot(
  p_r2_total_factor,
  filename = "variance_explained_total_per_factor",
  width = 7.5,
  height = 5
)

p_var <- plot_variance_explained(
  model,
  plot_total = TRUE,
  min_r2 = 0,
  max_r2 = 5,
  factors = "all"
)

var_theme <- theme(
  axis.text.x = element_text(size = 18, angle = 90, vjust = 0.5, hjust = 1),
  axis.text.y = element_text(size = 18),
  axis.title.x = element_text(size = 20, face = "bold"),
  axis.title.y = element_text(size = 18, face = "bold"),
  plot.title = element_text(size = 16, face = "bold"),
  legend.title = element_text(size = 13, face = "bold"),
  legend.text = element_text(size = 12),
  legend.key.size = unit(0.6, "cm")
)

p_var_by_view <- p_var[[1]] + var_theme
p_var_total <- p_var[[2]] + var_theme

print(p_var_by_view)
print(p_var_total)

save_plot(
  p_var_by_view,
  filename = "variance_explained_by_view",
  width = 6,
  height = 6
)

save_plot(
  p_var_total,
  filename = "variance_explained_total",
  width = 4,
  height = 5.5
)

model <- calculate_contribution_scores(model)

calculate_variance_explained_per_sample(
  model,
  views = "all",
  groups = "all",
  factors = "all"
)

txt_file <- file.path(txt_dir, "variance_explained_R2_summary.txt")

old_width <- getOption("width")
old_digits <- getOption("digits")
old_scipen <- getOption("scipen")

options(width = 300, digits = 10, scipen = 999)

sink(txt_file)

cat("> head(model@cache$variance_explained$r2_total[[1]])\n")
print(head(r2_total))

cat("\n\n# total variance per-factor and per-view R2\n")
cat("> model@cache$variance_explained$r2_per_factor[[1]]\n")
print(r2)

cat("\n\n> round(rowSums(r2), 3)\n")
print(round(rowSums(r2), 3))

sink()

options(
  width = old_width,
  digits = old_digits,
  scipen = old_scipen
)

message("Saved variance summary: ", txt_file)

# Factor-covariate correlations
p_corr_host_trait_pval <- correlate_factors_with_covariates(
  model,
  factors = "all",
  covariates = c(
    "Plasma Glucose",
    "Plasma L-Lactate",
    "Plasma D-Lactate",
    "HSI",
    "Unmodified C",
    "5-mC",
    "5-hmC"
  ),
  transpose = FALSE,
  #return_data = TRUE,
  alpha = 0.001,
  plot = "log_pval" # use "r" for coeff or "log_pval"
)

p_corr_host_trait_pval

save_gtable_plot(
  p_corr_host_trait_pval,
  filename = "factor_covariate_correlations_pval",
  width = 3.5,
  height = 6
)

saveRDS(
  p_corr_host_trait_pval,
  file.path(rds_dir, "factor_covariate_correlations_pval.rds")
)

# Factor 2 versus Factor 3
p_f2_f3 <- plot_factors(
  model,
  factors = c(2, 3),
  color_by = "Diet",
  shape_by = "Temporal Phase",
  dot_size = 4
) +
  scale_color_manual(values = cols2diet, name = "Diet") +
  scale_fill_manual(values = cols2diet, name = "Diet") +
  theme_bw(base_size = 18) +
  theme(
    strip.background = element_rect(
      fill = "grey85",
      colour = "grey50",
      linewidth = 0.6
    ),
    strip.text = element_text(size = 18, colour = "black"),
    axis.text = element_text(size = 16, colour = "black"),
    axis.title = element_text(size = 20),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.key.size = unit(0.5, "cm"),
    panel.grid = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.7,
    colour = "black"
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.7,
    colour = "black"
  ) +
  geom_density_2d(
    aes(color = color_by),
    linewidth = 0.5,
    alpha = 0.7,
    show.legend = FALSE
  ) +
  guides(
    shape = guide_legend(title = NULL),
    color = guide_legend(title = NULL)
  )

print(p_f2_f3)

save_plot(
  p_f2_f3,
  filename = "F2_F3",
  width = 6,
  height = 4
)

# UMAP
model <- run_umap(model, factor = "all")

p_umap <- plot_dimred(
  model,
  method = "UMAP",
  shape_by = "Temporal Phase",
  color_by = "Diet",
  label = FALSE,
  legend = TRUE,
  dot_size = 5
) +
  scale_color_manual(values = cols2diet, name = "Diet") +
  scale_fill_manual(values = cols2diet, name = "Diet") +
  geom_density_2d(
    aes(color = .data[["color_by"]]),
    show.legend = FALSE,
    linewidth = 0.6,
    alpha = 0.7
  ) +
  labs(
    x = "UMAP1",
    y = "UMAP2"
  ) +
  theme_bw(base_size = 18) +
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.key.size = unit(0.6, "cm"),
    panel.grid = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) +
  guides(
    shape = guide_legend(title = NULL),
    color = guide_legend(
      title = "Diet",
      override.aes = list(
        shape = 21,
        fill = unname(cols2diet),
        colour = unname(cols2diet),
        size = 5,
        stroke = 0.8
      )
    ),
    fill = guide_legend(title = "Diet")
  )

print(p_umap)

save_plot(
  p_umap,
  filename = "UMAP_diet",
  width = 6,
  height = 4
)

saveRDS(
  model,
  file.path(rds_dir, "mofa_model_with_umap_and_metadata.rds")
)


# Annotation colour palettes
diet_cols <- c(
  "NC" = "#44a7c4",
  "HC" = "#f99943"
)

day_cols <- c(
  "Day 1"  = "#1b9e77",
  "Day 2"  = "#d95f02",
  "Day 3"  = "#7570b3",
  "Day 4"  = "#e7298a",
  "Day 10" = "#66a61e",
  "Day 15" = "#e6ab02",
  "Day 22" = "#a6761d"
)

phase_cols <- c(
  "Early" = "#7b6fd6",
  "Late"  = "#4daf4a"
)

dtg_cols <- c(
  "NC - Early" = "#6baed6",
  "NC - Late"  = "#2171b5",
  "HC - Early" = "#fdae6b",
  "HC - Late"  = "#d94801"
)

annotation_cols <- list(
  Diet = diet_cols,
  Day = day_cols,
  `Temporal Phase` = phase_cols,
  `Diet x Temporal Phase` = dtg_cols
)

heatmap_cols <- grDevices::colorRampPalette(
  c("#2166AC", "#F7F7F7", "#B2182B")
)(100)
# Protein annotations for plotting
stopifnot(file.exists(protein_info_file))

protein_info <- readr::read_tsv(
  protein_info_file,
  show_col_types = FALSE
) |>
  dplyr::rename(string_protein_id = `#string_protein_id`) |>
  dplyr::mutate(
    protein_id = stringr::str_remove(
      string_protein_id,
      "^STRG0A55HWH\\."
    ),
    annotation_clean = dplyr::if_else(
      is.na(annotation) | annotation == "",
      preferred_name,
      annotation
    ),
    annotation_clean = stringr::str_squish(annotation_clean),
    protein_label = paste(
      protein_id,
      annotation_clean,
      sep = "_"
    )
  ) |>
  dplyr::select(
    protein_id,
    protein_label
  )

write.csv(
  protein_info,
  file.path(tab_dir, "protein_annotation_labels.csv"),
  row.names = FALSE
)

saveRDS(
  protein_info,
  file.path(rds_dir, "protein_annotation_labels.rds")
)

add_protein_labels <- function(df, protein_info) {
  df |>
    mutate(
      protein_id = clean_protein_id(feature),
      protein_label = protein_info$protein_label[
        match(protein_id, protein_info$protein_id)
      ],
      feature_plot = ifelse(
        is.na(protein_label),
        protein_id,
        protein_label
      )
    )
}
plot_weights_dot <- function(df, title = "", use_feature_plot = FALSE) {
  df <- df |>
    mutate(
      sign = factor(sign, levels = c("-", "+")),
      value_signed = ifelse(sign == "-", -value, value),
      feature_clean = if (use_feature_plot && "feature_plot" %in% colnames(df)) {
        feature_plot
      } else {
        gsub("^Digesta_g__|^Midgut_|^MidGut_|^Liver_|^Mucus_g__", "", feature)
      }
    ) |>
    arrange(sign, value_signed) |>
    mutate(
      feature_clean = factor(feature_clean, levels = unique(feature_clean))
    )

  ggplot(df, aes(x = value_signed, y = feature_clean, color = sign)) +
    geom_segment(
      aes(
        x = 0,
        xend = value_signed,
        y = feature_clean,
        yend = feature_clean
      ),
      linewidth = 1.4,
      alpha = 0.8
    ) +
    geom_point(size = 7) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.7
    ) +
    scale_color_manual(
      values = c(
        "-" = "#1f4e79",
        "+" = "#8b0000"
      )
    ) +
    labs(
      x = "Weight",
      y = NULL,
      color = "Sign",
      title = title
    ) +
    theme_bw(base_size = 14) +
    theme(
      axis.text.y = element_text(size = 12),
      axis.text.x = element_text(size = 16),
      axis.title.x = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 16),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 16),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
}
# STRING enrichment annotations
stopifnot(file.exists(string_terms_file))

ann <- read_tsv(
  string_terms_file,
  comment = "#",
  col_names = c(
    "string_protein_id",
    "category",
    "term",
    "description"
  ),
  show_col_types = FALSE
) |>
  mutate(
    protein_id = sub("^[^.]+\\.", "", string_protein_id)
  )

ann_bp <- ann |>
  filter(category == "Biological Process (Gene Ontology)")

ann_rp <- ann |>
  filter(category == "Reactome Pathways")

ann_kegg <- ann |>
  filter(stringr::str_detect(category, "KEGG"))

t2g_bp <- ann_bp |>
  transmute(term = term, gene = protein_id) |>
  distinct()

t2n_bp <- ann_bp |>
  select(term, name = description) |>
  distinct()

t2g_kegg <- ann_kegg |>
  transmute(term = term, gene = protein_id) |>
  distinct()

t2n_kegg <- ann_kegg |>
  select(term, name = description) |>
  distinct()

t2g_rp <- ann_rp |>
  transmute(term = term, gene = protein_id) |>
  distinct()

t2n_rp <- ann_rp |>
  select(term, name = description) |>
  distinct()

saveRDS(
  list(
    ann = ann,
    ann_bp = ann_bp,
    ann_kegg = ann_kegg,
    ann_reactome = ann_rp,
    t2g_bp = t2g_bp,
    t2n_bp = t2n_bp,
    t2g_kegg = t2g_kegg,
    t2n_kegg = t2n_kegg,
    t2g_reactome = t2g_rp,
    t2n_reactome = t2n_rp
  ),
  file.path(rds_dir, "STRING_enrichment_annotation_tables.rds")
)

# Gene-set matrices
make_gs_matrix_named <- function(t2g, t2n, model_genes_raw) {
  feature_map <- tibble::tibble(
    gene = clean_protein_id(model_genes_raw),
    model_feature = model_genes_raw
  )

  out <- t2g |>
    inner_join(t2n, by = "term") |>
    mutate(
      label = paste(term, name, sep = ": ")
    ) |>
    inner_join(feature_map, by = "gene") |>
    distinct(label, model_feature) |>
    mutate(value = 1L) |>
    tidyr::pivot_wider(
      names_from = model_feature,
      values_from = value,
      values_fill = 0
    ) |>
    as.data.frame()

  rownames(out) <- out$label
  out$label <- NULL

  as.matrix(out)
}

prot_feats_midgut_raw <- features_names(model)[["Midgut Proteome"]]
prot_feats_liver_raw <- features_names(model)[["Liver Proteome"]]

gs_bp_mat_midgut <- make_gs_matrix_named(
  t2g_bp,
  t2n_bp,
  prot_feats_midgut_raw
)

gs_kegg_mat_midgut <- make_gs_matrix_named(
  t2g_kegg,
  t2n_kegg,
  prot_feats_midgut_raw
)

gs_rp_mat_midgut <- make_gs_matrix_named(
  t2g_rp,
  t2n_rp,
  prot_feats_midgut_raw
)

gs_bp_mat_liver <- make_gs_matrix_named(
  t2g_bp,
  t2n_bp,
  prot_feats_liver_raw
)

gs_kegg_mat_liver <- make_gs_matrix_named(
  t2g_kegg,
  t2n_kegg,
  prot_feats_liver_raw
)

gs_rp_mat_liver <- make_gs_matrix_named(
  t2g_rp,
  t2n_rp,
  prot_feats_liver_raw
)

saveRDS(
  list(
    midgut = list(
      GO_BP = gs_bp_mat_midgut,
      KEGG = gs_kegg_mat_midgut,
      Reactome = gs_rp_mat_midgut
    ),
    liver = list(
      GO_BP = gs_bp_mat_liver,
      KEGG = gs_kegg_mat_liver,
      Reactome = gs_rp_mat_liver
    )
  ),
  file.path(rds_dir, "gene_set_matrices_for_MOFA_enrichment.rds")
)

# GO BP enrichment
res_bp_pos_midgut <- run_enrichment(
  model,
  view = "Midgut Proteome",
  factors = "all",
  p.adj.method = "BH",
  min.size = 10,
  verbose = TRUE,
  set.statistic = "rank.sum",
  feature.sets = gs_bp_mat_midgut,
  sign = "positive",
  statistical.test = "parametric"
)

res_bp_neg_midgut <- run_enrichment(
  model,
  view = "Midgut Proteome",
  factors = "all",
  p.adj.method = "BH",
  min.size = 10,
  verbose = TRUE,
  set.statistic = "rank.sum",
  feature.sets = gs_bp_mat_midgut,
  sign = "negative",
  statistical.test = "parametric"
)

res_bp_pos_liver <- run_enrichment(
  model,
  view = "Liver Proteome",
  factors = "all",
  p.adj.method = "BH",
  min.size = 10,
  verbose = TRUE,
  set.statistic = "rank.sum",
  feature.sets = gs_bp_mat_liver,
  sign = "positive",
  statistical.test = "parametric"
)

res_bp_neg_liver <- run_enrichment(
  model,
  view = "Liver Proteome",
  factors = "all",
  p.adj.method = "BH",
  min.size = 10,
  verbose = TRUE,
  set.statistic = "rank.sum",
  feature.sets = gs_bp_mat_liver,
  sign = "negative",
  statistical.test = "parametric"
)

gsea_results <- list(
  GO_BP = list(
    midgut_positive = res_bp_pos_midgut,
    midgut_negative = res_bp_neg_midgut,
    liver_positive = res_bp_pos_liver,
    liver_negative = res_bp_neg_liver
  )
)

saveRDS(
  gsea_results,
  file.path(rds_dir, "MOFA_GSEA_results_GO_BP_all_factors.rds")
)

make_combined_enrichment_plot <- function(p_pos, p_neg, plot_title) {
  df_pos <- p_pos$data |>
    transmute(
      pathway = pathway,
      logp = logp,
      pvalue = pvalues,
      sign = "+"
    )

  df_neg <- p_neg$data |>
    transmute(
      pathway = pathway,
      logp = logp,
      pvalue = pvalues,
      sign = "-"
    )

  df <- bind_rows(df_pos, df_neg) |>
    mutate(
      sign = factor(sign, levels = c("-", "+")),
      logp_signed = ifelse(sign == "-", -logp, logp),
      pathway_clean = gsub("^GO:[0-9]+: ", "", pathway),
      pathway_clean = stringr::str_trunc(pathway_clean, width = 60),
      pathway_clean = forcats::fct_reorder(pathway_clean, logp_signed)
    )

  ggplot(df, aes(x = logp_signed, y = pathway_clean, color = sign)) +
    geom_segment(
      aes(
        x = 0,
        xend = logp_signed,
        y = pathway_clean,
        yend = pathway_clean
      ),
      linewidth = 0.8,
      alpha = 0.8
    ) +
    geom_point(size = 5) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.7
    ) +
    scale_color_manual(
      values = c(
        "+" = "#8b0000",
        "-" = "#1f4e79"
      )
    ) +
    labs(
      x = expression("Signed " * -log[10] * "(adjp)"),
      y = NULL,
      color = "Sign",
      title = plot_title
    ) +
    theme_bw(base_size = 14) +
    theme(
      axis.text.y = element_text(size = 16),
      axis.text.x = element_text(size = 16),
      axis.title.x = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 16),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20),
      legend.key.size = unit(1.2, "cm"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "transparent", colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.box.background = element_rect(fill = "transparent", colour = NA),
      legend.key = element_rect(fill = "transparent", colour = NA)
    )
}
# FSEA detailed enrichment plots
protein_info_plot <- readr::read_tsv(
  protein_info_file,
  show_col_types = FALSE
) |>
  dplyr::rename(string_protein_id = `#string_protein_id`) |>
  dplyr::mutate(
    protein_id = stringr::str_remove(
      string_protein_id,
      "^STRG0A55HWH\\."
    ),
    annotation_clean = dplyr::if_else(
      is.na(annotation) | annotation == "",
      preferred_name,
      annotation
    ),
    annotation_clean = stringr::str_squish(annotation_clean),
    protein_label = paste(
      protein_id,
      annotation_clean,
      sep = " "
    )
  ) |>
  dplyr::select(
    protein_id,
    protein_label
  )

stopifnot(!any(stringr::str_detect(protein_info_plot$protein_label, "\n")))
stopifnot(!any(stringr::str_detect(protein_info_plot$protein_label, "___")))

make_clean_enrichment_lollipop <- function(p_enrich,
                                           protein_info,
                                           title = "",
                                           n_pathways = 5,
                                           n_genes = 5,
                                           y_text_size = 7,
                                           strip_text_size = 10,
                                           point_color = "#1f4e79") {
  df <- p_enrich$data |>
    mutate(
      feature_original = as.character(feature),
      protein_id = clean_protein_id(feature_original),
      protein_label = protein_info$protein_label[
        match(protein_id, protein_info$protein_id)
      ],
      protein_label = dplyr::if_else(
        is.na(protein_label),
        protein_id,
        protein_label
      ),
      protein_label = stringr::str_squish(protein_label),
      pathway_clean = pathway |>
        stringr::str_remove("^GO:[0-9]+: ") |>
        stringr::str_wrap(width = 45)
    )

  cat("Mapped proteins:", sum(!is.na(match(df$protein_id, protein_info$protein_id))), "\n")
  cat("Unmapped proteins:", sum(is.na(match(df$protein_id, protein_info$protein_id))), "\n")

  pathway_tbl <- df |>
    group_by(pathway, pathway_clean) |>
    summarise(
      pvalue_pathway = min(pvalue, na.rm = TRUE),
      n_total_genes = n_distinct(feature_original),
      .groups = "drop"
    ) |>
    arrange(pvalue_pathway) |>
    slice_head(n = n_pathways) |>
    mutate(
      pathway_label = paste0(
        pathway_clean,
        "\nN = ",
        n_total_genes,
        "; p = ",
        signif(pvalue_pathway, 2)
      ),
      pathway_rank = row_number()
    )

  df_plot <- df |>
    inner_join(
      pathway_tbl,
      by = c("pathway", "pathway_clean")
    ) |>
    group_by(pathway, pathway_label, pathway_rank) |>
    slice_max(
      order_by = abs(feature.statistic),
      n = n_genes,
      with_ties = FALSE
    ) |>
    ungroup() |>
    arrange(pathway_rank, feature.statistic) |>
    group_by(pathway_label) |>
    mutate(
      gene_rank = row_number()
    ) |>
    ungroup() |>
    mutate(
      pathway_label = factor(
        pathway_label,
        levels = rev(pathway_tbl$pathway_label)
      ),
      y_id = paste0("pathway_", pathway_rank, "_gene_", gene_rank),
      y_id = factor(
        y_id,
        levels = rev(unique(y_id))
      )
    )

  y_labels <- df_plot$protein_label
  names(y_labels) <- as.character(df_plot$y_id)

  stopifnot(!any(stringr::str_detect(y_labels, "___")))
  stopifnot(!any(stringr::str_detect(y_labels, "GO:")))
  stopifnot(!any(stringr::str_detect(y_labels, "N =")))

  p <- ggplot(
    df_plot,
    aes(
      x = feature.statistic,
      y = y_id
    )
  ) +
    geom_segment(
      aes(
        x = 0,
        xend = feature.statistic,
        y = y_id,
        yend = y_id
      ),
      linewidth = 0.6,
      color = "grey45"
    ) +
    geom_point(
      size = 2.8,
      color = point_color
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.5
    ) +
    facet_grid(
      pathway_label ~ .,
      scales = "free_y",
      space = "free_y"
    ) +
    scale_y_discrete(
      labels = y_labels
    ) +
    theme_bw(base_size = 14) +
    labs(
      title = title,
      x = "Weight",
      y = NULL
    ) +
    theme(
      plot.title = element_text(
        size = 16,
        face = "bold",
        hjust = 0.5
      ),
      strip.text.y = element_text(
        size = strip_text_size,
        face = "bold",
        angle = 0
      ),
      strip.background = element_rect(
        fill = "grey90",
        color = "black",
        linewidth = 0.4
      ),
      axis.text.y = element_text(size = y_text_size),
      axis.text.x = element_text(size = 11, angle = 90),
      axis.title.x = element_text(
        size = 13,
        face = "bold"
      ),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.y = grid::unit(0.6, "lines")
    )

  list(
    plot = p,
    data = df_plot
  )
}


# Selected-factor variance explained by view
selected_factors <- c(2, 3, 4, 11, 13)

selected_r2 <- as.data.frame(r2_dt) |>
  dplyr::mutate(
    factor_number = as.integer(as.character(factor))
  ) |>
  dplyr::filter(
    factor_number %in% selected_factors
  ) |>
  dplyr::mutate(
    factor = factor(
      paste0("Factor ", factor_number),
      levels = paste0("Factor ", selected_factors)
    ),
    view = factor(
      as.character(view),
      levels = rev(unique(as.character(r2_dt$view)))
    )
  )

write.csv(
  selected_r2,
  file.path(
    tab_dir,
    "variance_explained_by_view_selected_factors.csv"
  ),
  row.names = FALSE
)

saveRDS(
  selected_r2,
  file.path(
    rds_dir,
    "variance_explained_by_view_selected_factors.rds"
  )
)

p_selected_r2 <- ggplot(
  selected_r2,
  aes(x = r2, y = view)
) +
  geom_col(
    fill = "#00627d",
    colour = "black",
    linewidth = 0.4,
    width = 0.72
  ) +
  geom_text(
    aes(label = sprintf("%.2f", r2)),
    hjust = -0.15,
    size = 4
  ) +
  facet_wrap(
    ~ factor,
    ncol = 2,
    scales = "free_x"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.20)),
    breaks = scales::pretty_breaks(n = 4)
  ) +
  labs(
    x = "Variance explained (%)",
    y = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.y = element_text(
      size = 13,
      colour = "black"
    ),
    axis.text.x = element_text(
      size = 12,
      colour = "black"
    ),
    axis.title.x = element_text(
      size = 15,
      face = "bold",
      margin = margin(t = 8)
    ),
    strip.background = element_rect(
      fill = "grey90",
      colour = "black",
      linewidth = 0.4
    ),
    strip.text = element_text(
      size = 15,
      face = "bold",
      colour = "black",
      margin = margin(5, 5, 5, 5)
    ),
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5),
    panel.spacing = grid::unit(1, "lines"),
    plot.margin = margin(10, 15, 10, 10)
  )

print(p_selected_r2)

save_plot(
  p_selected_r2,
  filename = "variance_explained_by_view_selected_factors",
  width = 7,
  height = 7
)

# Shared factor-score theme
factor_day_theme <- theme_bw(base_size = 18) +
  theme(
    strip.background = element_rect(
      fill = "grey85",
      colour = "grey50",
      linewidth = 0.6
    ),
    strip.text = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(size = 14, angle = 90),
    axis.title.y = element_text(size = 20),
    axis.title.x = element_blank(),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.key.size = unit(0.5, "cm"),
    panel.grid = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

# Run the complete Factor 3 analysis layout for one factor
run_factor_analysis <- function(factor_id) {

  factor_label <- paste0("Factor", factor_id)
  factor_dir <- file.path(out_dir, paste0("Factor_", factor_id))
  factor_fig_dir <- file.path(factor_dir, "figures")
  factor_tab_dir <- file.path(factor_dir, "tables")
  factor_rds_dir <- file.path(factor_dir, "rds")

  dir.create(
    factor_fig_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  dir.create(
    factor_tab_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  dir.create(
    factor_rds_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  message("Running analysis for ", factor_label)

  # Factor scores grouped by Day and coloured by metadata
  score_covariates <- c(
    "Diet",
    "Plasma Glucose",
    "HSI",
    "Plasma L-Lactate",
    "Plasma D-Lactate",
    "Unmodified C",
    "5-mC",
    "5-hmC"
  )

  score_filenames <- c(
    "Diet" = "Diet",
    "Plasma Glucose" = "Plasma_Glucose",
    "HSI" = "HSI",
    "Plasma L-Lactate" = "Plasma_L_Lactate",
    "Plasma D-Lactate" = "Plasma_D_Lactate",
    "Unmodified C" = "Unmodified_C",
    "5-mC" = "5_mC",
    "5-hmC" = "5_hmC"
  )

  score_plots <- vector(
    mode = "list",
    length = length(score_covariates)
  )
  names(score_plots) <- score_covariates

  for (color_variable in score_covariates) {
    p_factor_score <- plot_factor(
      model,
      factors = factor_id,
      group_by = "Day",
      color_by = color_variable,
      add_dots = TRUE,
      dot_size = 4,
      dot_alpha = 0.8,
      add_violin = FALSE,
      add_boxplot = FALSE,
      show_missing = FALSE,
      scale = TRUE,
      dodge = FALSE,
      legend = TRUE,
      rasterize = FALSE
    ) +
      factor_day_theme

    if (identical(color_variable, "Diet")) {
      p_factor_score <- p_factor_score +
        scale_color_manual(values = cols2diet, name = "Diet") +
        scale_fill_manual(values = cols2diet, name = "Diet")
    }

    print(p_factor_score)

    save_plot_to(
      p_factor_score,
      filename = paste0(
        factor_label,
        "_scores_by_Day_",
        score_filenames[[color_variable]]
      ),
      width = if (identical(color_variable, "Plasma Glucose")) 8 else 7,
      height = 4,
      output_dir = factor_fig_dir
    )

    score_plots[[color_variable]] <- p_factor_score
  }

  # Factor data heatmaps
  save_heatmap_pdf_to(
    plot_data_heatmap(
      model,
      view = "Midgut Proteome",
      factor = factor_id,
      features = 20,
      show_colnames = FALSE,
      cluster_cols = TRUE,
      denoise = TRUE,
      scale = "row",
      min.value = -2,
      max.value = 2,
      color = heatmap_cols,
      breaks = seq(-2, 2, length.out = 101),
      annotation_samples = c(
        "Diet",
        "Day",
        "Temporal Phase",
        "Diet x Temporal Phase"
      ),
      annotation_colors = annotation_cols,
      annotation_legend = TRUE,
      silent = TRUE,
      border_color = NA
    ),
    filename = paste0(
      "heatmap_",
      factor_label,
      "_midgut_proteome"
    ),
    width = 10,
    height = 6,
    output_dir = factor_fig_dir
  )

  save_heatmap_pdf_to(
    plot_data_heatmap(
      model,
      view = "Liver Proteome",
      factor = factor_id,
      features = 20,
      show_colnames = FALSE,
      cluster_cols = TRUE,
      denoise = TRUE,
      scale = "row",
      min.value = -2,
      max.value = 2,
      color = heatmap_cols,
      breaks = seq(-2, 2, length.out = 101),
      annotation_samples = c(
        "Diet",
        "Day",
        "Temporal Phase",
        "Diet x Temporal Phase"
      ),
      annotation_colors = annotation_cols,
      annotation_legend = TRUE,
      silent = TRUE,
      border_color = NA
    ),
    filename = paste0(
      "heatmap_",
      factor_label,
      "_liver_proteome"
    ),
    width = 10,
    height = 6,
    output_dir = factor_fig_dir
  )

  save_heatmap_pdf_to(
    plot_data_heatmap(
      model,
      view = "Liver hPTMs",
      factor = factor_id,
      features = 20,
      show_colnames = FALSE,
      cluster_cols = TRUE,
      denoise = TRUE,
      scale = "row",
      min.value = -2,
      max.value = 2,
      color = heatmap_cols,
      breaks = seq(-2, 2, length.out = 101),
      annotation_samples = c(
        "Diet",
        "Day",
        "Temporal Phase",
        "Diet x Temporal Phase"
      ),
      annotation_colors = annotation_cols,
      annotation_legend = TRUE,
      silent = TRUE,
      border_color = NA
    ),
    filename = paste0(
      "heatmap_",
      factor_label,
      "_liver_hPTM"
    ),
    width = 10,
    height = 6,
    output_dir = factor_fig_dir
  )

  save_heatmap_pdf_to(
    plot_data_heatmap(
      model,
      view = "Midgut hPTMs",
      factor = factor_id,
      features = 10,
      show_colnames = FALSE,
      cluster_cols = TRUE,
      cluster_rows = FALSE,
      denoise = TRUE,
      scale = "row",
      min.value = -2,
      max.value = 2,
      color = heatmap_cols,
      breaks = seq(-2, 2, length.out = 101),
      annotation_samples = c(
        "Diet",
        "Day",
        "Temporal Phase",
        "Diet x Temporal Phase"
      ),
      annotation_colors = annotation_cols,
      annotation_legend = TRUE,
      silent = TRUE,
      border_color = NA
    ),
    filename = paste0(
      "heatmap_",
      factor_label,
      "_midgut_hPTM"
    ),
    width = 10,
    height = 5,
    output_dir = factor_fig_dir
  )

  save_heatmap_pdf_to(
    plot_data_heatmap(
      model,
      view = "Mucus Microbiota",
      factor = factor_id,
      features = 10,
      show_colnames = FALSE,
      cluster_cols = TRUE,
      denoise = TRUE,
      scale = "row",
      min.value = -2,
      max.value = 2,
      color = heatmap_cols,
      breaks = seq(-2, 2, length.out = 101),
      annotation_samples = c(
        "Diet",
        "Day",
        "Temporal Phase",
        "Diet x Temporal Phase"
      ),
      annotation_colors = annotation_cols,
      annotation_legend = TRUE,
      silent = TRUE,
      border_color = NA
    ),
    filename = paste0(
      "heatmap_",
      factor_label,
      "_mucus_microbiota"
    ),
    width = 10,
    height = 5,
    output_dir = factor_fig_dir
  )

  save_heatmap_pdf_to(
    plot_data_heatmap(
      model,
      view = "Digesta Microbiota",
      factor = factor_id,
      features = 10,
      show_colnames = FALSE,
      cluster_cols = TRUE,
      denoise = TRUE,
      scale = "row",
      min.value = -2,
      max.value = 2,
      color = heatmap_cols,
      breaks = seq(-2, 2, length.out = 101),
      annotation_samples = c(
        "Diet",
        "Day",
        "Temporal Phase",
        "Diet x Temporal Phase"
      ),
      annotation_colors = annotation_cols,
      annotation_legend = TRUE,
      silent = TRUE,
      border_color = NA
    ),
    filename = paste0(
      "heatmap_",
      factor_label,
      "_digesta_microbiota"
    ),
    width = 10,
    height = 5,
    output_dir = factor_fig_dir
  )

  # Top weights
  p_dag <- plot_top_weights(
    model,
    view = "Digesta Microbiota",
    factors = factor_id,
    scale = TRUE,
    nfeatures = 10
  )

  p_mam <- plot_top_weights(
    model,
    view = "Mucus Microbiota",
    factors = factor_id,
    scale = TRUE,
    nfeatures = 10
  )

  p_hptm_midgut <- plot_top_weights(
    model,
    view = "Midgut hPTMs",
    factors = factor_id,
    scale = TRUE,
    nfeatures = 10
  )

  p_hptm_liver <- plot_top_weights(
    model,
    view = "Liver hPTMs",
    factors = factor_id,
    scale = TRUE,
    nfeatures = 10
  )

  p_prot_midgut <- plot_top_weights(
    model,
    view = "Midgut Proteome",
    factors = factor_id,
    scale = TRUE,
    nfeatures = 10
  )

  p_prot_liver <- plot_top_weights(
    model,
    view = "Liver Proteome",
    factors = factor_id,
    scale = TRUE,
    nfeatures = 10
  )

  df_dag <- p_dag$data
  df_mam <- p_mam$data
  df_hptm_midgut <- p_hptm_midgut$data
  df_hptm_liver <- p_hptm_liver$data
  df_prot_midgut <- p_prot_midgut$data
  df_prot_liver <- p_prot_liver$data

  df_prot_midgut <- add_protein_labels(
    df_prot_midgut,
    protein_info
  )

  df_prot_liver <- add_protein_labels(
    df_prot_liver,
    protein_info
  )

  protein_label_qc <- bind_rows(
    df_prot_midgut |>
      mutate(view = "Midgut Proteome") |>
      select(view, feature, protein_id, protein_label),
    df_prot_liver |>
      mutate(view = "Liver Proteome") |>
      select(view, feature, protein_id, protein_label)
  )

  write.csv(
    protein_label_qc,
    file.path(
      factor_tab_dir,
      paste0(
        factor_label,
        "_top_weight_protein_label_QC.csv"
      )
    ),
    row.names = FALSE
  )

  p_dag_dot <- plot_weights_dot(
    df_dag,
    paste0("Digesta Microbiota - ", factor_label)
  ) +
    mono_bold_title

  p_mam_dot <- plot_weights_dot(
    df_mam,
    paste0("Mucus Microbiota - ", factor_label)
  ) +
    mono_bold_title

  p_hptm_midgut_dot <- plot_weights_dot(
    df_hptm_midgut,
    paste0("Midgut hPTMs - ", factor_label)
  ) +
    mono_bold_title

  p_hptm_liver_dot <- plot_weights_dot(
    df_hptm_liver,
    paste0("Liver hPTMs - ", factor_label)
  ) +
    mono_bold_title

  p_prot_midgut_dot <- plot_weights_dot(
    df_prot_midgut,
    paste0("Midgut Proteome - ", factor_label),
    use_feature_plot = TRUE
  ) +
    mono_bold_title

  p_prot_liver_dot <- plot_weights_dot(
    df_prot_liver,
    paste0("Liver Proteome - ", factor_label),
    use_feature_plot = TRUE
  ) +
    mono_bold_title

  print(p_dag_dot)
  print(p_mam_dot)
  print(p_hptm_midgut_dot)
  print(p_hptm_liver_dot)
  print(p_prot_midgut_dot)
  print(p_prot_liver_dot)

  top_weight_data <- list(
    digesta = df_dag,
    mucus = df_mam,
    midgut_hPTM = df_hptm_midgut,
    liver_hPTM = df_hptm_liver,
    midgut_proteome = df_prot_midgut,
    liver_proteome = df_prot_liver
  )

  saveRDS(
    top_weight_data,
    file.path(
      factor_rds_dir,
      paste0(factor_label, "_top_weight_data.rds")
    )
  )

  save_plot_to(
    p_dag_dot,
    filename = paste0(factor_label, "_digesta"),
    width = 8.5,
    height = 3.5,
    output_dir = factor_fig_dir
  )

  save_plot_to(
    p_mam_dot,
    filename = paste0(factor_label, "_mucus"),
    width = 8.5,
    height = 3.5,
    output_dir = factor_fig_dir
  )

  save_plot_to(
    p_hptm_midgut_dot,
    filename = paste0(factor_label, "_hptm_midgut"),
    width = 6.5,
    height = 3.5,
    output_dir = factor_fig_dir
  )

  save_plot_to(
    p_hptm_liver_dot,
    filename = paste0(factor_label, "_hptm_liver"),
    width = 6.5,
    height = 3.5,
    output_dir = factor_fig_dir
  )

  save_plot_to(
    p_prot_midgut_dot,
    filename = paste0(
      factor_label,
      "_proteome_midgut_annotated"
    ),
    width = 8.5,
    height = 3.5,
    output_dir = factor_fig_dir
  )

  save_plot_to(
    p_prot_liver_dot,
    filename = paste0(
      factor_label,
      "_proteome_liver_annotated"
    ),
    width = 8.5,
    height = 3.5,
    output_dir = factor_fig_dir
  )

  # GSEA plots
  p_bp_pos_midgut <- plot_enrichment(
    res_bp_pos_midgut,
    factor = factor_id,
    text_size = 1.5,
    dot_size = 5,
    alpha = 0.1,
    max.pathways = 10
  )

  p_bp_neg_midgut <- plot_enrichment(
    res_bp_neg_midgut,
    factor = factor_id,
    text_size = 1.5,
    dot_size = 5,
    alpha = 0.1,
    max.pathways = 10
  )

  p_bp_pos_liver <- plot_enrichment(
    res_bp_pos_liver,
    factor = factor_id,
    text_size = 1.5,
    dot_size = 5,
    alpha = 0.1,
    max.pathways = 10
  )

  p_bp_neg_liver <- plot_enrichment(
    res_bp_neg_liver,
    factor = factor_id,
    text_size = 1.5,
    dot_size = 5,
    alpha = 0.1,
    max.pathways = 10
  )

  saveRDS(
    list(
      midgut_positive = p_bp_pos_midgut$data,
      midgut_negative = p_bp_neg_midgut$data,
      liver_positive = p_bp_pos_liver$data,
      liver_negative = p_bp_neg_liver$data
    ),
    file.path(
      factor_rds_dir,
      paste0(
        "MOFA_GSEA_plot_data_GO_BP_",
        factor_label,
        ".rds"
      )
    )
  )

  p_bp_combined_midgut <- make_combined_enrichment_plot(
    p_bp_pos_midgut,
    p_bp_neg_midgut,
    paste0(
      "Midgut GSEA (GO BP) - ",
      factor_label
    )
  ) +
    mono_bold_title

  p_bp_combined_liver <- make_combined_enrichment_plot(
    p_bp_pos_liver,
    p_bp_neg_liver,
    paste0(
      "Liver GSEA (GO BP) - ",
      factor_label
    )
  ) +
    mono_bold_title

  print(p_bp_combined_midgut)
  print(p_bp_combined_liver)

  save_plot_to(
    p_bp_combined_midgut,
    filename = paste0(
      "midgut_enrichment_",
      factor_label,
      "_combined"
    ),
    width = 10,
    height = 5,
    output_dir = factor_fig_dir
  )

  save_plot_to(
    p_bp_combined_liver,
    filename = paste0(
      "liver_enrichment_",
      factor_label,
      "_combined"
    ),
    width = 9,
    height = 5,
    output_dir = factor_fig_dir
  )

  # FSEA detailed enrichment plots
  p_enrich_neg_midgut <- plot_enrichment_detailed(
    res_bp_neg_midgut,
    alpha = 0.1,
    factor = factor_id,
    max.genes = 5,
    max.pathways = 5,
    text_size = 4
  )

  p_enrich_neg_liver <- plot_enrichment_detailed(
    res_bp_neg_liver,
    alpha = 0.1,
    factor = factor_id,
    max.genes = 5,
    max.pathways = 5,
    text_size = 4
  )

  p_enrich_pos_midgut <- plot_enrichment_detailed(
    res_bp_pos_midgut,
    alpha = 0.1,
    factor = factor_id,
    max.genes = 5,
    max.pathways = 5,
    text_size = 4
  )

  p_enrich_pos_liver <- plot_enrichment_detailed(
    res_bp_pos_liver,
    alpha = 0.1,
    factor = factor_id,
    max.genes = 5,
    max.pathways = 5,
    text_size = 4
  )

  fsea_neg_midgut <- make_clean_enrichment_lollipop(
    p_enrich_neg_midgut,
    protein_info_plot,
    title = paste0(
      "Midgut FSEA - ",
      factor_label,
      " (-)"
    ),
    n_pathways = 5,
    n_genes = 5,
    y_text_size = 12,
    strip_text_size = 12,
    point_color = "#1f4e79"
  )

  fsea_neg_liver <- make_clean_enrichment_lollipop(
    p_enrich_neg_liver,
    protein_info_plot,
    title = paste0(
      "Liver FSEA - ",
      factor_label,
      " (-)"
    ),
    n_pathways = 5,
    n_genes = 5,
    y_text_size = 12,
    strip_text_size = 12,
    point_color = "#1f4e79"
  )

  fsea_pos_midgut <- make_clean_enrichment_lollipop(
    p_enrich_pos_midgut,
    protein_info_plot,
    title = paste0(
      "Midgut FSEA - ",
      factor_label,
      " (+)"
    ),
    n_pathways = 5,
    n_genes = 5,
    y_text_size = 12,
    strip_text_size = 12,
    point_color = "#8b0000"
  )

  fsea_pos_liver <- make_clean_enrichment_lollipop(
    p_enrich_pos_liver,
    protein_info_plot,
    title = paste0(
      "Liver FSEA - ",
      factor_label,
      " (+)"
    ),
    n_pathways = 5,
    n_genes = 5,
    y_text_size = 12,
    strip_text_size = 12,
    point_color = "#8b0000"
  )

  p_clean_neg_midgut <- fsea_neg_midgut$plot + mono_bold_title
  p_clean_neg_liver <- fsea_neg_liver$plot + mono_bold_title
  p_clean_pos_midgut <- fsea_pos_midgut$plot + mono_bold_title
  p_clean_pos_liver <- fsea_pos_liver$plot + mono_bold_title

  print(p_clean_neg_midgut)
  print(p_clean_neg_liver)
  print(p_clean_pos_midgut)
  print(p_clean_pos_liver)

  fsea_results <- list(
    negative_midgut = fsea_neg_midgut$data,
    negative_liver = fsea_neg_liver$data,
    positive_midgut = fsea_pos_midgut$data,
    positive_liver = fsea_pos_liver$data
  )

  saveRDS(
    list(
      raw_gsea_results = gsea_results,
      detailed_plot_objects = list(
        negative_midgut = p_enrich_neg_midgut,
        negative_liver = p_enrich_neg_liver,
        positive_midgut = p_enrich_pos_midgut,
        positive_liver = p_enrich_pos_liver
      ),
      cleaned_fsea_plot_data = fsea_results
    ),
    file.path(
      factor_rds_dir,
      paste0(
        "MOFA_FSEA_results_GO_BP_",
        factor_label,
        ".rds"
      )
    )
  )

  write.csv(
    fsea_neg_midgut$data,
    file.path(
      factor_tab_dir,
      paste0(
        factor_label,
        "_negative_FSEA_midgut_plot_data.csv"
      )
    ),
    row.names = FALSE
  )

  write.csv(
    fsea_neg_liver$data,
    file.path(
      factor_tab_dir,
      paste0(
        factor_label,
        "_negative_FSEA_liver_plot_data.csv"
      )
    ),
    row.names = FALSE
  )

  write.csv(
    fsea_pos_midgut$data,
    file.path(
      factor_tab_dir,
      paste0(
        factor_label,
        "_positive_FSEA_midgut_plot_data.csv"
      )
    ),
    row.names = FALSE
  )

  write.csv(
    fsea_pos_liver$data,
    file.path(
      factor_tab_dir,
      paste0(
        factor_label,
        "_positive_FSEA_liver_plot_data.csv"
      )
    ),
    row.names = FALSE
  )

  save_plot_to(
    p_clean_neg_midgut,
    filename = paste0(
      factor_label,
      "_negative_FSEA_midgut"
    ),
    width = 10,
    height = 6,
    output_dir = factor_fig_dir
  )

  save_plot_to(
    p_clean_neg_liver,
    filename = paste0(
      factor_label,
      "_negative_FSEA_liver"
    ),
    width = 10,
    height = 6,
    output_dir = factor_fig_dir
  )

  save_plot_to(
    p_clean_pos_midgut,
    filename = paste0(
      factor_label,
      "_positive_FSEA_midgut"
    ),
    width = 10,
    height = 6,
    output_dir = factor_fig_dir
  )

  save_plot_to(
    p_clean_pos_liver,
    filename = paste0(
      factor_label,
      "_positive_FSEA_liver"
    ),
    width = 12,
    height = 6,
    output_dir = factor_fig_dir
  )

  saveRDS(
    score_plots,
    file.path(
      factor_rds_dir,
      paste0(factor_label, "_score_plots.rds")
    )
  )

  message("Completed analysis for ", factor_label)

  invisible(
    list(
      factor = factor_id,
      directory = factor_dir,
      figures = factor_fig_dir,
      tables = factor_tab_dir,
      rds = factor_rds_dir
    )
  )
}

# Run selected factors
factor_results <- vector(
  mode = "list",
  length = length(selected_factors)
)

names(factor_results) <- paste0(
  "Factor",
  selected_factors
)

for (i in seq_along(selected_factors)) {
  factor_results[[i]] <- run_factor_analysis(
    factor_id = selected_factors[[i]]
  )
}

saveRDS(
  factor_results,
  file.path(
    rds_dir,
    "selected_factor_output_directories.rds"
  )
)


# Merge selected factor-score plots coloured by Diet

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(purrr)
})

# Parameters to adjust
selected_factors <- c(2, 3, 4, 11, 13)

merged_ncol <- 2
merged_width <- 12
merged_height <- 10

out_dir <- "mofa_results"
fig_dir <- file.path(out_dir, "figures")
rds_dir <- file.path(out_dir, "rds")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

# Locate the score-plot RDS generated inside the factor loop
score_rds_files <- file.path(
  out_dir,
  paste0("Factor_", selected_factors),
  "rds",
  paste0("Factor", selected_factors, "_score_plots.rds")
)

# Check that all files exist
missing_files <- score_rds_files[!file.exists(score_rds_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing score-plot RDS files:\n",
    paste(missing_files, collapse = "\n")
  )
}

# Extract the Diet plot from each factor RDS
diet_factor_plots <- purrr::map2(
  score_rds_files,
  selected_factors,
  function(rds_file, factor_id) {
    
    score_plots <- readRDS(rds_file)
    
    if (!"Diet" %in% names(score_plots)) {
      stop(
        "Diet plot not found in: ",
        rds_file
      )
    }
    
    score_plots[["Diet"]] +
      theme(
        strip.text = element_text(
          size = 12,
          face = "bold",
          colour = "black"
        ),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        
        axis.text.y = element_text(
          size = 11,
          colour = "black"
        ),
        axis.title.y = element_text(
          size = 15,
          face = "bold"
        ),
        legend.title = element_text(
          size = 14,
          face = "bold"
        ),
        legend.text = element_text(size = 13),
        legend.key.size = grid::unit(0.5, "cm"),
        plot.margin = margin(5, 5, 5, 5)
      )
  }
)

names(diet_factor_plots) <- paste0(
  "Factor",
  selected_factors
)

# Merge plots and collect one common legend
p_selected_factors_diet <- (
  patchwork::wrap_plots(
    diet_factor_plots,
    ncol = merged_ncol
  ) +
    patchwork::plot_layout(
      guides = "collect"
    )
) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal"
  )

# Print in console
print(p_selected_factors_diet)


ggsave(
  filename = file.path(
    fig_dir,
    "selected_factors_scores_by_Day_Diet_merged.pdf"
  ),
  plot = p_selected_factors_diet,
  width = merged_width,
  height = merged_height,
  units = "in",
  device = "pdf",
  bg = "transparent",
  useDingbats = FALSE
)

ggsave(
  filename = file.path(
    fig_dir,
    "selected_factors_scores_by_Day_Diet_merged.tiff"
  ),
  plot = p_selected_factors_diet,
  width = merged_width,
  height = merged_height,
  units = "in",
  dpi = 1200,
  compression = "lzw",
  bg = "transparent"
)

saveRDS(
  p_selected_factors_diet,
  file.path(
    rds_dir,
    "selected_factors_scores_by_Day_Diet_merged.rds"
  )
)




# Session information
writeLines(
  capture.output(sessionInfo()),
  con = file.path(out_dir, "sessionInfo.txt")
)

message("MOFA downstream analysis complete.")
message("Results directory: ", out_dir)
