suppressPackageStartupMessages({
  library(MOFA2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

set.seed(42)

# Input and output paths
model_file <- "mofa_models/K20_setseed_123.hdf5"

out_dir <- "mofa_results"
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
rds_dir <- file.path(out_dir, "rds")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(model_file))

# Save PDF and TIFF versions of a ggplot or patchwork object
save_plot <- function(plot, filename, width, height, dpi = 1200) {
  
  transparent_theme <- theme(
    plot.background = element_rect(fill = NA, colour = NA),
    panel.background = element_rect(fill = NA, colour = NA),
    legend.background = element_rect(fill = NA, colour = NA),
    legend.box.background = element_rect(fill = NA, colour = NA),
    legend.key = element_rect(fill = NA, colour = NA)
  )
  
  # Apply to every component of a patchwork object
  if (inherits(plot, "patchwork")) {
    plot <- plot & transparent_theme
  } else {
    plot <- plot + transparent_theme
  }
  
  ggsave(
    filename = file.path(fig_dir, paste0(filename, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = "pdf",
    bg = "transparent",
    useDingbats = FALSE
  )
  
  ggsave(
    filename = file.path(fig_dir, paste0(filename, ".tiff")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    compression = "lzw",
    bg = "transparent"
  )
}
# Load model
model <- load_model(
  model_file,
  load_interpol_Z = FALSE,
  sort_factors = TRUE,
  verbose = TRUE
)


# Define the top-to-bottom order for all figure panels
view_order <- c(
  "Digesta Microbiota",
  "Mucus Microbiota",
  "Midgut hPTMs",
  "Midgut Proteome",
  "Liver hPTMs",
  "Liver Proteome"
)

view_order <- view_order[view_order %in% views_names(model)]

if (length(view_order) == 0) {
  stop("None of the expected view names were found in the model.")
}

# Extract data in long format
data_long <- get_data(model, as.data.frame = TRUE) |>
  mutate(
    view = as.character(view),
    sample = as.character(sample),
    feature = as.character(feature)
  ) |>
  filter(view %in% view_order)

# Get model sample order
sample_order <- samples_names(model)

if (is.list(sample_order)) {
  sample_order <- sample_order[[1]]
}

sample_order <- as.character(sample_order)
n_total_samples <- length(sample_order)

# Calculate the number of samples and features for each view
view_info <- data_long |>
  group_by(view) |>
  summarise(
    D = n_distinct(feature),
    N_view = n_distinct(sample[!is.na(value)]),
    .groups = "drop"
  ) |>
  mutate(
    view = factor(view, levels = view_order)
  ) |>
  arrange(view) |>
  mutate(
    view_label_text = paste0(
      as.character(view),
      "\nN=", N_view,
      ", D=", D
    ),
    view_label = factor(
      view_label_text,
      levels = rev(view_label_text)
    )
  )

# Build presence / absence matrix across samples and views
present_df <- data_long |>
  group_by(view, sample) |>
  summarise(
    present = any(!is.na(value)),
    .groups = "drop"
  )

data_overview_df <- tidyr::expand_grid(
  view = view_order,
  sample = sample_order
) |>
  left_join(
    present_df,
    by = c("view", "sample")
  ) |>
  mutate(
    present = replace_na(present, FALSE),
    sample_index = match(sample, sample_order)
  ) |>
  left_join(
    view_info |>
      mutate(view = as.character(view)) |>
      select(view, view_label, N_view, D),
    by = "view"
  ) |>
  mutate(
    status = ifelse(present, view, "Missing"),
    view_label = factor(
      view_label,
      levels = levels(view_info$view_label)
    )
  )

# View colours
view_cols <- c(
  "Digesta Microbiota" = "#2ca25f",
  "Mucus Microbiota" = "#7570b3",
  "Midgut hPTMs" = "#fb8072",
  "Midgut Proteome" = "#1f78b4",
  "Liver hPTMs" = "#e66101",
  "Liver Proteome" = "#ffd92f",
  "Missing" = "grey80"
)

# Panel A: data availability overview
pA <- ggplot(
  data_overview_df,
  aes(x = sample_index, y = view_label, fill = status)
) +
  geom_tile(width = 1, height = 0.95) +
  scale_fill_manual(values = view_cols, guide = "none") +
  annotate(
    "text",
    x = n_total_samples / 2,
    y = length(view_order) + 0.65,
    label = paste0("N=", n_total_samples),
    fontface = "bold",
    size = 5.5
  ) +
  coord_cartesian(clip = "off") +
  labs(
    tag = "A",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.text.y = element_text(size = 15, colour = "black"),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    plot.tag = element_text(size = 28, face = "bold"),
    plot.margin = margin(25, 5, 5, 5)
  )

# Panel B: total variance explained per view
r2_total <- model@cache$variance_explained$r2_total[[1]]

r2_total_df <- enframe(
  r2_total,
  name = "view",
  value = "r2_total"
) |>
  filter(view %in% view_order) |>
  left_join(
    view_info |>
      mutate(view = as.character(view)) |>
      select(view, view_label),
    by = "view"
  ) |>
  mutate(
    view_label = factor(
      view_label,
      levels = levels(view_info$view_label)
    )
  )

pB <- ggplot(
  r2_total_df,
  aes(x = r2_total, y = view_label)
) +
  geom_col(
    fill = "#00627d",
    colour = "black",
    linewidth = 0.35,
    width = 0.75
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    tag = "B",
    x = expression(R^2~"(%)"),
    y = NULL
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.text.x = element_text(size = 14, colour = "black"),
    axis.title.x = element_text(size = 16, face = "bold"),
    plot.tag = element_text(size = 28, face = "bold"),
    plot.margin = margin(25, 5, 5, 5)
  )

# Panel C: factor-wise variance explained per view
r2 <- model@cache$variance_explained$r2_per_factor[[1]]

factor_names <- rownames(r2)

if (is.null(factor_names)) {
  factor_names <- paste0("Factor", seq_len(nrow(r2)))
  rownames(r2) <- factor_names
}

factors_to_plot <- factor_names[
  seq_len(min(20, length(factor_names)))
]

r2_heatmap_df <- as.data.frame(
  r2[factors_to_plot, view_order, drop = FALSE]
) |>
  rownames_to_column("factor") |>
  pivot_longer(
    cols = -factor,
    names_to = "view",
    values_to = "r2"
  ) |>
  left_join(
    view_info |>
      mutate(view = as.character(view)) |>
      select(view, view_label),
    by = "view"
  ) |>
  mutate(
    factor_label = str_remove(factor, "^Factor"),
    factor_label = factor(
      factor_label,
      levels = str_remove(factors_to_plot, "^Factor")
    ),
    view_label = factor(
      view_label,
      levels = levels(view_info$view_label)
    )
  )

pC <- ggplot(
  r2_heatmap_df,
  aes(x = factor_label, y = view_label, fill = r2)
) +
  geom_tile(colour = "grey60", linewidth = 0.35) +
  scale_fill_gradient(
    low = "grey95",
    high = "navy",
    limits = c(0, 3),
    oob = scales::squish,
    name = expression(R^2~"(%)")
  ) +
  labs(
    tag = "C",
    x = "Factor",
    y = NULL
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.text.x = element_text(size = 14, colour = "black"),
    axis.title.x = element_text(size = 16, face = "bold"),
    plot.tag = element_text(size = 28, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    plot.margin = margin(25, 5, 5, 5)
  )

# Combine panels
fig_mofa_ABC <- pA + pB + pC +
  plot_layout(widths = c(1.35, 0.9, 1.75))

pA_notag <- pA + labs(tag = NULL)
pB_notag <- pB + labs(tag = NULL)
pC_notag <- pC + labs(tag = NULL)

fig_mofa_ABC_notag <- patchwork::wrap_plots(
  pA_notag,
  pB_notag,
  pC_notag,
  nrow = 1,
  widths = c(1.25, 0.9, 1.85)
)

print(fig_mofa_ABC)
print(fig_mofa_ABC_notag)

# Save combined figures
save_plot(
  fig_mofa_ABC,
  filename = "Figure5_MOFA_overview_aligned_ABC_with_tags",
  width = 16,
  height = 6
)

save_plot(
  fig_mofa_ABC_notag,
  filename = "Figure5_MOFA_overview_aligned_ABC_no_tags",
  width = 16,
  height = 6
)


# Export data used for the figure
write.csv(
  view_info |>
    mutate(view = as.character(view)) |>
    select(view, N_view, D),
  file.path(tab_dir, "MOFA_view_N_and_D_summary.csv"),
  row.names = FALSE
)

write.csv(
  r2_total_df |>
    select(view, r2_total),
  file.path(tab_dir, "MOFA_total_R2_per_view.csv"),
  row.names = FALSE
)

write.csv(
  r2_heatmap_df |>
    transmute(
      view,
      factor = as.character(factor_label),
      r2
    ),
  file.path(tab_dir, "MOFA_factor_R2_per_view.csv"),
  row.names = FALSE
)

# Save reusable figure objects and source data
saveRDS(
  list(
    view_info = view_info,
    data_overview = data_overview_df,
    r2_total = r2_total_df,
    r2_factor_by_view = r2_heatmap_df,
    plots = list(
      panel_A = pA,
      panel_B = pB,
      panel_C = pC,
      figure_ABC = fig_mofa_ABC,
      figure_ABC_no_tags = fig_mofa_ABC_notag
    )
  ),
  file.path(rds_dir, "MOFA_overview_figure_objects.rds")
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(out_dir, "sessionInfo_MOFA_overview.txt")
)

message("MOFA overview results saved in: ", out_dir)
