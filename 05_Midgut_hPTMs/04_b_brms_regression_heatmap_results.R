library(dplyr)
library(tidyr)
library(ggtree)
library(ggplot2)
library(cowplot)
library(ggnewscale)
library(stringr)
library(tibble)

# Directories
base_dir    <- "04_brms_regression/brms_best_family_final_batch_fixed_20260628_162843"
heatmap_dir <- file.path(base_dir, "summary_tables")

if (!dir.exists(heatmap_dir)) dir.create(heatmap_dir, recursive = TRUE)

# ------------------------------------------------------------
# Load contrasts from new best-family pipeline
# ------------------------------------------------------------
contrast_df <- read.csv(
  file.path(heatmap_dir, "all_features_contrast_diet_by_day_best_family_posterior.csv"),
  stringsAsFactors = FALSE
)

normalize_ptm <- function(x) {
  vapply(x, function(entry) {
    parts <- strsplit(entry, "\\|")[[1]]
    norm_parts <- vapply(parts, function(p) {
      g <- regmatches(p, regexec("^(H[34])(K|R)([0-9]{1,3})([A-Za-z0-9]+)$", p))[[1]]
      if (length(g) < 5 || is.na(g[1])) return(NA_character_)
      paste0(g[2], g[3], g[4], tolower(g[5]))
    }, character(1))
    if (any(is.na(norm_parts))) return(NA_character_)
    paste(norm_parts, collapse = "|")
  }, character(1), USE.NAMES = FALSE)
}

contrast_df$feature <- normalize_ptm(as.character(contrast_df$feature))

print(unique(contrast_df$feature))   # <-- check RIGHT HERE, same script run

contrast_df <- contrast_df %>%
  mutate(
    cri_excl_zero = (lower_95 > 0 | upper_95 < 0),
    Pd = ifelse(prob_gt0 < 0.5, 1 - prob_gt0, prob_gt0)
  )

# ------------------------------------------------------------------
# Filter to HC - NC and prep for clustering
# ------------------------------------------------------------------

heat_df <- contrast_df %>%
  filter(contrast == "HC - NC") %>%
  mutate(
    day_categ = as.character(day_categ),
    feature   = as.character(feature)
  )


TISSUE <- "Midgut"   # set explicitly here; use "Liver" when you run the liver hPTM version
LAYER  <- "hPTM"

set.seed(123)
boot_median_ci <- function(x, n_boot = 2000, conf = 0.95) {
  x <- x[!is.na(x)]
  if (length(x) < 2) {
    return(c(median = if (length(x) == 1) x else NA_real_,
             lower = NA_real_, upper = NA_real_))
  }
  boot_meds <- replicate(n_boot, median(sample(x, length(x), replace = TRUE)))
  alpha <- (1 - conf) / 2
  c(
    median = median(x),
    lower  = unname(quantile(boot_meds, alpha, na.rm = TRUE)),
    upper  = unname(quantile(boot_meds, 1 - alpha, na.rm = TRUE))
  )
}

activation_summary <- heat_df %>%
  mutate(day_num = as.numeric(gsub("day_", "", day_categ))) %>%
  group_by(day_num) %>%
  summarise(
    boot     = list(boot_median_ci(abs(estimate))),
    n_tested = n(),
    .groups  = "drop"
  ) %>%
  mutate(
    median_abs_delta = purrr::map_dbl(boot, "median"),
    ci_lower         = purrr::map_dbl(boot, "lower"),
    ci_upper         = purrr::map_dbl(boot, "upper")
  ) %>%
  select(-boot) %>%
  rename(day = day_num) %>%
  right_join(tibble(day = c(1, 2, 3, 4, 10, 15, 22)), by = "day") %>%
  mutate(
    tissue = TISSUE,
    layer  = LAYER
  ) %>%
  arrange(day) %>%
  mutate(row_id = paste(tissue, layer, day, sep = "_")) %>%
  column_to_rownames("row_id") %>%
  select(tissue, layer, day, median_abs_delta, ci_lower, ci_upper, n_tested)

saveRDS(
  activation_summary,
  file.path(heatmap_dir, paste0(tolower(LAYER), "_", tolower(TISSUE), "_median_abs_delta.rds"))
)

write.csv(
  activation_summary,
  file.path(heatmap_dir, paste0(tolower(LAYER), "_", tolower(TISSUE), "_median_abs_delta.csv")),
  row.names = TRUE
)

combined <- activation_summary %>%
  mutate(
    day_lab      = factor(paste0("Day ", day),
                          levels = paste0("Day ", c(1, 2, 3, 4, 10, 15, 22))),
    day_idx      = as.numeric(day_lab),
    tissue_layer = factor(paste(tissue, layer, sep = " - "))
  )

ggplot(combined, aes(x = day_idx, y = median_abs_delta, group = tissue_layer)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "#44a7c4", alpha = 0.4) +
  geom_line(linewidth = 1, color = "#2171b5") +
  geom_point(aes(size = n_tested), color = "#2171b5") +
  scale_x_continuous(breaks = 1:7, labels = levels(combined$day_lab)) +
  scale_size_continuous(range = c(2, 8), name = "n features tested") +
  theme_grey(base_size = 13) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  labs(x = NULL, y = "Median |estimate| (bootstrap 95% CI)",
       title = paste(unique(combined$layer), "-", unique(combined$tissue), "effect size dynamics"))


# For clustering: numeric day + wide matrix
heat_df2 <- heat_df %>%
  mutate(
    day_num   = as.numeric(gsub("day_", "", day_categ)),
    day_label = paste0("Day ", day_num)
  )

clust_df <- heat_df2 %>%
  dplyr::select(feature, day_num, estimate) %>%
  dplyr::arrange(day_num) %>%
  tidyr::pivot_wider(names_from = day_num, values_from = estimate)

row_names  <- clust_df$feature
clust_mat  <- as.matrix(clust_df[, -1, drop = FALSE])
rownames(clust_mat) <- row_names

# Hierarchical clustering on PTMs
hc        <- hclust(dist(clust_mat), method = "complete")
ptm_order <- rownames(clust_mat)[hc$order]

ddgram      <- as.dendrogram(hc)
ggtree_plot <- ggtree::ggtree(ddgram, branch.length = "none")

# ------------------------------------------------------------
# Build plotting data with Pd and HDI info
# ------------------------------------------------------------

heat_df2 <- heat_df %>%
  mutate(
    day_num   = as.numeric(gsub("day_", "", day_categ)),
    day_label = paste0("Day ", day_num),
    feature   = as.character(feature)
  ) %>%
  filter(feature %in% ptm_order)

# enforce PTM order from clustering
heat_df2$feature <- factor(heat_df2$feature, levels = ptm_order)

# nice day ordering
heat_df2$day_label <- factor(
  heat_df2$day_label,
  levels = paste0("Day ", sort(unique(heat_df2$day_num)))
)

# add direction, |effect|, and Pd label
heat_df2 <- heat_df2 %>%
  mutate(
    direction  = ifelse(estimate >= 0, "HC > NC", "HC < NC"),
    abs_est    = abs(estimate),
    prob_label = ifelse(Pd > 0.95, sprintf("%.2f", Pd), "")
  )

cols2 <- c(
  "HC < NC" = "#44a7c4",   # blue
  "HC > NC" = "#f99943"    # red
)

# ------------------------------------------------------------
# Dot-plot: fill = direction, size = |effect|
# Outline: green & thicker when HDI excludes zero
# ------------------------------------------------------------

p_dot <- ggplot(heat_df2, aes(x = day_label, y = feature)) +
  
  # base dot: fill and radius from |effect|
  geom_point(
    aes(size = abs_est, fill = direction),
    shape  = 21,
    colour = "black",
    alpha  = 0.9,
    stroke = 1
  ) +
  
  scale_fill_manual(
    name   = "Direction",
    values = cols2,
    labels = c("HC < NC", "HC > NC"),
    guide  = guide_legend(
      override.aes = list(size = 10),
      order = 1
    )
  ) +
  
  scale_size_continuous(
    name  = "|HC − NC| (log2 FC)",
    range = c(4, 20),
    guide = guide_legend(order = 2)
  ) +
  
  ggnewscale::new_scale("colour") +
  
  # overlay outline circle: HDI info
  geom_point(
    aes(
      size   = abs_est,
      colour = cri_excl_zero
    ),
    shape  = 21,
    fill   = NA,
    alpha  = 1,
    stroke = ifelse(heat_df2$cri_excl_zero, 2.2, 1)
  ) +
  
  scale_colour_manual(
    name   = "95% CrI",
    values = c(
      `FALSE` = "black",
      `TRUE`  = "green4"
    ),
    labels = c(
      `FALSE` = "includes 0",
      `TRUE`  = "excludes 0"
    ),
    guide = guide_legend(
      override.aes = list(size = 10, shape = 21, fill = NA),
      order = 4
    )
  ) +
  
  # Pd text on top of the dot
  geom_text(
    aes(label = prob_label),
    size   = 4.5,
    colour = "black",
    vjust  = 0.5
  ) +
  
  # Pd legend explanation
  geom_point(
    aes(shape = "Pd text"),
    x = NA, y = NA, size = 0
  ) +
  
  scale_shape_manual(
    name   = "Text labels",
    values = c("Pd text" = NA),
    labels = c("Pd (probability of direction)\nshown when Pd > 0.95"),
    guide  = guide_legend(order = 3)
  ) +
  
  labs(
    #title = "Midgut",
    x = NULL,
    y = NULL
  ) +
  
  theme_bw(base_size = 15) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size  = 18
    ),
    axis.text.y      = element_text(size = 18),
    plot.title = element_text(
      size = 22,
      face = "bold",
      family = "mono",
      hjust = 0.5
    ),
    legend.title     = element_text(size = 14, face = "bold"),
    legend.text      = element_text(size = 14),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white"),
    panel.background = element_rect(fill = "grey90")
  )

plot_grid <- plot_grid(
  ggtree_plot, p_dot, nrow = 1,
  rel_widths = c(0.25, 2.25), align = "h"
)

plot_grid

ggsave(
  filename = file.path(heatmap_dir, "dotplot_diet_effect_best_family.tiff"),
  plot     = plot_grid,
  device   = "tiff",
  width    = 10,
  height   = 12,
  dpi      = 1200
)

ggsave(
  filename = file.path(heatmap_dir, "dotplot_diet_effect_best_family.pdf"),
  plot     = plot_grid,
  device   = "pdf",
  width    = 10,
  height   = 12,
  dpi      = 1200
)
