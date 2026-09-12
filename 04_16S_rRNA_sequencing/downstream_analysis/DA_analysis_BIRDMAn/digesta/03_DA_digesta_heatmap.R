suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggtree)
  library(ggplot2)
  library(cowplot)
  library(ggnewscale)
  library(readr)
  library(stringr)
})

# OUTPUT DIRECTORY

PLOT_OUTDIR <- "DA_digesta_heatmap"
dir.create(PLOT_OUTDIR, recursive = TRUE, showWarnings = FALSE)


# INPUTS

f_contr_long <- "export_for_R/HCminusNC_contrasts_posterior_summary_LONG.tsv"
f_diag_long  <- file.path("export_for_R", "diagnostics", "delta_contrast_diagnostics_LONG.tsv")

day_levels <- c("day_01","day_02","day_03","day_04","day_10","day_15","day_22")

hdi_prob <- "0.95"
hdi_low_col  <- paste0("hdi_lower_", hdi_prob)
hdi_high_col <- paste0("hdi_upper_", hdi_prob)

pd_keep_thresh  <- 0.95
pd_label_thresh <- 0.95

# diagnostics thresholds (ASV-LEVEL!)
rhat_bad     <- 1.05
ess_tail_bad <- 400


# HELPERS

stop_if_missing <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0) stop("Missing file(s):\n- ", paste(missing, collapse = "\n- "))
}
assert_cols <- function(df, cols, df_name = "data") {
  miss <- setdiff(cols, colnames(df))
  if (length(miss) > 0) stop(df_name, " missing columns: ", paste(miss, collapse = ", "))
}
day_to_label <- function(x) {
  n <- suppressWarnings(as.integer(str_extract(x, "\\d+")))
  paste0("Day ", n)
}

# genus extraction: expects "g__Something" somewhere in feature string
extract_genus <- function(feature) {
  g <- stringr::str_match(feature, "g__([^;|\\s]+)")[,2]
  ifelse(is.na(g), NA_character_, g)
}
bad_genus_values <- c("","NA","N/A","none","None","unclassified","Unclassified","uncultured","Unknown")


# CHECKS + LOAD

stop_if_missing(c(f_contr_long, f_diag_long))

contrast_df <- read_tsv(f_contr_long, show_col_types = FALSE)
diag_df     <- read_tsv(f_diag_long,  show_col_types = FALSE)

assert_cols(contrast_df, c("day","feature","asv_id","mean","PD", hdi_low_col, hdi_high_col), "contrasts LONG")
assert_cols(diag_df,     c("day","asv_id","r_hat","ess_tail"), "diagnostics LONG")

contrast_df <- contrast_df %>%
  mutate(
    day     = as.character(day),
    asv_id  = str_trim(as.character(asv_id)),
    feature = str_trim(as.character(feature))
  )

diag_df <- diag_df %>%
  mutate(
    day    = as.character(day),
    asv_id = str_trim(as.character(asv_id))
  )


# DAY ORDER

obs_days <- sort(unique(contrast_df$day))
keep_days <- intersect(day_levels, obs_days)
if (length(keep_days) == 0) keep_days <- obs_days

day_labels <- day_to_label(keep_days)


# ASV-LEVEL DIAGNOSTIC FILTER

diag_asv <- diag_df %>%
  filter(day %in% keep_days) %>%
  group_by(asv_id) %>%
  summarise(
    worst_rhat   = max(r_hat, na.rm = TRUE),
    min_ess_tail = min(ess_tail, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    keep_diag = (worst_rhat < rhat_bad) & (min_ess_tail > ess_tail_bad)
  )

good_asv <- diag_asv %>% filter(keep_diag) %>% pull(asv_id)
bad_asv  <- diag_asv %>% filter(!keep_diag) %>% pull(asv_id)

cat("Removing ASVs due to diagnostics:", length(bad_asv), "\n")

contrast_df <- contrast_df %>%
  filter(asv_id %in% good_asv)


# GENUS FILTER (ASV-LEVEL VIA FEATURE STRING)

contrast_df <- contrast_df %>%
  mutate(genus = extract_genus(feature)) %>%
  filter(
    !is.na(genus),
    !(genus %in% bad_genus_values)
  ) %>%
  select(-genus)


# KEEP DAYS + JOIN DIAGNOSTICS (for completeness; not used for filtering anymore)

contrast_df <- contrast_df %>%
  filter(day %in% keep_days) %>%
  mutate(day = factor(day, levels = keep_days)) %>%
  left_join(
    diag_df %>%
      filter(day %in% keep_days) %>%
      mutate(day = factor(day, levels = keep_days)) %>%
      select(asv_id, day, r_hat, ess_tail),
    by = c("asv_id","day")
  )

# ---- Extraction for cross-layer activation summary (all originally-tested ASVs, no significance filter) ----
TISSUE <- "Midgut"
LAYER  <- "Digesta_microbiota"

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

activation_summary <- contrast_df %>%
  mutate(day_num = as.numeric(gsub("day_", "", as.character(day)))) %>%
  group_by(day_num) %>%
  summarise(
    boot     = list(boot_median_ci(abs(mean))),
    n_tested = n_distinct(asv_id),
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
  tibble::column_to_rownames("row_id") %>%
  select(tissue, layer, day, median_abs_delta, ci_lower, ci_upper, n_tested)

saveRDS(
  activation_summary,
  file.path(PLOT_OUTDIR, paste0(tolower(LAYER), "_median_abs_delta.rds"))
)
write.csv(
  activation_summary,
  file.path(PLOT_OUTDIR, paste0(tolower(LAYER), "_median_abs_delta.csv")),
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
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "#2ca25f", alpha = 0.4) +
  geom_line(linewidth = 1, color = "#006d2c") +
  geom_point(aes(size = n_tested), color = "#006d2c") +
  scale_x_continuous(breaks = 1:7, labels = levels(combined$day_lab)) +
  scale_size_continuous(range = c(2, 8), name = "n ASVs tested") +
  theme_grey(base_size = 13) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  labs(x = NULL, y = "Median |posterior mean| (bootstrap 95% CI)",
       title = paste(unique(combined$layer), "-", unique(combined$tissue), "effect size dynamics"))

# FLAGS

contrast_df <- contrast_df %>%
  mutate(
    hdi_excl_zero = (!!sym(hdi_low_col) > 0) | (!!sym(hdi_high_col) < 0),
    Pd = PD
  )


# SIGNIFICANCE FILTER (ASV kept if ANY day passes)

keep_asv <- contrast_df %>%
  group_by(asv_id) %>%
  summarise(keep = any(hdi_excl_zero & (Pd > pd_keep_thresh), na.rm = TRUE), .groups = "drop") %>%
  filter(keep) %>%
  pull(asv_id)

heat_df <- contrast_df %>%
  filter(asv_id %in% keep_asv) %>%
  mutate(
    day_categ = as.character(day),
    feature   = as.character(feature)
  )


# CLUSTERING (EFFECT PROFILE ACROSS DAYS)

heat_df2 <- heat_df %>%
  mutate(
    day_num   = as.numeric(gsub("day_", "", day_categ)),
    day_label = paste0("Day ", day_num)
  )

clust_df <- heat_df2 %>%
  select(feature, day_num, estimate = mean) %>%
  arrange(day_num) %>%
  pivot_wider(names_from = day_num, values_from = estimate)

row_names <- clust_df$feature
clust_mat <- as.matrix(clust_df[, -1, drop = FALSE])
rownames(clust_mat) <- row_names
clust_mat[is.na(clust_mat)] <- 0

hc        <- hclust(dist(clust_mat), method = "complete")
ptm_order <- rownames(clust_mat)[hc$order]

ddgram      <- as.dendrogram(hc)
ggtree_plot <- ggtree::ggtree(ddgram, branch.length = "none")

# PLOTTING DATA
heat_df2 <- heat_df %>%
  mutate(
    day_num   = as.numeric(gsub("day_", "", day_categ)),
    day_label = paste0("Day ", day_num),
    feature   = as.character(feature)
  ) %>%
  filter(feature %in% ptm_order)

heat_df2$feature <- factor(heat_df2$feature, levels = ptm_order)

heat_df2$day_label <- factor(
  heat_df2$day_label,
  levels = paste0("Day ", sort(unique(heat_df2$day_num)))
)

heat_df2 <- heat_df2 %>%
  mutate(
    direction  = ifelse(mean >= 0, "HC > NC", "HC < NC"),
    abs_est    = abs(mean),
    prob_label = ifelse(Pd > pd_label_thresh, sprintf("%.2f", Pd), "")
  )


# EXPORT: final ASV list used in the plot

final_asv <- heat_df2 %>%
  distinct(asv_id) %>%
  arrange(asv_id)

# Final ASV list used in the plot
final_asv_path <- file.path(PLOT_OUTDIR, "ASV_used_for_dotplot.tsv")

write_tsv(final_asv, final_asv_path)

cat(
  "Saved final ASV list: ", nrow(final_asv),
  " ASVs -> ", final_asv_path, "\n",
  sep = ""
)


# set colors
cols2 <- c("HC < NC" = "#44a7c4", "HC > NC" = "#f99943")
# plot heatmap
p_dot <- ggplot(heat_df2, aes(x = day_label, y = feature)) +
  geom_point(aes(size = abs_est, fill = direction),
             shape = 21, colour = "black", alpha = 0.9, stroke = 1) +
  scale_fill_manual(
    name   = "Direction",
    values = cols2,
    labels = c("HC < NC", "HC > NC"),
    guide  = guide_legend(override.aes = list(size = 6), order = 1)
  ) +
  scale_size_continuous(
    name  = "|HC − NC| (log2 FC)",
    range = c(4, 20),
    guide = guide_legend(order = 2)
  ) +
  ggnewscale::new_scale("colour") +
  geom_point(aes(size = abs_est, colour = hdi_excl_zero),
             shape = 21, fill = NA, alpha = 1,
             stroke = ifelse(heat_df2$hdi_excl_zero, 2.2, 1)) +
  scale_colour_manual(
    name   = "95% HDI",
    values = c(`FALSE` = "black", `TRUE` = "green4"),
    labels = c(`FALSE` = "includes 0", `TRUE` = "excludes 0"),
    guide  = guide_legend(override.aes = list(size = 6, shape = 21, fill = NA), order = 4)
  ) +
  geom_text(aes(label = prob_label), size = 4.5, colour = "black", vjust = 0.5) +
  geom_point(aes(shape = "Pd text"), x = NA, y = NA, size = 0) +
  scale_shape_manual(
    name   = "Text labels",
    values = c("Pd text" = NA),
    labels = c("Pd (probability of direction)\nshown when Pd > 0.95"),
    guide  = guide_legend(order = 3)
  ) +
  labs(
    #title = "Digesta-Associated Microbiota DA Analysis",
    x = NULL, y = NULL
  ) +
  theme_bw(base_size = 15) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18),
    axis.text.y = element_text(size = 22),
    #plot.title  = element_text(face = "bold", size = 22),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text  = element_text(size = 13),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white"),
    panel.background = element_rect(fill = "grey90")
  )
p_dot <- p_dot +
  labs(title = NULL)

plot_grid_out <- plot_grid(
  ggtree_plot,
  p_dot,
  nrow = 1,
  rel_widths = c(0.1, 2),
  align = "h"
)

plot_grid_out_title <- ggdraw() +
  draw_plot(
    plot_grid_out,
    x = 0,
    y = 0,
    width = 1,
    height = 0.94
  )

print(plot_grid_out_title)

ggsave(
  file.path(PLOT_OUTDIR, "dotplot_diet_effect_digesta_microbiome.tiff"),
  plot_grid_out_title,
  device = "tiff",
  width = 13,
  height = 12,
  dpi = 1200,
  compression = "lzw"
)

ggsave(
  file.path(PLOT_OUTDIR, "dotplot_diet_effect_digesta_microbiome.pdf"),
  plot_grid_out_title,
  device = "pdf",
  width = 13,
  height = 12,
  dpi = 1200
)

# QC PRINT: diagnostics for ASVs actually plotted
qc_df <- heat_df %>%
  distinct(asv_id) %>%
  left_join(
    diag_df %>%
      filter(day %in% keep_days) %>%
      group_by(asv_id) %>%
      summarise(
        worst_rhat   = max(r_hat, na.rm = TRUE),
        min_ess_tail = min(ess_tail, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "asv_id"
  ) %>%
  arrange(desc(worst_rhat), min_ess_tail)

cat("\n================ QC REPORT =================\n")
cat("ASVs retained and plotted:", nrow(qc_df), "\n")
cat("Rhat threshold:", rhat_bad, "\n")
cat("ESS_tail threshold:", ess_tail_bad, "\n")
cat("-------------------------------------------\n")

print(
  qc_df %>%
    mutate(
      worst_rhat   = round(worst_rhat, 5),
      min_ess_tail = round(min_ess_tail, 0)
    ),
  n = Inf
)

cat("-------------------------------------------\n")
cat("Summary:\n")
cat("  Max Rhat (retained): ", round(max(qc_df$worst_rhat, na.rm = TRUE), 5), "\n")
cat("  Min ESS_tail (retained): ", round(min(qc_df$min_ess_tail, na.rm = TRUE), 0), "\n")
cat("===========================================\n\n")

# hard safety check (optional but recommended)
stopifnot(
  max(qc_df$worst_rhat, na.rm = TRUE) < rhat_bad,
  min(qc_df$min_ess_tail, na.rm = TRUE) > ess_tail_bad
)

