suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggbeeswarm)
})

# ----------------------------- [1] INPUTS ---------------------------------

infile  <- "input_files/methylome_for_dirichlet.tsv"
out_dir <- "results/beeswarm_plot"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# ----------------------------- [2] LOAD DATA ------------------------------

df <- readr::read_tsv(infile, show_col_types = FALSE)

# Check required columns
need <- c(
  "diet", "day",
  "dC_rel", "mdC_rel", "hmdC_rel"
)

miss <- setdiff(need, names(df))

if (length(miss) > 0) {
  stop("Missing columns: ", paste(miss, collapse = ", "))
}

# Remove Fasted group and prepare variables
df <- df %>%
  filter(diet != "Fasted") %>%
  mutate(
    diet = factor(diet, levels = c("NC", "HC")),
    
    day = factor(
      day,
      levels = c(
        "day_01", "day_02", "day_03", "day_04",
        "day_10", "day_15", "day_22"
      ),
      labels = c(
        "Day 1", "Day 2", "Day 3", "Day 4",
        "Day 10", "Day 15", "Day 22"
      )
    ),
    
    dC_pct   = 100 * dC_rel,
    mdC_pct  = 100 * mdC_rel,
    hmdC_pct = 100 * hmdC_rel
  ) %>%
  droplevels()

# ----------------------------- [3] LONG FORMAT ----------------------------

df_long <- df %>%
  select(
    diet, day,
    dC_pct, mdC_pct, hmdC_pct
  ) %>%
  pivot_longer(
    cols = c(dC_pct, mdC_pct, hmdC_pct),
    names_to = "response",
    values_to = "value"
  ) %>%
  mutate(
    x_pos = case_when(
      diet == "NC" ~ 1,
      diet == "HC" ~ 3
    )
  )

# ----------------------------- [3b] SUMMARY TABLES ------------------------

df_long_sum <- df_long %>%
  mutate(
    response_label = recode(
      as.character(response),
      dC_pct   = "C",
      mdC_pct  = "5-mC",
      hmdC_pct = "5-hmC"
    )
  )

# Sample numbers per day × diet × cytosine form
n_per_group <- df_long_sum %>%
  group_by(response_label, day, diet) %>%
  summarise(
    n_total = n(),
    n_used = sum(!is.na(value)),
    .groups = "drop"
  ) %>%
  arrange(response_label, day, diet)

write_csv(
  n_per_group,
  file.path(out_dir, "n_per_group_for_raw_beeswarm_plots.csv")
)

print(n_per_group, n = Inf)

# NC and HC side-by-side sample-number table
n_per_group_wide <- n_per_group %>%
  select(response_label, day, diet, n_used) %>%
  pivot_wider(
    names_from = diet,
    values_from = n_used
  ) %>%
  arrange(response_label, day)

write_csv(
  n_per_group_wide,
  file.path(out_dir, "n_per_group_wide_for_raw_beeswarm_plots.csv")
)

print(n_per_group_wide, n = Inf)

# Summary statistics
summary_table <- df_long_sum %>%
  group_by(response_label, day, diet) %>%
  summarise(
    n_used = sum(!is.na(value)),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    q1 = quantile(value, 0.25, na.rm = TRUE),
    q3 = quantile(value, 0.75, na.rm = TRUE),
    min = min(value, na.rm = TRUE),
    max = max(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(response_label, day, diet)

write_csv(
  summary_table,
  file.path(out_dir, "summary_table_raw_beeswarm_data.csv")
)

print(summary_table, n = Inf)

# Global summary per cytosine form
global_summary <- df_long_sum %>%
  group_by(response_label) %>%
  summarise(
    n_used = sum(!is.na(value)),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    min = min(value, na.rm = TRUE),
    max = max(value, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  global_summary,
  file.path(out_dir, "global_summary_table_raw_beeswarm_data.csv")
)

print(global_summary, n = Inf)

# ----------------------------- [4] PLOT SETTINGS --------------------------

pal <- c(
  NC = "#44a7c4",
  HC = "#f99943"
)

make_beeswarm_plot <- function(dat, plot_title, ylab) {
  
  ggplot(dat, aes(x = x_pos, y = value)) +
    ggbeeswarm::geom_beeswarm(
      aes(color = diet, group = diet),
      method = "swarm",
      priority = "density",
      cex = 5,
      size = 3,
      alpha = 0.55,
      corral = "none",
      preserve.data.axis = TRUE,
      na.rm = TRUE
    ) +
    facet_wrap(
      ~ day,
      nrow = 1,
      scales = "free_x"
    ) +
    scale_color_manual(
      values = pal,
      name = "Diet"
    ) +
    scale_x_continuous(
      limits = c(0, 4),
      breaks = c(1, 3),
      labels = c("", "")
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.18))
    ) +
    labs(
      title = plot_title,
      x = NULL,
      y = ylab,
      color = "Diet"
    ) +
    theme_bw(base_size = 16) +
    theme(
      plot.title = element_text(
        size = 20,
        face = "bold",
        hjust = 0.5,
        family = "mono"
      ),
      strip.text = element_text(
        size = 12,
      ),
      legend.position = "right",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 18),
      axis.text.y = element_text(size = 14),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 16),
      panel.grid = element_blank()
    )
}

# ----------------------------- [5] CREATE PLOTS ---------------------------

p_dC <- make_beeswarm_plot(
  dat = filter(df_long, response == "dC_pct"),
  plot_title = "Global Liver Unmodified C",
  ylab = "%"
)
p_dC

p_mdC <- make_beeswarm_plot(
  dat = filter(df_long, response == "mdC_pct"),
  plot_title = "Global Liver 5-mC",
  ylab = "%"
)
p_mdC

p_hmdC <- make_beeswarm_plot(
  dat = filter(df_long, response == "hmdC_pct"),
  plot_title = "Global Liver 5-hmC",
  ylab = "%"
)
p_hmdC
# ----------------------------- [6] SAVE FIGURES ----------------------------

ggsave(
  filename = file.path(out_dir, "raw_beeswarm_dC.pdf"),
  plot = p_dC,
  device = "pdf",
  width = 6.7,
  height = 3.5,
  dpi = 1200
)

ggsave(
  filename = file.path(out_dir, "raw_beeswarm_dC.tiff"),
  plot = p_dC,
  device = "tiff",
  width = 6.7,
  height = 3.5,
  dpi = 1200,
  compression = "lzw"
)

ggsave(
  filename = file.path(out_dir, "raw_beeswarm_mdC.pdf"),
  plot = p_mdC,
  device = "pdf",
  width = 6.7,
  height = 3.5,
  dpi = 1200
)

ggsave(
  filename = file.path(out_dir, "raw_beeswarm_mdC.tiff"),
  plot = p_mdC,
  device = "tiff",
  width = 6.7,
  height = 3.5,
  dpi = 1200,
  compression = "lzw"
)

ggsave(
  filename = file.path(out_dir, "raw_beeswarm_hmdC.pdf"),
  plot = p_hmdC,
  device = "pdf",
  width = 6.7,
  height = 3.5,
  dpi = 1200
)

ggsave(
  filename = file.path(out_dir, "raw_beeswarm_hmdC.tiff"),
  plot = p_hmdC,
  device = "tiff",
  width = 6.7,
  height = 3.5,
  dpi = 1200,
  compression = "lzw"
)

# ----------------------------- [7] SESSION INFO ----------------------------

writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "sessionInfo_raw_beeswarm.txt")
)

message(
  "Done.\n",
  "Created in folder: ", out_dir, "\n",
  "- raw_beeswarm_dC.pdf / .tiff\n",
  "- raw_beeswarm_mdC.pdf / .tiff\n",
  "- raw_beeswarm_hmdC.pdf / .tiff\n",
  "- summary tables and session information\n"
)


library(cowplot)

# ----------------------------- [5b] COMBINE IN ONE COLUMN -----------------

# Combined-figure versions: no titles and specific y-axis labels
p_dC_combined <- p_dC +
  labs(
    title = NULL,
    y = "C (%)"
  ) +
  theme(
    plot.margin = margin(t = 4, r = 4, b = 4, l = 4)
  )

p_mdC_combined <- p_mdC +
  labs(
    title = NULL,
    y = "5-mC (%)"
  ) +
  theme(
    plot.margin = margin(t = 4, r = 4, b = 4, l = 4)
  )

p_hmdC_combined <- p_hmdC +
  labs(
    title = NULL,
    y = "5-hmC (%)"
  ) +
  theme(
    plot.margin = margin(t = 4, r = 4, b = 4, l = 4)
  )

# Extract one shared legend
shared_legend <- cowplot::get_legend(
  p_dC_combined +
    theme(legend.position = "right")
)

# Remove legends from panels
p_dC_no_legend <- p_dC_combined +
  theme(legend.position = "none")

p_mdC_no_legend <- p_mdC_combined +
  theme(legend.position = "none")

p_hmdC_no_legend <- p_hmdC_combined +
  theme(legend.position = "none")

# Stack panels vertically
plots_one_column <- cowplot::plot_grid(
  p_dC_no_legend,
  p_mdC_no_legend,
  p_hmdC_no_legend,
  ncol = 1,
  align = "v",
  axis = "lr"
)

# Add shared legend
p_methylome_combined <- cowplot::plot_grid(
  plots_one_column,
  shared_legend,
  ncol = 2,
  rel_widths = c(1, 0.12)
)

p_methylome_combined

# ----------------------------- [6b] SAVE COMBINED FIGURE ------------------

ggsave(
  filename = file.path(out_dir, "raw_beeswarm_methylome_combined.pdf"),
  plot = p_methylome_combined,
  device = "pdf",
  width = 7.5,
  height = 10,
  dpi = 1200
)

ggsave(
  filename = file.path(out_dir, "raw_beeswarm_methylome_combined.tiff"),
  plot = p_methylome_combined,
  device = "tiff",
  width = 7.5,
  height = 10,
  dpi = 1200,
  compression = "lzw"
)

