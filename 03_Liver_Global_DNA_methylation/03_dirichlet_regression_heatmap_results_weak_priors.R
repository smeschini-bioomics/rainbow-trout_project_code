suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(ggnewscale)
})

# Input/output directory for weak-prior sensitivity analysis

out_dir <- "results/dirichlet_regression_weak_priors"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Reload contrasts from weak-prior model

contr_summ <- readr::read_csv(
  file.path(
    out_dir,
    "contrasts_HC_minus_NC_by_day_response.csv"
  ),
  show_col_types = FALSE
)

# Clean and order sampling days

contr_summ <- contr_summ %>%
  mutate(
    day = str_trim(as.character(day)),
    day_num = as.integer(
      str_match(
        day,
        "(?i)day[_-]?(\\d+)"
      )[, 2]
    ),
    day_num = ifelse(
      is.na(day_num),
      match(day, unique(day)),
      day_num
    )
  ) %>%
  arrange(day_num) %>%
  mutate(
    day = factor(
      day,
      levels = unique(day)
    ),
    direction = factor(
      direction,
      levels = c(
        "decrease",
        "increase",
        "neutral"
      )
    )
  )

# Cross-layer summary: liver DNA 5-mC

activation_summary_methylation <- contr_summ %>%
  filter(response == "mdC_rel") %>%
  transmute(
    tissue = "Liver",
    layer = "DNA 5-mC",
    day = day_num,
    effect = estimate_pp,
    ci_lower = lower95_pp,
    ci_upper = upper95_pp,
    n_tested = 1L
  ) %>%
  arrange(day)

saveRDS(
  activation_summary_methylation,
  file.path(
    out_dir,
    "liver_5mc_effect_summary_weak_priors.rds"
  )
)

readr::write_csv(
  activation_summary_methylation,
  file.path(
    out_dir,
    "liver_5mc_effect_summary_weak_priors.csv"
  )
)

# Check plot: HC − NC difference in 5-mC

day_order <- c(
  1,
  2,
  3,
  4,
  10,
  15,
  22
)

combined <- activation_summary_methylation %>%
  mutate(
    day_idx = match(
      day,
      day_order
    )
  )

p_5mc <- ggplot(
  combined,
  aes(
    x = day_idx,
    y = effect,
    group = layer
  )
) +
  geom_hline(
    yintercept = 0,
    color = "grey60",
    linewidth = 0.5
  ) +
  geom_errorbar(
    aes(
      ymin = ci_lower,
      ymax = ci_upper
    ),
    width = 0.12,
    color = "#7B3294"
  ) +
  geom_line(
    linewidth = 1,
    color = "#7B3294"
  ) +
  geom_point(
    size = 3,
    color = "#7B3294"
  ) +
  scale_x_continuous(
    breaks = seq_along(day_order),
    labels = paste0(
      "Day ",
      day_order
    )
  ) +
  labs(
    x = NULL,
    y = paste0(
      "HC − NC difference in 5-mC\n",
      "(percentage points)"
    ),
    title = "Liver DNA 5-mC: weak-prior sensitivity"
  ) +
  theme_classic(
    base_size = 13
  )

p_5mc

ggsave(
  filename = file.path(
    out_dir,
    "liver_5mc_effect_weak_priors.pdf"
  ),
  plot = p_5mc,
  width = 7.5,
  height = 5.5
)

ggsave(
  filename = file.path(
    out_dir,
    "liver_5mc_effect_weak_priors.tiff"
  ),
  plot = p_5mc,
  width = 7.5,
  height = 5.5,
  dpi = 1200,
  compression = "lzw"
)

# Response labels and row order

contr_plot <- contr_summ %>%
  mutate(
    response_label = recode(
      response,
      dC_rel = "C (%)",
      mdC_rel = "5-mC (%)",
      hmdC_rel = "5-hmC (%)"
    ),
    response_label = factor(
      response_label,
      levels = c(
        "5-hmC (%)",
        "5-mC (%)",
        "C (%)"
      )
    )
  ) %>%
  filter(direction != "neutral") %>%
  droplevels()

# Pretty day labels

day_levels_orig <- levels(
  contr_plot$day
)

day_labels_pretty <- setNames(
  sub(
    "(?i)^day[_-]?0*(\\d+)$",
    "Day \\1",
    day_levels_orig,
    perl = TRUE
  ),
  day_levels_orig
)

# Heatmap

p_heat <- ggplot(
  contr_plot,
  aes(
    x = day,
    y = response_label
  )
) +
  geom_tile(
    fill = "grey95",
    color = "white",
    linewidth = 0.3
  ) +
  geom_point(
    data = subset(
      contr_plot,
      abs_est_pp_clean > 0
    ),
    aes(
      size = abs_est_pp_clean,
      fill = direction
    ),
    color = "grey20",
    shape = 21,
    stroke = 0.4,
    alpha = 0.95
  ) +
  scale_fill_manual(
    values = c(
      decrease = "#44a7c4",
      increase = "#f99943"
    ),
    labels = c(
      decrease = "HC < NC",
      increase = "HC > NC"
    ),
    name = "Direction",
    guide = guide_legend(
      override.aes = list(
        size = 6
      ),
      order = 1
    )
  ) +
  scale_size_continuous(
    range = c(
      5,
      12
    ),
    name = "|HC − NC| (%)",
    guide = guide_legend(
      order = 2
    )
  ) +
  ggnewscale::new_scale("colour") +
  geom_point(
    aes(
      size = abs_est_pp_clean,
      colour = sig95
    ),
    shape = 21,
    fill = NA,
    alpha = 1,
    stroke = ifelse(
      contr_plot$sig95,
      2.2,
      1
    )
  ) +
  scale_colour_manual(
    name = "95% CrI",
    values = c(
      `FALSE` = "black",
      `TRUE` = "green4"
    ),
    labels = c(
      `FALSE` = "includes 0",
      `TRUE` = "excludes 0"
    ),
    guide = guide_legend(
      override.aes = list(
        size = 6,
        shape = 21,
        fill = NA
      ),
      order = 4
    )
  ) +
  geom_text(
    data = subset(
      contr_plot,
      prob_certainty > 0.95
    ),
    aes(
      label = label_text
    ),
    color = "black",
    size = 7,
    fontface = "bold",
    vjust = -1.2
  ) +
  geom_point(
    aes(
      shape = "Pd text"
    ),
    x = NA,
    y = NA,
    size = 0
  ) +
  scale_shape_manual(
    name = "Text labels",
    values = c(
      "Pd text" = NA
    ),
    labels = c(
      "Pd text" =
        "Pd (probability of direction)\nshown when Pd > 0.95"
    ),
    guide = guide_legend(
      order = 3
    )
  ) +
  scale_x_discrete(
    labels = day_labels_pretty
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Weak-prior sensitivity analysis"
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_text(
      face = "bold",
      size = 14,
      angle = 90,
      hjust = 0.5
    ),
    axis.text.x = element_text(
      size = 16,
      angle = 90,
      hjust = 0.5
    ),
    plot.title = element_text(
      family = "mono",
      face = "bold",
      size = 20,
      hjust = 0.5
    ),
    legend.title = element_text(
      size = 12,
      face = "bold"
    ),
    legend.text = element_text(
      size = 12
    )
  )

p_heat

# Save heatmap in the weak-prior directory

ggsave(
  filename = file.path(
    out_dir,
    "heatmap_centered_probabilities_weak_priors.pdf"
  ),
  plot = p_heat,
  width = 7.5,
  height = 5.5
)

ggsave(
  filename = file.path(
    out_dir,
    "heatmap_centered_probabilities_weak_priors.tiff"
  ),
  plot = p_heat,
  width = 7.5,
  height = 5.5,
  dpi = 1200,
  compression = "lzw"
)