suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(ggnewscale)
})

# 1) Reload
contr_summ <- readr::read_csv(
  "results/dirichlet_regression/contrasts_HC_minus_NC_by_day_response.csv",
  show_col_types = FALSE
)

# 2) Clean & order day properly
contr_summ <- contr_summ %>%
  mutate(
    day = str_trim(as.character(day)),
    day_num = as.integer(str_match(day, "(?i)day[_-]?(\\d+)")[, 2]),
    day_num = ifelse(is.na(day_num), match(day, unique(day)), day_num)
  ) %>%
  arrange(day_num) %>%
  mutate(
    day = factor(day, levels = unique(day)),
    direction = factor(direction, levels = c("decrease", "increase", "neutral"))
  )


# ---- Extraction for cross-layer activation summary: liver DNA methylation ----
# Compositional data (C / 5-mC / 5-hmC sum to ~100%) — using only 5-mC
# as the single representative measure, since the three responses are
# not independent and 5-mC is the primary, most interpretable signal.

# Extraction for cross-layer summary: signed liver 5-mC effect

activation_summary_methylation <- contr_summ %>%
  filter(response == "mdC_rel") %>%
  transmute(
    tissue   = "Liver",
    layer    = "DNA 5-mC",
    day      = day_num,
    effect   = estimate_pp,
    ci_lower = lower95_pp,
    ci_upper = upper95_pp,
    n_tested = 1L
  ) %>%
  arrange(day)

saveRDS(
  activation_summary_methylation,
  "results/dirichlet_regression/liver_5mc_effect_summary.rds"
)

write.csv(
  activation_summary_methylation,
  "results/dirichlet_regression/liver_5mc_effect_summary.csv",
  row.names = FALSE
)

# Check plot

combined <- activation_summary_methylation %>%
  mutate(
    day_idx = match(
      day,
      c(1, 2, 3, 4, 10, 15, 22)
    )
  )

ggplot(
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
    breaks = 1:7,
    labels = paste0(
      "Day ",
      c(1, 2, 3, 4, 10, 15, 22)
    )
  ) +
  labs(
    x = NULL,
    y = "HC − NC difference in 5-mC\n(percentage points)",
    title = "Liver DNA 5-mC"
  ) +
  theme_classic(
    base_size = 13
  )

# 3) Response labels & row order
contr_plot <- contr_summ %>%
  mutate(
    response_label = recode(
      response,
      dC_rel   = "C (%)",
      mdC_rel  = "5-mC (%)",
      hmdC_rel = "5-hmC (%)"
    ),
    response_label = factor(
      response_label,
      levels = c("5-hmC (%)", "5-mC (%)", "C (%)")
    )
  ) %>%
  filter(direction != "neutral") %>%
  droplevels()

# 4) Pretty day labels
day_levels_orig <- levels(contr_plot$day)
day_labels_pretty <- setNames(
  sub("(?i)^day[_-]?0*(\\d+)$", "Day \\1", day_levels_orig, perl = TRUE),
  day_levels_orig
)

# 5) Plot with CrI outline logic
p_heat <- ggplot(contr_plot, aes(x = day, y = response_label)) +
  geom_tile(fill = "grey95", color = "white", linewidth = 0.3) +
  geom_point(
    aes(size = abs_est_pp_clean, fill = direction),
    color  = "grey20",
    shape  = 21,
    stroke = 0.4,
    alpha  = 0.95,
    data   = subset(contr_plot, abs_est_pp_clean > 0)
  ) +
  scale_fill_manual(
    values = c(
      decrease = "#44a7c4",
      increase = "#f99943"
    ),
    labels = c(
      "HC < NC",
      "HC > NC"
    ),
    name = "Direction",
    guide = guide_legend(
      override.aes = list(size = 6),
      order = 1
    )
  ) +
  scale_size_continuous(
    range = c(5, 12),
    name = "|HC − NC| (%)",
    guide = guide_legend(order = 2)
  ) +
  ggnewscale::new_scale("colour") +
  geom_point(
    aes(
      size   = abs_est_pp_clean,
      colour = sig95
    ),
    shape  = 21,
    fill   = NA,
    alpha  = 1,
    stroke = ifelse(contr_plot$sig95, 2.2, 1)
  ) +
  scale_colour_manual(
    name = "95% CrI",
    values = c(
      `FALSE` = "black",
      `TRUE`  = "green4"
    ),
    labels = c(
      `FALSE` = "includes 0",
      `TRUE`  = "excludes 0"
    ),
    guide = guide_legend(
      override.aes = list(size = 6, shape = 21, fill = NA),
      order = 4
    )
  ) +
  geom_text(
    aes(label = label_text),
    data     = subset(contr_plot, prob_certainty > 0.95),
    color    = "black",
    size     = 7,
    fontface = "bold",
    vjust    = -1.2
  ) +
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
  scale_x_discrete(labels = day_labels_pretty) +
  labs(
    x = NULL,
    y = NULL
  ) +
  labs(
    #title = "Global Liver DNA Methylation:\nDirichlet Regression"
  ) +
  theme(
    axis.text.y  = element_text(face = "bold", size = 14, angle = 90, hjust = 0.5),
    plot.title = element_text(
      family = "mono",
      face = "bold",
      size = 20,
      hjust = 0.5),
    axis.text.x  = element_text(size = 16, angle = 90, hjust = 0.5),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 12)
  )

p_heat

# Save
out_dir <- "results/dirichlet_regression"

ggsave(
  filename = file.path(out_dir, "heatmap_centered_probabilities.pdf"),
  plot = p_heat,
  width = 7.5,
  height = 5.5
)

ggsave(
  filename = file.path(out_dir, "heatmap_centered_probabilities.tiff"),
  plot = p_heat,
  width = 7.5,
  height = 5.5,
  dpi = 1200
)
