suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ggpubr)
})


# Output directory
out_dir <- "results/CV_analysis"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Figures will be saved in: ", normalizePath(out_dir))


# Colors
diet_cols <- c(
  "NC" = "#44a7c4",
  "HC" = "#f99943"
)


# Load saved ggplot objects
p_BW      <- readRDS("RDS/body_weight_plot.rds")
p_glucose <- readRDS("RDS/plasma_glucose_plot.rds")
p_HSI     <- readRDS("RDS/hepatosomatic_index_plot.rds")
p_L_La    <- readRDS("RDS/plasma_llactate_plot.rds")
p_D_La    <- readRDS("RDS/plasma_dlactate_plot.rds")

# Extract raw data from the ggplot objects
bw_df <- p_BW$data %>%
  select(diet, day, value = body_weight) %>%
  mutate(trait = "Body Weight")

glu_df <- p_glucose$data %>%
  select(diet, day, value = plasma_glucose) %>%
  mutate(trait = "Plasma Glucose")

hsi_df <- p_HSI$data %>%
  select(diet, day, value = HSI) %>%
  mutate(trait = "HSI")

lla_df <- p_L_La$data %>%
  select(diet, day, value = Plasma_L_Lactate) %>%
  mutate(trait = "Plasma L-Lactate")

dla_df <- p_D_La$data %>%
  select(diet, day, value = Plasma_D_Lactate) %>%
  mutate(trait = "Plasma D-Lactate")

traits_df <- bind_rows(
  bw_df,
  glu_df,
  hsi_df,
  lla_df,
  dla_df
) %>%
  filter(
    !is.na(value),
    !is.na(diet),
    !is.na(day)
  ) %>%
  mutate(
    trait = factor(
      trait,
      levels = c(
        "Body Weight",
        "Plasma Glucose",
        "HSI",
        "Plasma L-Lactate",
        "Plasma D-Lactate"
      )
    ),
    diet = factor(diet, levels = c("NC", "HC"))
  )


# Compute CV per day within each diet and trait
cv_day <- traits_df %>%
  group_by(trait, diet, day) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    sd_value = sd(value, na.rm = TRUE),
    CV_percent = 100 * sd_value / mean_value,
    .groups = "drop"
  ) %>%
  filter(is.finite(CV_percent))

print(cv_day)


# Jittered boxplots of day-specific CV values by diet
p_cv <- ggplot(
  cv_day,
  aes(x = diet, y = CV_percent, fill = diet)
) +
  geom_boxplot(
    width = 0.35,
    color = "black",
    linewidth = 0.4,
    outlier.shape = NA,
    alpha = 0.95
  ) +
  geom_point(
    position = position_jitter(width = 0.06, height = 0),
    shape = 21,
    size = 1.8,
    color = "black",
    stroke = 0.6
  ) +
  facet_wrap(~trait, scales = "free_y") +
  scale_fill_manual(
    name = "Diet",
    values = diet_cols
  ) +
  stat_compare_means(
    comparisons = list(c("NC", "HC")),
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = FALSE,
    size = 4.8,
    bracket.size = 0.8
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.18))
  ) +
  labs(
    x = "Diet",
    y = "Coefficient of variation (%)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 12, face = "bold"),
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_text(face = "bold", size = 16)
  )

print(p_cv)


# Save figure
ggsave(
  filename = file.path(out_dir, "host_traits_CV_boxplot.pdf"),
  plot = p_cv,
  width = 10,
  height = 6,
  dpi = 1200
)

ggsave(
  filename = file.path(out_dir, "host_traits_CV_boxplot.tiff"),
  plot = p_cv,
  width = 10,
  height = 6,
  dpi = 1200,
  compression = "lzw"
)

print(list.files(out_dir, full.names = TRUE))
