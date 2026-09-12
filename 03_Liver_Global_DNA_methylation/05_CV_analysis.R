suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(purrr)
})

# SETTING
infile  <- "input_files/methylome_for_dirichlet.tsv"
out_dir <- "results/CV_analysis"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

pal <- c(
  NC = "#44a7c4",
  HC = "#f99943"
)

# HELPER FUNCTION
boot_cv <- function(x, n_boot = 1000) {
  x <- x[is.finite(x) & !is.na(x)]
  
  if (length(x) < 2) {
    return(rep(NA_real_, n_boot))
  }
  
  replicate(n_boot, {
    xb <- sample(x, size = length(x), replace = TRUE)
    m  <- mean(xb, na.rm = TRUE)
    s  <- sd(xb, na.rm = TRUE)
    
    if (is.na(m) || m == 0) {
      NA_real_
    } else {
      100 * s / m
    }
  })
}

# LOAD DATA
df_cv <- readr::read_tsv(infile) %>%
  filter(diet != "Fasted") %>%
  mutate(
    diet = factor(diet, levels = c("NC", "HC")),
    day = factor(
      day,
      levels = c("day_01", "day_02", "day_03", "day_04",
                 "day_10", "day_15", "day_22"),
      labels = c("Day 1", "Day 2", "Day 3", "Day 4",
                 "Day 10", "Day 15", "Day 22")
    )
  )

df_long <- df_cv %>%
  pivot_longer(
    cols = c(dC_rel, mdC_rel, hmdC_rel),
    names_to = "response",
    values_to = "value"
  ) %>%
  mutate(
    response_label = factor(
      response,
      levels = c("dC_rel", "mdC_rel", "hmdC_rel"),
      labels = c("C", "5-mC", "5-hmC")
    )
  )

# CV COMPUTATION
cv_day <- df_long %>%
  filter(!is.na(value)) %>%
  group_by(response_label, diet, day) %>%
  summarise(
    n = n(),
    mean_value = mean(value, na.rm = TRUE),
    sd_value   = sd(value, na.rm = TRUE),
    CV_percent = 100 * sd_value / mean_value,
    .groups = "drop"
  ) %>%
  filter(is.finite(CV_percent))

# Save table
readr::write_csv(
  cv_day,
  file.path(out_dir, "cv_percent_relative_all_data_by_day_diet.csv")
)

print(cv_day, n = Inf)

# P-VALUES
stat.test <- compare_means(
  CV_percent ~ diet,
  data = cv_day,
  group.by = "response_label",
  method = "wilcox.test"
) %>%
  filter(p <= 0.05) %>%
  left_join(
    cv_day %>%
      group_by(response_label) %>%
      summarise(
        y.position = max(CV_percent, na.rm = TRUE) * 1.08,
        .groups = "drop"
      ),
    by = "response_label"
  )

# PLOT
p_cv <- ggplot(cv_day, aes(x = diet, y = CV_percent, fill = diet)) +
  geom_boxplot(
    width = 0.35,
    color = "black",
    linewidth = 0.6,
    outlier.shape = NA,
    alpha = 0.95
  ) +
  geom_point(
    position = position_jitter(width = 0.06),
    shape = 21,
    size = 2.8,
    color = "black",
    stroke = 0.6
  ) +
  facet_wrap(~ response_label, scales = "free_y") +
  scale_fill_manual(values = pal, guide = "none") +
  labs(
    x = NULL,
    y = "Coefficient of variation (%)"
  ) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 16, face = "bold"),
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_text(face = "bold", size = 16)
  )

p_final <- p_cv +
  stat_pvalue_manual(
    stat.test,
    label = "p.signif",
    tip.length = 0.01,
    inherit.aes = FALSE
  )

print(p_final)

# SAVE GLOBAL PLOT
ggsave(
  file.path(out_dir, "CV_boxjitter_relative_all_data.pdf"),
  p_final,
  width = 5,
  height = 3
)

ggsave(
  file.path(out_dir, "CV_boxjitter_relative_all_data.tiff"),
  p_final,
  width = 5,
  height = 3,
  dpi = 1200,
  compression = "lzw"
)

# N TABLE
n_table <- df_long %>%
  filter(!is.na(value)) %>%
  group_by(response_label, diet) %>%
  summarise(n = n(), .groups = "drop")

readr::write_csv(
  n_table,
  file.path(out_dir, "n_CV_by_diet_methylome_relative_all_data.csv")
)

