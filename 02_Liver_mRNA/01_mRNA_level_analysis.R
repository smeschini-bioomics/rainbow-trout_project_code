# Libraries
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(glmmTMB)
  library(emmeans)
  library(DHARMa)
  library(ggsignif)
  library(ggbeeswarm)
})

# Load
file_path <- "input_files/<file name.tsv>"

analysis_name <- tools::file_path_sans_ext(basename(file_path))
file_base <- analysis_name

results_dir <- file.path("results", analysis_name)
rds_dir <- "RDS"

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

# inspet file
data <- read_tsv(file_path)
data <- as.data.frame(data)
str(data)

# Convert variables to appropriate factor levels
data <- data %>%
  mutate(
    day = factor(day, levels = c("0", "1", "2", "3", "4", "10", "15", "22")),
    diet = factor(diet, levels = c("Fasted", "NC", "HC"))
  )

# Filter out Fasted group
data_no_fasted <- data %>%
  filter(diet != "Fasted") %>%
  droplevels()

str(data_no_fasted)

# -------------------------------------------------------------------------
# Minimum replication requirement for inferential analysis
# Require at least 3 non-missing RC values in every planned Diet × Day cell.
# -------------------------------------------------------------------------
analysis_days <- c("1", "2", "3", "4", "10", "15", "22")

data_no_fasted <- data_no_fasted %>%
  filter(
    !is.na(RC),
    !is.na(day),
    !is.na(diet)
  )

group_counts <- tidyr::expand_grid(
  day = analysis_days,
  diet = c("NC", "HC")
) %>%
  left_join(
    data_no_fasted %>%
      transmute(
        day = as.character(day),
        diet = as.character(diet)
      ) %>%
      count(day, diet, name = "n"),
    by = c("day", "diet")
  ) %>%
  mutate(
    n = dplyr::coalesce(n, 0L),
    passes_minimum_n = n >= 3L
  )

write_csv(
  group_counts,
  file.path(results_dir, paste0(file_base, "_group_counts.csv"))
)

insufficient_groups <- group_counts %>%
  filter(!passes_minimum_n)

run_inference <- nrow(insufficient_groups) == 0L

if (!run_inference) {
  
  write_csv(
    insufficient_groups,
    file.path(
      results_dir,
      paste0(file_base, "_insufficient_group_sizes.csv")
    )
  )
  
  writeLines(
    c(
      "Inferential analysis skipped.",
      "At least one Diet × Day group had fewer than 3 non-missing RC values.",
      "",
      capture.output(print(insufficient_groups))
    ),
    file.path(
      results_dir,
      paste0(file_base, "_inferential_analysis_skipped.txt")
    )
  )
  
  message(
    "Inferential analysis skipped for ", file_base,
    ": at least one NC/HC × day group has n < 3. ",
    "The plot without statistical annotations will still be generated."
  )
}

# Data Exploratory Visualization

ggplot(data_no_fasted, aes(x = RC)) +
  # Histogram of RC (scaled to density)
  geom_histogram(aes(y = ..density..), bins = 30, fill = "steelblue", alpha = 0.5) +
  
  # Kernel density estimate overlay
  geom_density(color = "red", linewidth = 1.2) +
  
  # Vertical line for the mean
  geom_vline(aes(xintercept = mean(RC, na.rm = TRUE)), 
             color = "darkgreen", linetype = "dashed", linewidth = 1) +
  
  # Vertical line for the median
  geom_vline(aes(xintercept = median(RC, na.rm = TRUE)), 
             color = "orange", linetype = "dotted", linewidth = 1) +
  
  # Labels and theme
  labs(
    title = "Distribution of mRNA Relative Concentration (RC)",
    subtitle = "Histogram with density curve, mean (dashed), and median (dotted)",
    x = "GOI Relative Concentration",
    y = "Density"
  ) +
  
  theme_minimal(base_size = 14)



# Pre-calculate values
mean_rc <- mean(data_no_fasted$RC, na.rm = TRUE)
median_rc <- median(data_no_fasted$RC, na.rm = TRUE)
logmean_rc <- exp(mean(log(data_no_fasted$RC), na.rm = TRUE))

# Plot
ggplot(data_no_fasted, aes(x = RC)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "steelblue", alpha = 0.5) +
  geom_density(color = "red", linewidth = 1.2) +
  
  # Now use the pre-calculated values
  geom_vline(xintercept = mean_rc, color = "darkgreen", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = median_rc, color = "orange", linetype = "dotted", linewidth = 1) +
  geom_vline(xintercept = logmean_rc, color = "purple", linetype = "dotdash", linewidth = 1) +
  
  geom_rug(alpha = 0.2) +
  facet_wrap(~ diet, scales = "free") +
  labs(
    title = "Distribution of GOI Relative Concentration by diet group",
    subtitle = "Histogram with density, mean (dashed), median (dotted), log-mean (dotdash), and individual data (rug)",
    x = "mRNA RC / rna18s eef1a",
    y = "Density"
  ) +
  theme_minimal(base_size = 14)

# Fit Candidate Models

# Gaussian Models
model_gaussian_RC <- glmmTMB(RC ~ day * diet,
                             data = data_no_fasted,
                             family = gaussian(link = "identity"))
# Gaussian model fitted to log-transformed RC
model_gaussian_log_RC <- glmmTMB(log(RC) ~ day * diet,
                             data = data_no_fasted, family = gaussian(link = "identity"))
#logNormal family
model_logNormal_RC <- glmmTMB(RC ~ day * diet,
                          data = data_no_fasted, family = lognormal(link = "log"))
# t_family models
model_student_RC <- glmmTMB(RC ~ day * diet,
                          data = data_no_fasted, family = t_family(link = "identity"))
# Gamma Models
model_gamma_RC <- glmmTMB(RC ~ day * diet,
                          data = data_no_fasted, family = Gamma(link = "log"))


# Residual Diagnostics

res_gaussian_RC <- simulateResiduals(model_gaussian_RC,
                                     n = 1000,
                                     refit = FALSE)

res_gaussian_log_RC <- simulateResiduals(model_gaussian_log_RC,
                                         n = 1000,
                                         refit = FALSE)

res_logNormal_RC <- simulateResiduals(model_logNormal_RC,
                                      n = 1000,
                                      refit = FALSE)

res_student_RC <- simulateResiduals(model_student_RC,
                                    n = 1000,
                                    refit = FALSE)


res_gamma_RC <- simulateResiduals(model_gamma_RC,
                                  n = 1000,
                                  refit = FALSE)

plot(res_gaussian_RC)
plot(res_gaussian_log_RC)
plot(res_logNormal_RC)
plot(res_student_RC)
plot(res_gamma_RC)

# Pairwise contrasts

# -------------------------------------------------------------------------
# Select the final model after reviewing diagnostics for the current gene.
# Valid choices: "gaussian", "gaussian_log", "lognormal", "student", "gamma".
#
# `model_on_log_scale` controls whether exp(estimate) is a valid HC/NC ratio:
# - TRUE: Gaussian on log(RC), log-normal, or Gamma-log model
# - FALSE: Gaussian or Student-t model with identity link
# -------------------------------------------------------------------------
final_model_id <- "student"

final_model <- switch(
  final_model_id,
  gaussian     = model_gaussian_RC,
  gaussian_log = model_gaussian_log_RC,
  lognormal    = model_logNormal_RC,
  student      = model_student_RC,
  gamma        = model_gamma_RC,
  stop(
    "Unknown final_model_id. Use one of: ",
    "gaussian, gaussian_log, lognormal, student, gamma."
  )
)

model_on_log_scale <- final_model_id %in% c(
  "gaussian_log",
  "lognormal",
  "gamma"
)

contrast_interpretation <- if (model_on_log_scale) {
  "estimate is on the log scale; exp(estimate) is the HC/NC ratio."
} else {
  "estimate is the additive HC - NC difference in relative-concentration units."
}

# Use the model/link scale explicitly. For log-scale models, this makes
# exponentiation of the HC - NC contrast a valid HC/NC ratio.
emm_RC <- emmeans(final_model, ~ diet | day, type = "link")
plot(emm_RC, comparisons = TRUE)
print(emm_RC)

# Reversed pairwise contrast gives HC - NC. BH adjustment is applied across days.
all_contrasts_RC <- contrast(
  emm_RC,
  method = "revpairwise",
  adjust = "none"
) %>%
  as.data.frame() %>%
  mutate(
    p_adj = p.adjust(p.value, method = "BH"),
    selected_model = final_model_id,
    contrast_scale = if_else(model_on_log_scale, "log scale", "identity scale"),
    effect_interpretation = if_else(
      model_on_log_scale,
      "HC/NC ratio after exponentiation",
      "HC - NC difference in RC units"
    ),
    # Valid only when the contrast is on a log scale.
    fold_change = if (model_on_log_scale) exp(estimate) else rep(NA_real_, dplyr::n()),
    sig = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

# Prepare position for asterisks above the boxplots
sig_labels_RC <- data_no_fasted %>%
  group_by(day) %>%
  summarise(y_pos = max(RC, na.rm = TRUE) + 0.005) %>%
  left_join(all_contrasts_RC %>% filter(contrast == "HC - NC") %>% select(day, sig), by = "day")

# Create bracket positions
# Filter for significant comparisons only
sig_brackets_RC_sig <- sig_labels_RC %>%
  filter(sig %in% c("*", "**", "***")) %>%
  mutate(
    x = as.numeric(day),
    xstart = x - 0.2,
    xend = x + 0.2,
    yend = y_pos + 0.1,
    ytext = yend + 0.05
  )

# Prepare data for faceted layout
day_levels <- sort(unique(data_no_fasted$day))

data_plot_RC <- data_no_fasted %>%
  mutate(
    day = factor(day, levels = day_levels),
    diet = factor(diet, levels = c("NC", "HC")),
    x_key = diet
  )

# Significance dataframe
sig_df_RC <- data_plot_RC %>%
  group_by(day) %>%
  summarise(
    y_position = max(RC, na.rm = TRUE) + 0.08 * diff(range(data_plot_RC$RC, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  left_join(
    all_contrasts_RC %>%
      filter(contrast == "HC - NC") %>%
      mutate(
        day = factor(day, levels = day_levels),
        annotations = sig
      ) %>%
      select(day, p_adj, annotations),
    by = "day"
  ) %>%
  mutate(
    xmin = "NC",
    xmax = "HC"
  ) %>%
  filter(!is.na(p_adj), p_adj < 0.05, annotations != "")

plot_title <- tools::file_path_sans_ext(basename(file_path))

lab_map <- c(
  "NC" = "NC",
  "HC" = "HC"
)
day_labs <- c(
  "1"  = "Day 1",
  "2"  = "Day 2",
  "3"  = "Day 3",
  "4"  = "Day 4",
  "10" = "Day 10",
  "15" = "Day 15",
  "22" = "Day 22"
)

plot <- ggplot(data_plot_RC, aes(x = x_key, y = RC)) +
  ggbeeswarm::geom_beeswarm(
    aes(color = diet),
    cex = 5,
    method = "compactswarm",
    priority = "density",
    size = 2,
    alpha = 0.8,
    corral = "gutter",
    corral.width = 0.75,
    preserve.data.axis = TRUE
  )  +
  scale_x_discrete(
    expand = expansion(add = c(2, 2))
  ) +
  facet_wrap(
    ~ day,
    nrow = 1,
    scales = "free_x",
    labeller = labeller(day = day_labs)
  ) +
  scale_color_manual(values = c(NC = "#44a7c4", HC = "#f99943")) +
  scale_x_discrete(labels = lab_map) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  labs(
    title = bquote(bold("Liver") ~ bolditalic(.(plot_title))),
    x = NULL,
    y = expression(
      italic("mRNA RUC") / 
        "(" * italic("18S rRNA") * ", " * italic("ef1a") * ")"
    ),
    color = "Diet"
  ) +
  {
    if (nrow(sig_df_RC) > 0) {
      ggsignif::geom_signif(
        data = sig_df_RC,
        aes(
          xmin = xmin,
          xmax = xmax,
          annotations = annotations,
          y_position = y_position
        ),
        manual = TRUE,
        inherit.aes = FALSE,
        tip_length = 0.02,
        textsize = 7,
        vjust = 0.4
      )
    } else {
      NULL
    }
  } +
  theme(
    plot.title = element_text(size = 20,
                              face = "italic",
                              hjust = 0.5,
                              family = "mono"),
    strip.text = element_text(size = 12),
    legend.position = "right",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16)
  )

plot


saveRDS(
  plot,
  file.path("RDS", paste0(plot_title, "_plot.rds"))
)


# Save plots
pdf_path  <- file.path(results_dir, paste0(file_base, "_facet_plot.pdf"))
tiff_path <- file.path(results_dir, paste0(file_base, "_facet_plot.tiff"))

ggsave(pdf_path, plot, width = 6, height = 4)
ggsave(tiff_path, plot, width = 6, height = 4, dpi = 1200, compression = "lzw")

# DHARMa residual diagnostic plot
res_final <- DHARMa::simulateResiduals(final_model, n = 1000, refit = FALSE)

pdf(
  file.path(results_dir, paste0(file_base, "_DHARMa_residuals.pdf")),
  width = 7,
  height = 7
)
plot(res_final, quantreg = FALSE)
dev.off()

tiff(
  file.path(results_dir, paste0(file_base, "_DHARMa_residuals.tiff")),
  width = 7,
  height = 7,
  units = "in",
  res = 300,
  compression = "lzw"
)
plot(res_final, quantreg = FALSE)
dev.off()

# Save model summary
sink(file.path(results_dir, paste0(file_base, "_model_summary.txt")))
cat("Model summary for:", file_base, "\n")
cat("Selected model:", final_model_id, "\n")
cat("Contrast interpretation:", contrast_interpretation, "\n\n")
print(summary(final_model))
sink()

# Save contrasts
write.csv(
  all_contrasts_RC,
  file.path(results_dir, paste0(file_base, "_contrasts.csv")),
  row.names = FALSE
)

# Save session information
writeLines(
  capture.output(sessionInfo()),
  file.path(results_dir, "sessionInfo.txt")
)






# plot data without statistical annotations (only for gcka and gckb)
data_plot_RC <- data_no_fasted %>%
  mutate(
    day = factor(day, levels = day_levels),
    diet = factor(diet, levels = c("NC", "HC")),
    x_key = diet
  )
plot_title <- tools::file_path_sans_ext(basename(file_path))
plot_no_stat <- ggplot(data_plot_RC, aes(x = x_key, y = RC)) +
  ggbeeswarm::geom_beeswarm(
    aes(color = diet),
    cex = 5,
    method = "compactswarm",
    priority = "density",
    size = 2,
    alpha = 0.8,
    corral = "gutter",
    corral.width = 0.75,
    preserve.data.axis = TRUE
  ) +
  scale_x_discrete(expand = expansion(add = c(2, 2))) +
  facet_wrap(
    ~ day,
    nrow = 1,
    scales = "free_x",
    labeller = labeller(day = day_labs)
  ) +
  scale_color_manual(values = c(NC = "#44a7c4", HC = "#f99943")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  labs(
    title = bquote(bold("Liver") ~ bolditalic(.(plot_title))),
    x = NULL,
    y = expression(
      italic("mRNA RUC") /
        "(" * italic("18S rRNA") * ", " * italic("ef1a") * ")"
    ),
    color = "Diet"
  ) +
  theme(
    plot.title = element_text(
      size = 20,
      face = "italic",
      hjust = 0.5,
      family = "mono"
    ),
    strip.text = element_text(size = 12),
    legend.position = "right",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16)
  )

plot_no_stat

saveRDS(
  plot_no_stat,
  file.path("RDS", paste0(plot_title, "_plot_no_stat.rds"))
)

ggsave(
  file.path(results_dir, paste0(file_base, "_facet_plot_no_stat.pdf")),
  plot_no_stat,
  width = 6,
  height = 4,
  dpi = 1200
)

ggsave(
  file.path(results_dir, paste0(file_base, "_facet_plot_no_stat.tiff")),
  plot_no_stat,
  width = 6,
  height = 4,
  dpi = 1200,
  compression = "lzw"
)

