# Libraries
library(readr)
library(dplyr)
library(ggplot2)
library(glmmTMB)
library(emmeans)
library(DHARMa)
library(performance)
library(insight)
library(purrr)
library(ggsignif)
library(ggbeeswarm)

# Load file
file_path <- "input_files/plasma_llactate.tsv"
# setting output directory
analysis_name <- tools::file_path_sans_ext(basename(file_path))
file_base <- analysis_name
results_dir <- file.path("results", analysis_name)
rds_dir <- "RDS"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

# Inspect file
data <- read_tsv(file_path)
str(data)

# Convert factors to the correct format
data$day <- as.factor(data$day)
data$diet <- as.factor(data$diet)

levels(data$diet)
levels(data$day)

# Convert to data frame and inspect structure
data <- as.data.frame(data)
str(data)

# Filter out Fasted group
data_no_fasted <- data %>%
  filter(diet != "Fasted") %>%
  droplevels()

str(data_no_fasted)

# Factor ordering
data_no_fasted$diet <- factor(data_no_fasted$diet, levels = c("NC","HC"))
data_no_fasted <- data_no_fasted %>%
  mutate(day_num = parse_number(as.character(day)))
# check
table(data_no_fasted$day, data_no_fasted$day_num)

# Ensure day is a factor and rename levels for clarity
data_no_fasted$day <- factor(
  data_no_fasted$day,
  levels = c("day_01", "day_02", "day_03", "day_04", "day_10", "day_15", "day_22"),
  labels = c("Day 1", "Day 2", "Day 3", "Day 4", "Day 10", "Day 15", "Day 22")
)


# Data Exploratory Visualization

ggplot(data_no_fasted, aes(x = Plasma_L_Lactate)) +
  # Histogram of Plasma L-Lactate (scaled to density)
  geom_histogram(aes(y = ..density..), bins = 30, fill = "steelblue", alpha = 0.5) +
  
  # Kernel density estimate overlay
  geom_density(color = "red", linewidth = 1.2) +
  
  # Vertical line for the mean
  geom_vline(aes(xintercept = mean(Plasma_L_Lactate, na.rm = TRUE)), 
             color = "darkgreen", linetype = "dashed", linewidth = 1) +
  
  # Vertical line for the median
  geom_vline(aes(xintercept = median(Plasma_L_Lactate, na.rm = TRUE)), 
             color = "orange", linetype = "dotted", linewidth = 1) +
  
  # Labels and theme
  labs(
    title = "Distribution of Plasma L-Lactate (mM)",
    subtitle = "Histogram with density curve, mean (dashed), and median (dotted)",
    x = "Plasma L-Lactate (mM)",
    y = "Density"
  ) +
  
  theme_minimal(base_size = 14)


# Compute per-diet statistics
stats_by_diet <- data_no_fasted %>%
  group_by(diet) %>%
  summarise(
    mean_lla   = mean(Plasma_L_Lactate, na.rm = TRUE),
    median_lla = median(Plasma_L_Lactate, na.rm = TRUE),
    logmean_lla = exp(mean(log(Plasma_L_Lactate[Plasma_L_Lactate > 0]), na.rm = TRUE)),
    .groups = "drop"
  )

# Plot with per-diet summary lines
ggplot(data_no_fasted, aes(x = Plasma_L_Lactate)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "steelblue", alpha = 0.5) +
  geom_density(color = "red", linewidth = 1.2) +
  # Facet-specific lines for each diet
  geom_vline(data = stats_by_diet, aes(xintercept = mean_lla, color = "Mean"), 
             linetype = "dashed", linewidth = 1) +
  geom_vline(data = stats_by_diet, aes(xintercept = median_lla, color = "Median"), 
             linetype = "dotted", linewidth = 1) +
  geom_vline(data = stats_by_diet, aes(xintercept = logmean_lla, color = "Log-mean"), 
             linetype = "dotdash", linewidth = 1) +
  geom_rug(alpha = 0.2) +
  facet_wrap(~ diet, scales = "free") +
  scale_color_manual(
    name = "Statistics",
    values = c("Mean" = "darkgreen", "Median" = "orange", "Log-mean" = "purple")
  ) +
  labs(
    title = "Distribution of Plasma L-Lactate (mM)",
    x = "Plasma L-Lactate (mM)",
    y = "Density"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")


# Fit candidate models

data_no_fasted <- data_no_fasted %>%
  mutate(
    diet = factor(diet, levels = c("NC", "HC"))
  )

# Optional check for positive values, needed for log and Gamma models
if (any(data_no_fasted$Plasma_L_Lactate <= 0, na.rm = TRUE)) {
  stop("Plasma_L_Lactate contains zero or negative values. Log-Gaussian and Gamma models require positive values.")
}

model_gaussian <- glmmTMB(
  Plasma_L_Lactate ~ day * diet,
  data = data_no_fasted,
  family = gaussian(link = "identity")
)

model_gaussian_disp <- glmmTMB(
  Plasma_L_Lactate ~ day * diet,
  dispformula = ~ day * diet,
  data = data_no_fasted,
  family = gaussian(link = "identity")
)

model_gaussian_log <- glmmTMB(
  log(Plasma_L_Lactate) ~ day * diet,
  data = data_no_fasted,
  family = gaussian(link = "identity")
)

model_gaussian_log_disp <- glmmTMB(
  log(Plasma_L_Lactate) ~ day * diet,
  dispformula = ~ day * diet,
  data = data_no_fasted,
  family = gaussian(link = "identity")
)

model_student <- glmmTMB(
  Plasma_L_Lactate ~ day * diet,
  data = data_no_fasted,
  family = t_family(link = "identity")
)

model_student_disp <- glmmTMB(
  Plasma_L_Lactate ~ day * diet,
  dispformula = ~ day * diet,
  data = data_no_fasted,
  family = t_family(link = "identity")
)

model_gamma <- glmmTMB(
  Plasma_L_Lactate ~ day * diet,
  data = data_no_fasted,
  family = Gamma(link = "log")
)

model_gamma_disp <- glmmTMB(
  Plasma_L_Lactate ~ day * diet,
  dispformula = ~ day * diet,
  data = data_no_fasted,
  family = Gamma(link = "log")
)


# Model comparison

models <- list(
  gaussian = model_gaussian,
  gaussian_disp = model_gaussian_disp,
  gaussian_log = model_gaussian_log,
  gaussian_log_disp = model_gaussian_log_disp,
  student = model_student,
  student_disp = model_student_disp,
  gamma = model_gamma,
  gamma_disp = model_gamma_disp
)

model_summaries <- purrr::map_df(names(models), function(name) {
  m <- models[[name]]
  fam <- m$modelInfo$family
  
  data.frame(
    model = name,
    family = fam$family,
    link = fam$link,
    formula = deparse(formula(m)),
    dispformula = deparse(m$modelInfo$allForm$dispformula),
    AIC = AIC(m),
    BIC = BIC(m),
    logLik = as.numeric(logLik(m)),
    df_resid = df.residual(m),
    convergence = m$fit$convergence,
    message = m$fit$message,
    pdHess = m$sdr$pdHess,
    stringsAsFactors = FALSE
  )
})

print(model_summaries[order(model_summaries$AIC), ])

# Compare nested models within the same family / same response scale
anova(model_gaussian, model_gaussian_disp)
anova(model_gaussian_log, model_gaussian_log_disp)
anova(model_student, model_student_disp)
anova(model_gamma, model_gamma_disp)


# DHARMa diagnostics for candidate models

check_dharma <- function(model, model_name) {
  cat("\n=============================\n")
  cat(model_name, "\n")
  cat("=============================\n")
  
  res <- DHARMa::simulateResiduals(
    model,
    n = 1000,
    refit = FALSE
  )
  
  plot(res, quantreg = FALSE)
  print(DHARMa::testUniformity(res))
  print(DHARMa::testDispersion(res))
  print(DHARMa::testOutliers(res))
  
  invisible(res)
}

res_gaussian          <- check_dharma(model_gaussian, "Gaussian")
res_gaussian_disp     <- check_dharma(model_gaussian_disp, "Gaussian + dispersion")
res_gaussian_log      <- check_dharma(model_gaussian_log, "Log Gaussian")
res_gaussian_log_disp <- check_dharma(model_gaussian_log_disp, "Log Gaussian + dispersion")
res_student           <- check_dharma(model_student, "Student")
res_student_disp      <- check_dharma(model_student_disp, "Student + dispersion")
res_gamma             <- check_dharma(model_gamma, "Gamma")
res_gamma_disp        <- check_dharma(model_gamma_disp, "Gamma + dispersion")


# Final model

final_model <- model_gamma_disp

summary(final_model)

performance::check_convergence(final_model)
insight::is_converged(final_model)

final_model$fit$convergence
final_model$fit$message
final_model$sdr$pdHess

res_final <- DHARMa::simulateResiduals(
  final_model,
  n = 1000,
  refit = FALSE
)

plot(res_final, quantreg = FALSE)
DHARMa::testUniformity(res_final)
DHARMa::testDispersion(res_final)
DHARMa::testOutliers(res_final)


# Pairwise contrasts

emm <- emmeans(final_model, ~ diet | day)

print(summary(emm, type = "response"))

all_contrasts <- contrast(
  emm,
  method = "revpairwise",
  adjust = "none"
) %>%
  as.data.frame() %>%
  mutate(
    p_adj = p.adjust(p.value, method = "BH"),
    fold_change = exp(estimate),
    sig = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

all_contrasts$day <- factor(
  as.character(all_contrasts$day),
  levels = levels(data_no_fasted$day)
)


# Plot data

data_plot <- data_no_fasted %>%
  mutate(
    diet = factor(diet, levels = c("NC", "HC")),
    x_pos = case_when(
      diet == "NC" ~ 1,
      diet == "HC" ~ 3
    )
  )

yrng <- range(data_plot$Plasma_L_Lactate, na.rm = TRUE)
pad  <- diff(yrng) * 0.08

sig_df <- data_plot %>%
  group_by(day) %>%
  summarise(
    y_position = max(Plasma_L_Lactate, na.rm = TRUE) + pad,
    .groups = "drop"
  ) %>%
  left_join(
    all_contrasts %>%
      filter(contrast %in% c("HC - NC", "NC - HC")) %>%
      select(day, p_adj, sig),
    by = "day"
  ) %>%
  filter(!is.na(p_adj), p_adj < 0.05, sig != "") %>%
  mutate(
    xmin = 1,
    xmax = 3,
    annotations = sig
  )

# 10. Final Plasma L-Lactate plot

p_L_La <- ggplot(data_plot, aes(x = x_pos, y = Plasma_L_Lactate)) +
  ggbeeswarm::geom_beeswarm(
    aes(color = diet, group = diet),
    method = "swarm",
    priority = "density",
    cex = 4,
    size = 1.5,
    alpha = 0.55,
    corral = "none",
    preserve.data.axis = TRUE
  ) +
  facet_wrap(
    ~ day,
    nrow = 1,
    scales = "free_x"
  ) +
  scale_color_manual(values = c(NC = "#44a7c4", HC = "#f99943")) +
  scale_x_continuous(
    limits = c(0, 4),
    breaks = c(1, 3),
    labels = c("", "")
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.18))
  ) +
  labs(
    title = "Plasma L-Lactate",
    x = NULL,
    y = "mM",
    color = "Diet"
  ) +
  {
    if (nrow(sig_df) > 0) {
      ggsignif::geom_signif(
        data = sig_df,
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
    plot.title = element_text(
      size = 20,
      face = "plain",
      hjust = 0.5,
      family = "mono"
    ),
    strip.text = element_text(size = 12),
    legend.position = "right",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16)
  )

p_L_La
saveRDS(
  p_L_La,
  file.path(rds_dir, paste0(analysis_name, "_plot.rds"))
)



ggsave(
  file.path(results_dir, paste0(file_base, "_facet_plot.pdf")),
  p_L_La,
  width = 6,
  height = 3.5,
  dpi = 1200
)

ggsave(
  file.path(results_dir, paste0(file_base, "_facet_plot.tiff")),
  p_L_La,
  width = 6,
  height = 3.5,
  dpi = 1200
)

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
  res = 300
)
plot(res_final, quantreg = FALSE)
dev.off()

sink(file.path(results_dir, paste0(file_base, "_model_summary.txt")))
cat("Final model for:", file_base, "\n\n")
print(summary(final_model))

cat("\n\nConvergence checks:\n")
cat("performance::check_convergence:\n")
print(performance::check_convergence(final_model))

cat("\ninsight::is_converged:\n")
print(insight::is_converged(final_model))

cat("\nfit$convergence:\n")
print(final_model$fit$convergence)

cat("\nfit$message:\n")
print(final_model$fit$message)

cat("\nsdr$pdHess:\n")
print(final_model$sdr$pdHess)

cat("\n\nDHARMa diagnostics:\n")
print(DHARMa::testUniformity(res_final))
print(DHARMa::testDispersion(res_final))
print(DHARMa::testOutliers(res_final))
sink()

write.csv(
  model_summaries,
  file.path(results_dir, paste0(file_base, "_model_comparison.csv")),
  row.names = FALSE
)

write.csv(
  all_contrasts,
  file.path(results_dir, paste0(file_base, "_contrasts.csv")),
  row.names = FALSE
)


writeLines(
  capture.output(sessionInfo()),
  file.path(results_dir, "sessionInfo.txt")
)
