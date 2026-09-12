# Load packages
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(brms)
  library(posterior)
  library(broom.mixed)
  library(loo)
  library(scales)
})

# Define input file and heatmap label type
infile <- "input_files/methylome_for_dirichlet.tsv"
label_mode <- "centered_prob"

# Load methylome data
df <- readr::read_tsv(infile)

# Remove the Fasted group
df <- df %>%
  filter(diet != "Fasted") %>%
  droplevels()

# Calculate descriptive statistics by diet and day
summary_tbl <- df %>%
  group_by(diet, day) %>%
  summarise(
    n = n(),
    mean_dC = mean(dC_rel, na.rm = TRUE),
    sd_dC = sd(dC_rel, na.rm = TRUE),
    cv_dC = 100 * sd_dC / mean_dC,
    mean_mdC = mean(mdC_rel, na.rm = TRUE),
    sd_mdC = sd(mdC_rel, na.rm = TRUE),
    cv_mdC = 100 * sd_mdC / mean_mdC,
    mean_hmdC = mean(hmdC_rel, na.rm = TRUE),
    sd_hmdC = sd(hmdC_rel, na.rm = TRUE),
    cv_hmdC = 100 * sd_hmdC / mean_hmdC,
    .groups = "drop"
  ) %>%
  arrange(diet, day)

print(summary_tbl)

# Save all sensitivity-analysis outputs in a separate directory
out_dir <- file.path("results", "dirichlet_regression_weak_priors")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sample_size_dir <- file.path(out_dir, "sample_size")
dir.create(sample_size_dir, recursive = TRUE, showWarnings = FALSE)

readr::write_csv(
  summary_tbl,
  file.path(sample_size_dir, "methylome_summary_by_diet_day.csv")
)

# This sensitivity analysis does not use the observed NC Day 1 composition
# to centre the intercept priors.

# Check that all required columns are available
need <- c("diet", "day", "dC_rel", "mdC_rel", "hmdC_rel")
miss <- setdiff(need, names(df))

if (length(miss) > 0) {
  stop("Missing columns: ", paste(miss, collapse = ", "))
}

# Convert grouping variables to factors and remove missing values
df <- df %>%
  mutate(
    diet = factor(diet),
    day = factor(day)
  ) %>%
  filter(if_all(all_of(need), ~ !is.na(.))) %>%
  droplevels()

# Check that methylation fractions form valid compositions
comp <- as.matrix(df %>% select(dC_rel, mdC_rel, hmdC_rel))
rs <- rowSums(comp)

stopifnot(
  all(is.finite(comp)),
  all(comp > 0),
  all(comp < 1),
  all(abs(rs - 1) < 1e-8)
)

# Set NC as the diet reference level
if (all(c("NC", "HC") %in% levels(df$diet))) {
  df$diet <- factor(df$diet, levels = c("NC", "HC"))
}

# Define the Dirichlet regression model
form_dir <- bf(
  cbind(dC_rel, mdC_rel, hmdC_rel) ~ diet * day
)

# Create default priors
prior_tbl <- brms::default_prior(
  form_dir,
  data = df,
  family = dirichlet()
)

# Set conventional weakly informative, non-empirical intercept priors
# on the log-ratio scale relative to the reference cytosine component.
prior_tbl$prior[
  prior_tbl$class == "Intercept" &
    prior_tbl$dpar == "mumdCrel"
] <- "normal(0, 2.5)"

prior_tbl$prior[
  prior_tbl$class == "Intercept" &
    prior_tbl$dpar == "muhmdCrel"
] <- "normal(0, 2.5)"

# Set priors for diet, day, and interaction effects
prior_tbl$prior[prior_tbl$class == "b"] <- "normal(0, 0.5)"

# Set the Dirichlet precision prior
if (any(prior_tbl$class == "phi")) {
  prior_tbl$prior[prior_tbl$class == "phi"] <- "exponential(0.25)"
}

# Save the exact prior specification used in this sensitivity analysis
readr::write_csv(
  prior_tbl,
  file.path(out_dir, "prior_specification_weak_nonempirical.csv")
)

# Fit the Bayesian Dirichlet regression model
fit_dir <- brm(
  form_dir,
  data = df,
  family = dirichlet(),
  prior = prior_tbl,
  chains = 4,
  iter = 3000,
  warmup = 1000,
  seed = 123,
  init = 0,
  control = list(
    adapt_delta = 0.9995,
    max_treedepth = 15
  ),
  save_pars = save_pars(all = TRUE)
)

fit_best <- fit_dir

# Save the fitted model
saveRDS(
  fit_best,
  file.path(out_dir, "fit_dirichlet_diet_day.rds")
)

# Summarise posterior draws
summarize_draws <- function(x) {
  tibble(
    estimate = mean(x),
    lower95 = quantile(x, 0.025),
    upper95 = quantile(x, 0.975),
    prob_pos = mean(x > 0)
  )
}

# Calculate HC minus NC contrasts for each day and component
day_levels <- levels(df$day)
components <- c("dC_rel", "mdC_rel", "hmdC_rel")
contrast_list <- list()

for (dy in day_levels) {
  
  nd_HC <- data.frame(
    diet = factor("HC", levels = levels(df$diet)),
    day = factor(dy, levels = levels(df$day))
  )
  
  nd_NC <- data.frame(
    diet = factor("NC", levels = levels(df$diet)),
    day = factor(dy, levels = levels(df$day))
  )
  
  mu_HC_all <- posterior_epred(
    fit_best,
    newdata = nd_HC,
    re_formula = NA
  )
  
  mu_NC_all <- posterior_epred(
    fit_best,
    newdata = nd_NC,
    re_formula = NA
  )
  
  for (j in seq_along(components)) {
    
    delta <- mu_HC_all[, 1, j] - mu_NC_all[, 1, j]
    
    contrast_list[[length(contrast_list) + 1]] <-
      summarize_draws(delta) %>%
      mutate(
        response = components[j],
        day = dy,
        abs_est = abs(estimate),
        direction = case_when(
          estimate > 0 ~ "increase",
          estimate < 0 ~ "decrease",
          TRUE ~ "neutral"
        ),
        sig95 = lower95 > 0 | upper95 < 0
      )
  }
}

contr_summ <- bind_rows(contrast_list) %>%
  mutate(day = factor(day, levels = day_levels))

# Calculate percentage-point effects and posterior certainty
contr_summ <- contr_summ %>%
  mutate(
    estimate_pp = 100 * estimate,
    lower95_pp = 100 * lower95,
    upper95_pp = 100 * upper95,
    abs_est_pp = 100 * abs_est,
    abs_est_clean = ifelse(abs(estimate) < 1e-6, 0, abs(estimate)),
    abs_est_pp_clean = ifelse(abs_est_pp < 1e-4, 0, abs_est_pp),
    prob_certainty = ifelse(prob_pos < 0.5, 1 - prob_pos, prob_pos),
    label_text = case_when(
      label_mode == "centered_prob" ~ sprintf("%.2f", prob_certainty),
      label_mode == "delta_pp" ~ sprintf("%+.2f pp", estimate_pp),
      TRUE ~ sprintf("%.2f", prob_certainty)
    )
  )

# Save HC minus NC contrasts
readr::write_csv(
  contr_summ,
  file.path(out_dir, "contrasts_HC_minus_NC_by_day_response.csv")
)

# Create all diet and day prediction combinations
newdata <- expand.grid(
  diet = levels(df$diet),
  day = levels(df$day),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

newdata$diet <- factor(newdata$diet, levels = levels(df$diet))
newdata$day <- factor(newdata$day, levels = levels(df$day))

# Calculate posterior predictions
mu_all <- posterior_epred(
  fit_best,
  newdata = newdata,
  re_formula = NA
)

stopifnot(length(dim(mu_all)) == 3)

# Summarise posterior predictions for each cytosine form
comp_names <- c("dC_rel", "mdC_rel", "hmdC_rel")

pred_summary <- lapply(seq_along(comp_names), function(j) {
  
  mat <- mu_all[, , j]
  
  tibble(
    cell = seq_len(ncol(mat)),
    mean = colMeans(mat),
    lower = apply(mat, 2, quantile, 0.025),
    upper = apply(mat, 2, quantile, 0.975),
    response = comp_names[j]
  )
  
}) %>%
  bind_rows() %>%
  mutate(
    diet = newdata$diet[cell],
    day = newdata$day[cell]
  ) %>%
  select(day, diet, response, mean, lower, upper)

# Convert observed methylation fractions to long format
df_long <- df %>%
  select(diet, day, dC_rel, mdC_rel, hmdC_rel) %>%
  pivot_longer(
    c(dC_rel, mdC_rel, hmdC_rel),
    names_to = "response",
    values_to = "frac"
  ) %>%
  mutate(response = factor(response, levels = comp_names))

# Define positions for observed and predicted values
pd <- position_dodge(width = 0.6)
pj <- position_jitterdodge(
  jitter.width = 0.15,
  dodge.width = 0.6
)

# Plot observed fractions and posterior credible intervals
p_jitter_ci <- ggplot() +
  geom_point(
    data = df_long,
    aes(x = day, y = frac, color = diet),
    position = pj,
    alpha = 0.35,
    size = 1.4
  ) +
  geom_pointrange(
    data = pred_summary,
    aes(
      x = day,
      y = mean,
      ymin = lower,
      ymax = upper,
      color = diet
    ),
    position = pd,
    size = 0.6,
    fatten = 1.4
  ) +
  facet_wrap(~response, ncol = 1, scales = "free_y") +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 0.1)
  ) +
  labs(
    x = "Day",
    y = "Relative fraction",
    title = "Raw cytosine fractions with posterior mean ± 95% CrI",
    subtitle = "Dirichlet model: diet × day"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

# Save the posterior interval figure
ggsave(
  file.path(out_dir, "jitter_with_posterior_CI_dirichlet.pdf"),
  p_jitter_ci,
  width = 8,
  height = 8
)

# Summarise the fitted model
sum_dir <- summary(fit_best)
print(sum_dir)

# Save the model summary
writeLines(
  capture.output(sum_dir),
  file.path(out_dir, "diagnostics_summary.txt")
)

# Calculate leave-one-out cross-validation
loo_dir <- loo(fit_best)
print(loo_dir)

# Save leave-one-out results
saveRDS(
  loo_dir,
  file.path(out_dir, "loo_dirichlet.rds")
)

writeLines(
  capture.output(loo_dir),
  file.path(out_dir, "loo_dirichlet.txt")
)

# Calculate fitted values for the observed samples
mu_all_obs <- posterior_epred(
  fit_best,
  re_formula = NA
)

stopifnot(dim(mu_all_obs)[2] == nrow(df))

# Create fitted versus observed data
fit_obs <- bind_rows(
  tibble(
    response = "dC_rel",
    y_obs = df$dC_rel,
    y_hat = colMeans(mu_all_obs[, , 1])
  ),
  tibble(
    response = "mdC_rel",
    y_obs = df$mdC_rel,
    y_hat = colMeans(mu_all_obs[, , 2])
  ),
  tibble(
    response = "hmdC_rel",
    y_obs = df$hmdC_rel,
    y_hat = colMeans(mu_all_obs[, , 3])
  )
)

# Plot fitted versus observed fractions
p_fitobs <- ggplot(fit_obs, aes(x = y_hat, y = y_obs)) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_point(alpha = 0.6) +
  facet_wrap(~response, scales = "free") +
  labs(
    x = "Posterior mean (fitted)",
    y = "Observed",
    title = "Fitted vs observed (Dirichlet: diet × day)"
  ) +
  theme_bw(base_size = 12)

# Save the fitted versus observed figure
ggsave(
  file.path(out_dir, "fitted_vs_observed_dirichlet.pdf"),
  p_fitobs,
  width = 7.5,
  height = 5,
  dpi = 600
)

# Save fixed model effects
fixed_effects <- as.data.frame(brms::fixef(fit_best, summary = TRUE))

fixed_effects$parameter <- rownames(fixed_effects)
rownames(fixed_effects) <- NULL

readr::write_csv(
  fixed_effects,
  file.path(out_dir, "brms_fixed_effects_dirichlet.csv")
)

# Save R session information
writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "sessionInfo.txt")
)

# Print the main output location
message("Results saved in: ", out_dir)