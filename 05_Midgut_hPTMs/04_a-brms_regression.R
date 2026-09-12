suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(brms)
  library(loo)
  library(posterior)
  library(tibble)
  library(emmeans)
  library(tidybayes)
  library(ggplot2)
  library(bayesplot)
  library(readr)
})

theme_set(theme_bw())
available_cores <- parallel::detectCores(logical = TRUE)

options(mc.cores = max(1, available_cores - 1))

cat("Detected logical cores:", available_cores, "\n")
cat("R allowed cores:", getOption("mc.cores"), "\n")

# Short writable temporary directory
tmp_dir <- file.path(getwd(), "tmp_brms")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

Sys.setenv(TMPDIR = normalizePath(tmp_dir))
options(tmpdir = normalizePath(tmp_dir))

cat("Using TMPDIR:", Sys.getenv("TMPDIR"), "\n")
cat("tempdir():", tempdir(), "\n")

# Load processed matrix and metadata
histone_data <- readRDS(
  "02_process_matrix/H3_H4_5PTM_midgut.rds"
) %>%
  as.data.frame()

histone_metadata <- data.table::fread(
  "STPN2309_metadata_all.tsv"
) %>%
  as.data.frame()

# First metadata column contains the sample IDs
rownames(histone_metadata) <- histone_metadata[[1]]
histone_metadata <- histone_metadata[, -1, drop = FALSE]



# Settings
# ---------------------------------------------------------
# [2] Settings and output directory
# ---------------------------------------------------------
BATCH_FIXED_COL <- "gut_hptm_batch"

DAY_LEVELS <- c(
  "day_01", "day_02", "day_03", "day_04",
  "day_10", "day_15", "day_22"
)

output_root <- "04_brms_regression"
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

output_prefix <- file.path(
  output_root,
  "midgut_hPTM_batch_adjusted"
)



# Align sample IDs between PTM matrix and metadata
ids_ptm  <- rownames(histone_data)
ids_meta <- rownames(histone_metadata)

common_ids <- sort(intersect(ids_ptm, ids_meta))

if (length(common_ids) == 0) {
  stop("No shared sample IDs between the hPTM matrix and metadata.")
}

histone_data <- histone_data[common_ids, , drop = FALSE]

histone_metadata <- histone_metadata[
  common_ids,
  ,
  drop = FALSE
]

stopifnot(
  identical(rownames(histone_data), rownames(histone_metadata))
)


# Keep NC and HC samples only
histone_metadata <- histone_metadata %>%
  filter(diet %in% c("NC", "HC"))

if (nrow(histone_metadata) == 0) {
  stop("No NC/HC samples remain after metadata filtering.")
}

histone_data <- histone_data[
  rownames(histone_metadata),
  ,
  drop = FALSE
]

stopifnot(
  identical(rownames(histone_data), rownames(histone_metadata))
)


# Metadata preparation
required_metadata <- c(
  "diet",
  "day_categ",
  BATCH_FIXED_COL
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(histone_metadata)
)

if (length(missing_metadata) > 0) {
  stop(
    "Metadata is missing: ",
    paste(missing_metadata, collapse = ", ")
  )
}

histone_metadata <- histone_metadata %>%
  mutate(
    diet = factor(diet, levels = c("NC", "HC")),
    day_categ = factor(day_categ, levels = DAY_LEVELS),
    gut_hptm_batch = factor(.data[[BATCH_FIXED_COL]])
  )

histone_metadata$diet <- droplevels(histone_metadata$diet)
histone_metadata$day_categ <- droplevels(histone_metadata$day_categ)
histone_metadata$gut_hptm_batch <- droplevels(
  histone_metadata$gut_hptm_batch
)

if (any(is.na(histone_metadata$day_categ))) {
  stop("Unexpected values detected in day_categ.")
}

if (nlevels(histone_metadata$gut_hptm_batch) < 2) {
  stop("gut_hptm_batch has fewer than two levels.")
}

cat("\nSamples retained:", nrow(histone_metadata), "\n")

cat("\nDiet counts:\n")
print(table(histone_metadata$diet))

cat("\nDiet × day table:\n")
print(table(histone_metadata$diet, histone_metadata$day_categ))

cat("\nBatch counts:\n")
print(table(histone_metadata$gut_hptm_batch))

histone_features <- colnames(histone_data)

# Candidate families
candidate_fams <- list(
  gaussian    = gaussian(),
  student     = student(),
  skew_normal = skew_normal()
)

# Priors for family comparison
priors_common <- c(
  prior(normal(0, 2),       class = "b"),
  prior(student_t(3, 0, 5), class = "Intercept"),
  prior(exponential(1),     class = "sigma")
)

get_priors_for_family <- function(family_name) {
  family_name <- tolower(family_name)
  
  if (family_name == "student") {
    return(
      c(
        priors_common,
        prior(gamma(2, 0.2), class = "nu")
      )
    )
  }
  
  if (family_name %in% c("gaussian", "skew_normal")) {
    return(priors_common)
  }
  
  stop("Unsupported family: ", family_name)
}

get_brms_family <- function(family_name) {
  family_name <- tolower(family_name)
  
  if (family_name == "student") return(student())
  if (family_name == "skew_normal") return(skew_normal())
  if (family_name == "gaussian") return(gaussian())
  
  stop("Unsupported family: ", family_name)
}

# [8] Function to fit one feature × one family
fit_brms_one <- function(feature,
                         fam_name,
                         fam_obj,
                         feature_df,
                         meta_df,
                         priors = get_priors_for_family(fam_name),
                         iter     = 1000,
                         warmup   = 500,
                         chains   = 4,
                         cores    = 4,
                         save_dir = paste0(output_prefix, "_fits")) {
  
  stopifnot(identical(rownames(feature_df), rownames(meta_df)))
  
  dat_feat <- meta_df %>%
    mutate(response = feature_df[[feature]]) %>%
    filter(!is.na(response))
  
  n_after <- nrow(dat_feat)
  
  if (n_after < 6) {
    warning("Too few non-missing PTM values (n = ", n_after,
            ") for feature ", feature, " – skipping.")
    return(NULL)
  }
  
  dat_feat <- dat_feat %>%
    mutate(
      diet           = factor(diet, levels = c("NC", "HC")),
      day_categ      = factor(day_categ),
      gut_hptm_batch = factor(gut_hptm_batch)
    )
  
  if (nlevels(dat_feat$gut_hptm_batch) < 2) {
    warning("Feature ", feature, ": gut_hptm_batch has fewer than 2 levels after NA filtering – skipping.")
    return(NULL)
  }
  
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  safe_feat <- gsub("[^A-Za-z0-9]+", "_", feature)
  fit_file  <- file.path(save_dir, paste0("brms_", safe_feat, "_", fam_name))
  
  fit <- brm(
    formula    = bf(response ~ day_categ * diet + gut_hptm_batch),
    data       = dat_feat,
    family     = fam_obj,
    prior      = priors,
    iter       = iter,
    warmup     = warmup,
    chains     = chains,
    cores      = cores,
    file       = fit_file,
    file_refit = "on_change",
    silent     = 0,
    refresh    = 100
  )
  
  loo_fit  <- loo(fit, moment_match = TRUE)
  R2_vals  <- as.numeric(bayes_R2(fit))
  rhat_vec <- brms::rhat(fit)
  neff_rat <- brms::neff_ratio(fit)
  
  diag_row <- tibble(
    feature      = feature,
    family       = fam_name,
    formula_used = "response ~ day_categ * diet + gut_hptm_batch",
    looic        = loo_fit$estimates["looic", "Estimate"],
    looic_se     = loo_fit$estimates["looic", "SE"],
    elpd_loo     = loo_fit$estimates["elpd_loo", "Estimate"],
    elpd_loo_se  = loo_fit$estimates["elpd_loo", "SE"],
    R2_mean      = mean(R2_vals, na.rm = TRUE),
    R2_sd        = sd(R2_vals, na.rm = TRUE),
    Rhat_max     = max(rhat_vec, na.rm = TRUE),
    n_eff_min    = min(neff_rat, na.rm = TRUE),
    n_obs        = n_after
  )
  
  list(
    fit  = fit,
    loo  = loo_fit,
    diag = diag_row
  )
}

# Loop over all features and families
all_diag_list <- list()
all_fit_list  <- list()
error_log     <- list()

for (feat in histone_features) {
  message("=== Feature: ", feat, " ===")
  
  fam_results <- list()
  fam_diags   <- list()
  error_log[[feat]] <- list()
  
  for (fam_name in names(candidate_fams)) {
    message("  Family: ", fam_name)
    
    res <- tryCatch(
      fit_brms_one(
        feature    = feat,
        fam_name   = fam_name,
        fam_obj    = candidate_fams[[fam_name]],
        feature_df = histone_data,
        meta_df    = histone_metadata
      ),
      error = function(e) {
        error_log[[feat]][[fam_name]] <<- conditionMessage(e)
        warning("FAILED for feature ", feat, " family ", fam_name,
                " | error: ", conditionMessage(e))
        return(NULL)
      }
    )
    
    if (is.null(res)) next
    
    fam_results[[fam_name]] <- res$fit
    fam_diags[[fam_name]]   <- res$diag
  }
  
  if (length(fam_diags) == 0) next
  
  all_diag_list[[feat]] <- bind_rows(fam_diags)
  all_fit_list[[feat]]  <- fam_results
}

all_diagnostics <- bind_rows(all_diag_list)

if (nrow(all_diagnostics) == 0) {
  stop("No successful models were fitted. Check 'error_log'.")
}

write.csv(
  all_diagnostics,
  paste0(output_prefix, "_family_diagnostics.csv"),
  row.names = FALSE
)

saveRDS(error_log, paste0(output_prefix, "_family_error_log.rds"))

print(all_diagnostics)

cat("\n===== Diagnostic summaries =====\n")
print(summary(all_diagnostics$Rhat_max))
print(summary(all_diagnostics$n_eff_min))

# Select best family per feature
family_priority <- c("student", "gaussian", "skew_normal")

best_models <- all_diagnostics %>%
  group_by(feature) %>%
  arrange(looic, .by_group = TRUE) %>%
  mutate(
    best_looic    = min(looic, na.rm = TRUE),
    delta_looic   = looic - best_looic,
    best_looic_se = looic_se[which.min(looic)][1]
  ) %>%
  mutate(
    equivalent  = delta_looic < 2 * best_looic_se,
    family_rank = match(family, family_priority)
  ) %>%
  filter(equivalent) %>%
  arrange(family_rank, looic, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

write.csv(
  best_models,
  paste0(output_prefix, "_best_family_per_feature.csv"),
  row.names = FALSE
)

print(best_models)

# FINAL MODEL PER FEATURE USING THE BEST FAMILY SELECTED
# WITH gut_hptm_batch as fixed effect

# Read best-family table
model_selection_file <- paste0(output_prefix, "_best_family_per_feature.csv")

model_selection <- read_csv(model_selection_file, show_col_types = FALSE) %>%
  as.data.frame()

required_cols <- c("feature", "family")
missing_cols <- setdiff(required_cols, colnames(model_selection))
if (length(missing_cols) > 0) {
  stop("Model selection table is missing required columns: ",
       paste(missing_cols, collapse = ", "))
}

model_selection <- model_selection %>%
  distinct(feature, .keep_all = TRUE)

family_lookup <- setNames(model_selection$family, model_selection$feature)

# Output directories
run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
base_dir <- file.path(
  output_root,
  paste0("brms_best_family_final_batch_fixed_", run_id)
)
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)

writeLines(
  capture.output(sessionInfo()),
  file.path(base_dir, "session_info.txt")
)

summary_dir <- file.path(base_dir, "summary_tables")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  model_selection,
  file = file.path(base_dir, "model_selection_used.csv"),
  row.names = FALSE
)

# Function to fit final model for one feature
fit_best_final_one <- function(feature,
                               feature_df   = histone_data,
                               meta_df      = histone_metadata,
                               family_lookup,
                               iter         = 4000,
                               warmup       = 2000,
                               chains       = 4,
                               cores        = 4,
                               base_dir     = base_dir) {
  message("=== Final best-family model for feature: ", feature, " ===")
  
  stopifnot(identical(rownames(feature_df), rownames(meta_df)))
  
  if (!feature %in% names(family_lookup)) {
    warning("No selected family found for feature ", feature, " – skipping.")
    return(NULL)
  }
  
  selected_family_name <- family_lookup[[feature]]
  brms_family_obj      <- get_brms_family(selected_family_name)
  priors_this          <- get_priors_for_family(selected_family_name)
  
  dat_feat <- meta_df %>%
    mutate(response = feature_df[[feature]]) %>%
    filter(!is.na(response))
  
  n_after <- nrow(dat_feat)
  
  if (n_after < 6) {
    warning("Too few non-missing values (n = ", n_after,
            ") for feature ", feature, " – skipping.")
    return(NULL)
  }
  
  dat_feat <- dat_feat %>%
    mutate(
      diet           = factor(diet, levels = c("NC", "HC")),
      day_categ      = factor(day_categ),
      gut_hptm_batch = factor(gut_hptm_batch)
    )
  
  if (nlevels(dat_feat$gut_hptm_batch) < 2) {
    warning("Feature ", feature, ": gut_hptm_batch has fewer than 2 levels after NA filtering – skipping.")
    return(NULL)
  }
  
  safe_feat <- gsub("[^A-Za-z0-9]+", "_", feature)
  feat_dir  <- file.path(base_dir, safe_feat)
  dir.create(feat_dir, recursive = TRUE, showWarnings = FALSE)
  
  fit_file <- file.path(feat_dir, paste0("brms_", selected_family_name, "_", safe_feat))
  
  writeLines(
    c(
      paste0("feature: ", feature),
      paste0("family: ", selected_family_name),
      "formula: response ~ day_categ * diet + gut_hptm_batch"
    ),
    con = file.path(feat_dir, paste0(safe_feat, "_model_info.txt"))
  )
  
  fit <- brm(
    formula    = bf(response ~ day_categ * diet + gut_hptm_batch),
    data       = dat_feat,
    family     = brms_family_obj,
    prior      = priors_this,
    iter       = iter,
    warmup     = warmup,
    chains     = chains,
    cores      = cores,
    file       = fit_file,
    file_refit = "on_change",
    silent     = 0,
    refresh    = 100,
    control    = list(adapt_delta = 0.95, max_treedepth = 12)
  )
  
  sink(file.path(feat_dir, paste0(safe_feat, "_summary.txt")))
  cat("Feature: ", feature, "\n", sep = "")
  cat("Selected family: ", selected_family_name, "\n\n", sep = "")
  print(summary(fit))
  sink()
  
  p_pp <- pp_check(fit, ndraws = 100) +
    ggtitle(paste(feature, "-", selected_family_name, "- posterior predictive check"))
  
  ggsave(
    filename = file.path(feat_dir, paste0(safe_feat, "_ppcheck.png")),
    plot     = p_pp,
    width    = 6, height = 4, dpi = 1200
  )
  
  p_trace <- mcmc_plot(fit, type = "trace", regex_pars = "^b_")
  ggsave(
    filename = file.path(feat_dir, paste0(safe_feat, "_trace_fixed.png")),
    plot     = p_trace,
    width    = 8, height = 6, dpi = 1200
  )
  
  ce_list <- conditional_effects(fit, "day_categ:diet", re_formula = NA)
  p_ce <- plot(ce_list, plot = FALSE)[[1]] +
    ggtitle(paste(feature, "-", selected_family_name, "- day_categ × diet"))
  
  ggsave(
    filename = file.path(feat_dir, paste0(safe_feat, "_cond_eff_day_diet.png")),
    plot     = p_ce,
    width    = 7, height = 5, dpi = 1200
  )
  
  emm   <- emmeans(
    fit,
    ~ diet | day_categ,
    re_formula = NA
  )
  contr <- contrast(emm, method = "revpairwise")
  
  emm_df <- as.data.frame(emm) %>%
    mutate(feature = feature, family = selected_family_name)
  
  contr_df <- as.data.frame(contr) %>%
    mutate(feature = feature, family = selected_family_name)
  
  write.csv(
    emm_df,
    file = file.path(feat_dir, paste0(safe_feat, "_emmeans_diet_by_day.csv")),
    row.names = FALSE
  )
  
  write.csv(
    contr_df,
    file = file.path(feat_dir, paste0(safe_feat, "_contrast_diet_by_day.csv")),
    row.names = FALSE
  )
  
  contr_draws <- gather_emmeans_draws(contr)
  
  contr_summ <- contr_draws %>%
    mutate(feature = feature, family = selected_family_name) %>%
    group_by(feature, family, day_categ, contrast) %>%
    summarise(
      estimate = mean(.value),
      lower_95 = quantile(.value, 0.025),
      upper_95 = quantile(.value, 0.975),
      prob_gt0 = mean(.value > 0),
      .groups  = "drop"
    )
  
  write.csv(
    contr_summ,
    file = file.path(feat_dir, paste0(safe_feat, "_contrast_diet_by_day_posterior.csv")),
    row.names = FALSE
  )
  
  loo_res <- tryCatch(
    loo(fit, moment_match = TRUE),
    error = function(e) NULL
  )
  
  if (!is.null(loo_res)) {
    loo_df <- data.frame(
      feature     = feature,
      family      = selected_family_name,
      looic       = loo_res$estimates["looic", "Estimate"],
      looic_se    = loo_res$estimates["looic", "SE"],
      elpd_loo    = loo_res$estimates["elpd_loo", "Estimate"],
      elpd_loo_se = loo_res$estimates["elpd_loo", "SE"]
    )
    
    write.csv(
      loo_df,
      file = file.path(feat_dir, paste0(safe_feat, "_loo_summary.csv")),
      row.names = FALSE
    )
  } else {
    loo_df <- NULL
  }
  
  list(
    fit              = fit,
    emmeans_summary  = emm_df,
    contrast_summary = contr_summ,
    contrast_table   = contr_df,
    loo_summary      = loo_df,
    family           = selected_family_name
  )
}


# Loop over all features
all_emmeans   <- list()
all_contrasts <- list()
all_contr_raw <- list()
all_loo       <- list()
all_models    <- list()

for (feat in histone_features) {
  message("Running final best-family model for feature: ", feat)
  
  res <- fit_best_final_one(
    feature       = feat,
    family_lookup = family_lookup,
    base_dir      = base_dir
  )
  
  if (!is.null(res)) {
    all_models[[feat]]    <- data.frame(feature = feat, family = res$family)
    all_emmeans[[feat]]   <- res$emmeans_summary
    all_contrasts[[feat]] <- res$contrast_summary
    all_contr_raw[[feat]] <- res$contrast_table
    if (!is.null(res$loo_summary)) all_loo[[feat]] <- res$loo_summary
  }
}

emmeans_df      <- bind_rows(all_emmeans)
contrast_df     <- bind_rows(all_contrasts)
contrast_raw_df <- bind_rows(all_contr_raw)
loo_df_all      <- bind_rows(all_loo)
models_df       <- bind_rows(all_models)

write.csv(
  emmeans_df,
  file.path(summary_dir, "all_features_emmeans_diet_by_day_best_family.csv"),
  row.names = FALSE
)

write.csv(
  contrast_df,
  file.path(summary_dir, "all_features_contrast_diet_by_day_best_family_posterior.csv"),
  row.names = FALSE
)

write.csv(
  contrast_raw_df,
  file.path(summary_dir, "all_features_contrast_diet_by_day_best_family.csv"),
  row.names = FALSE
)

write.csv(
  loo_df_all,
  file.path(summary_dir, "all_features_final_loo_summary_best_family.csv"),
  row.names = FALSE
)

write.csv(
  models_df,
  file.path(summary_dir, "all_features_selected_family_used.csv"),
  row.names = FALSE
)

