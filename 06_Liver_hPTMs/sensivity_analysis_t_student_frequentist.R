suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(readr)
  library(data.table)
  library(glmmTMB)
  library(emmeans)
  library(broom.mixed)
  library(stringr)
  library(ggplot2)
  library(ggtree)
  library(cowplot)
})

# ------------------------------------------------------------
# Input data
# ------------------------------------------------------------

histone_data <- readRDS(
  "02_process_matrix/H3_H4_8PTM_liver.rds"
) %>%
  as.data.frame()

histone_metadata <- data.table::fread(
  "STPN2309_metadata_all.tsv"
) %>%
  as.data.frame()

rownames(histone_metadata) <- histone_metadata[[1]]
histone_metadata <- histone_metadata[, -1, drop = FALSE]

# ------------------------------------------------------------
# Settings and output directories
# ------------------------------------------------------------

BATCH_FIXED_COL <- "trout_liver_batch"

DAY_LEVELS <- c(
  "day_01", "day_02", "day_03", "day_04",
  "day_10", "day_15", "day_22"
)

output_root <- "05_glmmTMB_student_t_regression"

base_dir <- file.path(
  output_root,
  "trout_liver_batch"
)

fits_dir <- file.path(
  base_dir,
  "fits"
)

summary_dir <- file.path(
  base_dir,
  "summary_tables"
)

dir.create(fits_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Align samples
# ------------------------------------------------------------

common_ids <- sort(
  intersect(
    rownames(histone_data),
    rownames(histone_metadata)
  )
)

if (length(common_ids) == 0) {
  stop("No shared sample IDs between the hPTM matrix and metadata.")
}

histone_data <- histone_data[
  common_ids,
  ,
  drop = FALSE
]

histone_metadata <- histone_metadata[
  common_ids,
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(histone_data),
    rownames(histone_metadata)
  )
)

histone_metadata <- histone_metadata %>%
  filter(
    diet %in% c("NC", "HC")
  )

if (nrow(histone_metadata) == 0) {
  stop("No NC/HC samples remain after filtering.")
}

histone_data <- histone_data[
  rownames(histone_metadata),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(histone_data),
    rownames(histone_metadata)
  )
)

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
    diet = factor(
      diet,
      levels = c("NC", "HC")
    ),
    day_categ = factor(
      day_categ,
      levels = DAY_LEVELS
    ),
    hptm_batch = factor(
      .data[[BATCH_FIXED_COL]]
    )
  )

histone_metadata$diet <- droplevels(
  histone_metadata$diet
)

histone_metadata$day_categ <- droplevels(
  histone_metadata$day_categ
)

histone_metadata$hptm_batch <- droplevels(
  histone_metadata$hptm_batch
)

if (any(is.na(histone_metadata$day_categ))) {
  stop("Unexpected values detected in day_categ.")
}

if (nlevels(histone_metadata$diet) != 2) {
  stop("Both NC and HC must be present.")
}

if (nlevels(histone_metadata$hptm_batch) < 2) {
  stop("The liver hPTM batch variable has fewer than two levels.")
}

cat("\nSamples retained:", nrow(histone_metadata), "\n")

cat("\nDiet counts:\n")
print(table(histone_metadata$diet))

cat("\nDiet x day table:\n")
print(
  table(
    histone_metadata$diet,
    histone_metadata$day_categ
  )
)

cat("\nBatch counts:\n")
print(table(histone_metadata$hptm_batch))

histone_features <- colnames(histone_data)

# ------------------------------------------------------------
# Fit one Student-t model
# ------------------------------------------------------------

fit_student_one <- function(
    feature,
    feature_df,
    meta_df,
    save_dir
) {

  message(
    "Fitting Student-t model: ",
    feature
  )

  stopifnot(
    identical(
      rownames(feature_df),
      rownames(meta_df)
    )
  )

  dat_feat <- meta_df

  dat_feat$response <- as.numeric(
    feature_df[
      rownames(meta_df),
      feature
    ]
  )

  dat_feat <- dat_feat %>%
    filter(
      !is.na(response)
    ) %>%
    mutate(
      diet = factor(
        diet,
        levels = c("NC", "HC")
      ),
      day_categ = factor(
        day_categ,
        levels = DAY_LEVELS
      ),
      hptm_batch = factor(
        hptm_batch
      )
    )

  dat_feat$diet <- droplevels(
    dat_feat$diet
  )

  dat_feat$day_categ <- droplevels(
    dat_feat$day_categ
  )

  dat_feat$hptm_batch <- droplevels(
    dat_feat$hptm_batch
  )

  n_after <- nrow(dat_feat)

  if (n_after < 6) {
    warning(
      "Too few non-missing values for ",
      feature,
      " - skipping."
    )
    return(NULL)
  }

  if (nlevels(dat_feat$diet) != 2) {
    warning(
      "Both diets are not represented for ",
      feature,
      " - skipping."
    )
    return(NULL)
  }

  if (nlevels(dat_feat$hptm_batch) < 2) {
    warning(
      "Fewer than two batch levels for ",
      feature,
      " - skipping."
    )
    return(NULL)
  }

  fit <- glmmTMB(
    response ~ day_categ * diet + hptm_batch,
    data = dat_feat,
    family = t_family(
      link = "identity"
    ),
    REML = FALSE,
    control = glmmTMBControl(
      optCtrl = list(
        iter.max = 10000,
        eval.max = 10000
      )
    )
  )

  safe_feature <- gsub(
    "[^A-Za-z0-9]+",
    "_",
    feature
  )

  feature_dir <- file.path(
    save_dir,
    safe_feature
  )

  dir.create(
    feature_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  saveRDS(
    fit,
    file.path(
      feature_dir,
      paste0(
        safe_feature,
        "_student_t_glmmTMB.rds"
      )
    )
  )

  writeLines(
    capture.output(
      summary(fit)
    ),
    file.path(
      feature_dir,
      paste0(
        safe_feature,
        "_model_summary.txt"
      )
    )
  )

  convergence_code <- fit$fit$convergence
  positive_hessian <- isTRUE(
    fit$sdr$pdHess
  )

  usable_for_inference <- (
    convergence_code == 0 &&
    positive_hessian
  )

  diagnostics <- tibble(
    feature = feature,
    family = "Student-t",
    formula_used =
      "response ~ day_categ * diet + hptm_batch",
    n_obs = n_after,
    convergence_code = convergence_code,
    convergence_message = fit$fit$message,
    positive_definite_hessian = positive_hessian,
    usable_for_inference = usable_for_inference,
    logLik = suppressWarnings(
      as.numeric(logLik(fit))
    ),
    AIC = suppressWarnings(
      AIC(fit)
    ),
    BIC = suppressWarnings(
      BIC(fit)
    ),
    residual_scale = suppressWarnings(
      sigma(fit)
    )
  )

  if (!usable_for_inference) {
    warning(
      "Model excluded from inference for ",
      feature,
      ": convergence code = ",
      convergence_code,
      ", positive-definite Hessian = ",
      positive_hessian,
      "."
    )

    return(
      list(
        diagnostics = diagnostics,
        emmeans = NULL,
        contrasts = NULL
      )
    )
  }

  fixed_effects <- broom.mixed::tidy(
    fit,
    effects = "fixed",
    component = "cond",
    conf.int = TRUE,
    conf.level = 0.95
  ) %>%
    mutate(
      feature = feature
    ) %>%
    relocate(
      feature
    )

  readr::write_csv(
    fixed_effects,
    file.path(
      feature_dir,
      paste0(
        safe_feature,
        "_fixed_effects.csv"
      )
    )
  )

  emm <- emmeans(
    fit,
    ~ diet | day_categ
  )

  emm_df <- as.data.frame(
    summary(
      emm,
      infer = c(TRUE, FALSE),
      level = 0.95
    )
  ) %>%
    mutate(
      feature = feature
    ) %>%
    relocate(
      feature
    )

  contrast_object <- contrast(
    emm,
    method = list(
      "HC - NC" = c(-1, 1)
    ),
    adjust = "none"
  )

  contrast_df <- as.data.frame(
    summary(
      contrast_object,
      infer = c(TRUE, TRUE),
      level = 0.95,
      adjust = "none"
    )
  )

  lower_column <- intersect(
    c(
      "lower.CL",
      "asymp.LCL"
    ),
    names(contrast_df)
  )

  upper_column <- intersect(
    c(
      "upper.CL",
      "asymp.UCL"
    ),
    names(contrast_df)
  )

  statistic_column <- intersect(
    c(
      "t.ratio",
      "z.ratio"
    ),
    names(contrast_df)
  )

  if (
    length(lower_column) != 1 ||
    length(upper_column) != 1
  ) {
    stop(
      "Could not identify confidence interval columns for ",
      feature,
      "."
    )
  }

  contrast_df$lower_95 <- contrast_df[[lower_column]]
  contrast_df$upper_95 <- contrast_df[[upper_column]]

  if (length(statistic_column) == 1) {
    contrast_df$statistic <- contrast_df[[statistic_column]]
  } else {
    contrast_df$statistic <- NA_real_
  }

  if (!"df" %in% colnames(contrast_df)) {
    contrast_df$df <- Inf
  }

  contrast_df <- contrast_df %>%
    transmute(
      feature = feature,
      day_categ = as.character(day_categ),
      contrast = as.character(contrast),
      estimate = estimate,
      SE = SE,
      df = df,
      lower_95 = lower_95,
      upper_95 = upper_95,
      statistic = statistic,
      p_value = p.value
    )

  readr::write_csv(
    emm_df,
    file.path(
      feature_dir,
      paste0(
        safe_feature,
        "_emmeans_diet_by_day.csv"
      )
    )
  )

  readr::write_csv(
    contrast_df,
    file.path(
      feature_dir,
      paste0(
        safe_feature,
        "_HC_minus_NC_by_day.csv"
      )
    )
  )

  list(
    diagnostics = diagnostics,
    emmeans = emm_df,
    contrasts = contrast_df
  )
}

# ------------------------------------------------------------
# Run all liver hPTMs
# ------------------------------------------------------------

results <- vector(
  mode = "list",
  length = length(histone_features)
)

names(results) <- histone_features

error_log <- list()

for (feature in histone_features) {

  results[[feature]] <- tryCatch(

    fit_student_one(
      feature = feature,
      feature_df = histone_data,
      meta_df = histone_metadata,
      save_dir = fits_dir
    ),

    error = function(e) {

      error_log[[feature]] <<- conditionMessage(e)

      warning(
        "Model failed for ",
        feature,
        ": ",
        conditionMessage(e)
      )

      NULL
    }
  )
}

successful_results <- purrr::compact(
  results
)

saveRDS(
  error_log,
  file.path(
    summary_dir,
    "student_t_model_error_log.rds"
  )
)

if (length(successful_results) == 0) {
  stop(
    "No Student-t models were returned. ",
    "Inspect student_t_model_error_log.rds."
  )
}

diagnostics_df <- bind_rows(
  map(
    successful_results,
    "diagnostics"
  )
)

usable_results <- successful_results[
  vapply(
    successful_results,
    function(x) {
      isTRUE(
        x$diagnostics$usable_for_inference[1]
      )
    },
    logical(1)
  )
]

if (length(usable_results) == 0) {
  readr::write_csv(
    diagnostics_df,
    file.path(
      summary_dir,
      "all_features_student_t_model_diagnostics.csv"
    )
  )

  stop(
    "No converged Student-t models with a positive-definite Hessian ",
    "were available for inference."
  )
}

emmeans_df <- bind_rows(
  map(
    usable_results,
    "emmeans"
  )
)

contrast_df <- bind_rows(
  map(
    usable_results,
    "contrasts"
  )
) %>%
  group_by(feature) %>%
  mutate(
    FDR = p.adjust(
      p_value,
      method = "BH"
    )
  ) %>%
  ungroup() %>%
  mutate(
    direction = if_else(
      estimate >= 0,
      "HC > NC",
      "HC < NC"
    ),
    abs_est = abs(
      estimate
    )
  )

readr::write_csv(
  diagnostics_df,
  file.path(
    summary_dir,
    "all_features_student_t_model_diagnostics.csv"
  )
)

readr::write_csv(
  emmeans_df,
  file.path(
    summary_dir,
    "all_features_student_t_emmeans_diet_by_day.csv"
  )
)

readr::write_csv(
  contrast_df,
  file.path(
    summary_dir,
    "all_features_student_t_HC_minus_NC_by_day.csv"
  )
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    base_dir,
    "session_info.txt"
  )
)

cat(
  "\nStudent-t frequentist analysis completed.\n",
  "Models returned: ",
  length(successful_results),
  " / ",
  length(histone_features),
  "\nModels retained for inference: ",
  length(usable_results),
  " / ",
  length(histone_features),
  "\nResults saved in: ",
  normalizePath(base_dir),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# Frequentist dot heatmap
# ------------------------------------------------------------

contrast_file <- file.path(
  summary_dir,
  "all_features_student_t_HC_minus_NC_by_day.csv"
)

if (!file.exists(contrast_file)) {
  stop(
    "Contrast file not found: ",
    contrast_file
  )
}

contrast_df <- readr::read_csv(
  contrast_file,
  show_col_types = FALSE
)

normalize_ptm <- function(x) {
  vapply(
    x,
    function(entry) {
      parts <- strsplit(
        entry,
        "\\|"
      )[[1]]

      normalized_parts <- vapply(
        parts,
        function(part) {
          matched <- regmatches(
            part,
            regexec(
              "^(H[34])(K|R)([0-9]{1,3})([A-Za-z0-9]+)$",
              part
            )
          )[[1]]

          if (
            length(matched) < 5 ||
            is.na(matched[1])
          ) {
            return(NA_character_)
          }

          paste0(
            matched[2],
            matched[3],
            matched[4],
            tolower(matched[5])
          )
        },
        character(1)
      )

      if (anyNA(normalized_parts)) {
        return(NA_character_)
      }

      paste(
        normalized_parts,
        collapse = "|"
      )
    },
    character(1),
    USE.NAMES = FALSE
  )
}

contrast_df <- contrast_df %>%
  mutate(
    feature = normalize_ptm(
      as.character(feature)
    ),
    day_categ = as.character(
      day_categ
    ),
    direction = factor(
      direction,
      levels = c(
        "HC < NC",
        "HC > NC"
      )
    ),
    abs_est = abs(
      estimate
    ),
    fdr_label = case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01 ~ "**",
      FDR < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  filter(
    !is.na(feature),
    contrast == "HC - NC"
  )

if (nrow(contrast_df) == 0) {
  stop("No valid HC - NC contrasts were found.")
}

# Cross-layer activation summary
TISSUE <- "Liver"
LAYER <- "hPTM"

set.seed(123)

boot_median_ci <- function(
    x,
    n_boot = 2000,
    conf = 0.95
) {
  x <- x[!is.na(x)]

  if (length(x) < 2) {
    return(
      c(
        median = if (length(x) == 1) x else NA_real_,
        lower = NA_real_,
        upper = NA_real_
      )
    )
  }

  bootstrap_medians <- replicate(
    n_boot,
    median(
      sample(
        x,
        length(x),
        replace = TRUE
      )
    )
  )

  alpha <- (1 - conf) / 2

  c(
    median = median(x),
    lower = unname(
      quantile(
        bootstrap_medians,
        alpha,
        na.rm = TRUE
      )
    ),
    upper = unname(
      quantile(
        bootstrap_medians,
        1 - alpha,
        na.rm = TRUE
      )
    )
  )
}

activation_summary <- contrast_df %>%
  mutate(
    day_num = as.numeric(
      gsub(
        "day_",
        "",
        day_categ
      )
    )
  ) %>%
  group_by(day_num) %>%
  summarise(
    boot = list(
      boot_median_ci(
        abs(estimate)
      )
    ),
    n_tested = n(),
    .groups = "drop"
  ) %>%
  mutate(
    median_abs_delta = map_dbl(
      boot,
      "median"
    ),
    ci_lower = map_dbl(
      boot,
      "lower"
    ),
    ci_upper = map_dbl(
      boot,
      "upper"
    )
  ) %>%
  select(-boot) %>%
  rename(day = day_num) %>%
  right_join(
    tibble(
      day = c(
        1,
        2,
        3,
        4,
        10,
        15,
        22
      )
    ),
    by = "day"
  ) %>%
  mutate(
    tissue = TISSUE,
    layer = LAYER
  ) %>%
  arrange(day) %>%
  mutate(
    row_id = paste(
      tissue,
      layer,
      day,
      sep = "_"
    )
  ) %>%
  column_to_rownames(
    "row_id"
  ) %>%
  select(
    tissue,
    layer,
    day,
    median_abs_delta,
    ci_lower,
    ci_upper,
    n_tested
  )

saveRDS(
  activation_summary,
  file.path(
    summary_dir,
    "hptm_liver_student_t_median_abs_delta.rds"
  )
)

readr::write_csv(
  activation_summary %>%
    rownames_to_column(
      "row_id"
    ),
  file.path(
    summary_dir,
    "hptm_liver_student_t_median_abs_delta.csv"
  )
)

combined <- activation_summary %>%
  rownames_to_column(
    "row_id"
  ) %>%
  mutate(
    day_lab = factor(
      paste0(
        "Day ",
        day
      ),
      levels = paste0(
        "Day ",
        c(
          1,
          2,
          3,
          4,
          10,
          15,
          22
        )
      )
    ),
    day_idx = as.numeric(
      day_lab
    ),
    tissue_layer = factor(
      paste(
        tissue,
        layer,
        sep = " - "
      )
    )
  )

p_activation <- ggplot(
  combined,
  aes(
    x = day_idx,
    y = median_abs_delta,
    group = tissue_layer
  )
) +
  geom_ribbon(
    aes(
      ymin = ci_lower,
      ymax = ci_upper
    ),
    fill = "#44a7c4",
    alpha = 0.4
  ) +
  geom_line(
    linewidth = 1,
    color = "#2171b5"
  ) +
  geom_point(
    aes(
      size = n_tested
    ),
    color = "#2171b5"
  ) +
  scale_x_continuous(
    breaks = 1:7,
    labels = levels(
      combined$day_lab
    )
  ) +
  scale_size_continuous(
    range = c(
      2,
      8
    ),
    name = "n features tested"
  ) +
  labs(
    x = NULL,
    y = "Median |HC - NC estimate|\n(bootstrap 95% CI)",
    title =
      "Frequentist Student-t: liver hPTM effect dynamics"
  ) +
  theme_grey(
    base_size = 13
  ) +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5
    )
  )

ggsave(
  filename = file.path(
    summary_dir,
    "student_t_effect_size_dynamics.pdf"
  ),
  plot = p_activation,
  width = 8,
  height = 5
)

ggsave(
  filename = file.path(
    summary_dir,
    "student_t_effect_size_dynamics.tiff"
  ),
  plot = p_activation,
  width = 8,
  height = 5,
  dpi = 1200,
  compression = "lzw"
)

# Cluster PTMs using HC - NC estimates across days
heat_df <- contrast_df %>%
  mutate(
    day_num = as.numeric(
      gsub(
        "day_",
        "",
        day_categ
      )
    ),
    day_label = paste0(
      "Day ",
      day_num
    )
  )

clust_df <- heat_df %>%
  select(
    feature,
    day_num,
    estimate
  ) %>%
  arrange(
    day_num
  ) %>%
  pivot_wider(
    names_from = day_num,
    values_from = estimate
  )

complete_rows <- complete.cases(
  clust_df[
    ,
    -1,
    drop = FALSE
  ]
)

clust_df <- clust_df[
  complete_rows,
  ,
  drop = FALSE
]

if (nrow(clust_df) < 2) {
  stop(
    "Fewer than two complete PTM profiles are available for clustering."
  )
}

row_names <- clust_df$feature

clust_mat <- as.matrix(
  clust_df[
    ,
    -1,
    drop = FALSE
  ]
)

rownames(clust_mat) <- row_names

hc <- hclust(
  dist(
    clust_mat
  ),
  method = "complete"
)

ptm_order <- rownames(
  clust_mat
)[hc$order]

dendrogram_object <- as.dendrogram(
  hc
)

ggtree_plot <- ggtree::ggtree(
  dendrogram_object,
  branch.length = "none"
)

heat_df2 <- heat_df %>%
  filter(
    feature %in% ptm_order
  ) %>%
  mutate(
    feature = factor(
      feature,
      levels = ptm_order
    ),
    day_label = factor(
      day_label,
      levels = paste0(
        "Day ",
        sort(
          unique(
            day_num
          )
        )
      )
    )
  )

direction_colors <- c(
  "HC < NC" = "#44a7c4",
  "HC > NC" = "#f99943"
)

p_dot <- ggplot(
  heat_df2,
  aes(
    x = day_label,
    y = feature
  )
) +
  geom_point(
    aes(
      size = abs_est,
      fill = direction
    ),
    shape = 21,
    colour = "black",
    alpha = 0.9,
    stroke = 1
  ) +
  scale_fill_manual(
    name = "Direction",
    values = direction_colors,
    labels = c(
      "HC < NC",
      "HC > NC"
    ),
    guide = guide_legend(
      override.aes = list(
        size = 10
      ),
      order = 1
    )
  ) +
  scale_size_continuous(
    name = "|HC - NC| (model scale)",
    range = c(
      4,
      20
    ),
    guide = guide_legend(
      order = 2
    )
  ) +
  geom_text(
    aes(
      label = fdr_label
    ),
    size = 6,
    fontface = "bold",
    colour = "black",
    vjust = 0.5
  ) +
  geom_point(
    aes(
      shape = "FDR symbols"
    ),
    x = NA,
    y = NA,
    size = 0
  ) +
  scale_shape_manual(
    name = "FDR-adjusted significance",
    values = c(
      "FDR symbols" = NA
    ),
    labels = c(
      "FDR symbols" =
        "* FDR < 0.05\n** FDR < 0.01\n*** FDR < 0.001\nwithin each PTM across 7 days"
    ),
    guide = guide_legend(
      order = 3
    )
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw(
    base_size = 15
  ) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 18
    ),
    axis.text.y = element_text(
      size = 18
    ),
    legend.title = element_text(
      size = 14,
      face = "bold"
    ),
    legend.text = element_text(
      size = 14
    ),
    panel.grid.major = element_line(
      color = "white"
    ),
    panel.grid.minor = element_line(
      color = "white"
    ),
    panel.background = element_rect(
      fill = "grey90"
    )
  )

plot_grid <- cowplot::plot_grid(
  ggtree_plot,
  p_dot,
  nrow = 1,
  rel_widths = c(
    0.25,
    2.25
  ),
  align = "h"
)

print(
  plot_grid
)

ggsave(
  filename = file.path(
    summary_dir,
    "dotplot_student_t_frequentist_FDR_stars.tiff"
  ),
  plot = plot_grid,
  device = "tiff",
  width = 10,
  height = 12,
  dpi = 1200,
  compression = "lzw"
)

ggsave(
  filename = file.path(
    summary_dir,
    "dotplot_student_t_frequentist_FDR_stars.pdf"
  ),
  plot = plot_grid,
  device = "pdf",
  width = 10,
  height = 12
)

