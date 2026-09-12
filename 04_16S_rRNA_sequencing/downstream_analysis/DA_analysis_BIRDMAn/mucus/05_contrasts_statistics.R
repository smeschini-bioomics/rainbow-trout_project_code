suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
})

# Configuration
TISSUE <- "Mucus"   # change to "Mucus" for the mucus analysis

F_METADATA <- "output_files/metadata_model_used_for_BIRDMAn.tsv"
F_CONTRAST <- "export_for_R/HCminusNC_contrasts_posterior_summary_LONG.tsv"
F_DIAGNOSTICS <- file.path(
  "export_for_R",
  "diagnostics",
  "delta_contrast_diagnostics_LONG.tsv"
)

OUTPUT_DIR <- paste0(
  "DA_",
  tolower(TISSUE),
  "_statistics"
)

DAY_LEVELS <- c(
  "day_01",
  "day_02",
  "day_03",
  "day_04",
  "day_10",
  "day_15",
  "day_22"
)

HDI_PROB <- "0.95"
HDI_LOW_COL <- paste0("hdi_lower_", HDI_PROB)
HDI_HIGH_COL <- paste0("hdi_upper_", HDI_PROB)

PD_THRESHOLD <- 0.95
RHAT_THRESHOLD <- 1.05
ESS_TAIL_THRESHOLD <- 400

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# Helpers
pick_column <- function(df, candidates, label) {
  hit <- candidates[candidates %in% names(df)]

  if (length(hit) == 0) {
    stop(
      "Could not find the ",
      label,
      " column. Expected one of: ",
      paste(candidates, collapse = ", ")
    )
  }

  hit[1]
}

normalise_day <- function(x) {
  day_number <- readr::parse_number(
    as.character(x)
  )

  ifelse(
    is.na(day_number),
    NA_character_,
    sprintf("day_%02d", day_number)
  )
}

extract_genus <- function(feature) {
  genus <- str_match(
    feature,
    "g__([^;|\\s]+)"
  )[, 2]

  ifelse(
    is.na(genus),
    NA_character_,
    genus
  )
}

bad_genus_values <- c(
  "",
  "NA",
  "N/A",
  "none",
  "None",
  "unclassified",
  "Unclassified",
  "uncultured",
  "Unknown"
)

# Load exported files
metadata <- read_tsv(
  F_METADATA,
  show_col_types = FALSE
)

contrast_df <- read_tsv(
  F_CONTRAST,
  show_col_types = FALSE
)

diagnostics_df <- read_tsv(
  F_DIAGNOSTICS,
  show_col_types = FALSE
)

required_contrast_columns <- c(
  "day",
  "feature",
  "asv_id",
  "mean",
  "PD",
  HDI_LOW_COL,
  HDI_HIGH_COL
)

missing_contrast_columns <- setdiff(
  required_contrast_columns,
  names(contrast_df)
)

if (length(missing_contrast_columns) > 0) {
  stop(
    "Contrast file missing column(s): ",
    paste(
      missing_contrast_columns,
      collapse = ", "
    )
  )
}

required_diagnostic_columns <- c(
  "day",
  "asv_id",
  "r_hat",
  "ess_tail"
)

missing_diagnostic_columns <- setdiff(
  required_diagnostic_columns,
  names(diagnostics_df)
)

if (length(missing_diagnostic_columns) > 0) {
  stop(
    "Diagnostics file missing column(s): ",
    paste(
      missing_diagnostic_columns,
      collapse = ", "
    )
  )
}

# Identify metadata columns
sample_col <- pick_column(
  metadata,
  c("SampleID", "sample", "sample_id"),
  "sample ID"
)

sample_type_col <- pick_column(
  metadata,
  c("sample_type", "SampleType", "compartment"),
  "sample type"
)

diet_col <- pick_column(
  metadata,
  c("Diet", "diet"),
  "diet"
)

day_col <- pick_column(
  metadata,
  c("Day", "day", "day_categ"),
  "day"
)

# Sample numbers from the exact metadata used by BIRDMAn
metadata_used <- metadata %>%
  transmute(
    sample_id = as.character(
      .data[[sample_col]]
    ),
    sample_type = as.character(
      .data[[sample_type_col]]
    ),
    diet = toupper(
      trimws(
        as.character(
          .data[[diet_col]]
        )
      )
    ),
    day = normalise_day(
      .data[[day_col]]
    )
  ) %>%
  filter(
    tolower(sample_type) ==
      tolower(TISSUE),
    diet %in% c("HC", "NC"),
    day %in% DAY_LEVELS
  ) %>%
  distinct(
    sample_id,
    .keep_all = TRUE
  )

sample_counts <- metadata_used %>%
  count(
    day,
    diet,
    name = "n"
  ) %>%
  complete(
    day = DAY_LEVELS,
    diet = c("HC", "NC"),
    fill = list(n = 0)
  ) %>%
  pivot_wider(
    names_from = diet,
    values_from = n
  ) %>%
  mutate(
    tissue = TISSUE,
    n_total = HC + NC,
    day = factor(
      day,
      levels = DAY_LEVELS
    )
  ) %>%
  arrange(day) %>%
  mutate(
    day = as.character(day)
  )

# Prepare complete posterior contrast output
contrast_all <- contrast_df %>%
  transmute(
    day = normalise_day(day),
    feature = str_trim(
      as.character(feature)
    ),
    asv_id = str_trim(
      as.character(asv_id)
    ),
    posterior_mean = as.numeric(mean),
    PD = as.numeric(PD),
    hdi_lower = as.numeric(
      .data[[HDI_LOW_COL]]
    ),
    hdi_upper = as.numeric(
      .data[[HDI_HIGH_COL]]
    )
  ) %>%
  filter(
    day %in% DAY_LEVELS
  ) %>%
  mutate(
    hdi_excludes_zero =
      hdi_lower > 0 |
      hdi_upper < 0,
    significant =
      hdi_excludes_zero &
      PD > PD_THRESHOLD,
    direction = case_when(
      posterior_mean > 0 ~ "Higher in HC",
      posterior_mean < 0 ~ "Lower in HC",
      TRUE ~ "No difference"
    ),
    genus = extract_genus(feature),
    genus_annotated =
      !is.na(genus) &
      !(genus %in% bad_genus_values)
  )

# ASV-level diagnostics, matching the heatmap script
diagnostic_summary <- diagnostics_df %>%
  transmute(
    day = normalise_day(day),
    asv_id = str_trim(
      as.character(asv_id)
    ),
    r_hat = as.numeric(r_hat),
    ess_tail = as.numeric(ess_tail)
  ) %>%
  filter(
    day %in% DAY_LEVELS
  ) %>%
  group_by(asv_id) %>%
  summarise(
    worst_rhat = max(
      r_hat,
      na.rm = TRUE
    ),
    min_ess_tail = min(
      ess_tail,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    diagnostic_pass =
      worst_rhat < RHAT_THRESHOLD &
      min_ess_tail > ESS_TAIL_THRESHOLD
  )

contrast_annotated <- contrast_all %>%
  left_join(
    diagnostic_summary,
    by = "asv_id"
  ) %>%
  mutate(
    diagnostic_pass =
      replace_na(
        diagnostic_pass,
        FALSE
      ),
    display_eligible =
      diagnostic_pass &
      genus_annotated
  )

# Counts by day
asv_counts_by_day <- contrast_annotated %>%
  group_by(day) %>%
  summarise(
    n_ASVs_tested =
      n_distinct(asv_id),
    n_ASVs_diagnostic_pass =
      n_distinct(
        asv_id[diagnostic_pass]
      ),
    n_ASVs_genus_annotated =
      n_distinct(
        asv_id[genus_annotated]
      ),
    n_ASVs_display_eligible =
      n_distinct(
        asv_id[display_eligible]
      ),
    n_significant_ASVs_all =
      n_distinct(
        asv_id[significant]
      ),
    n_higher_HC_all =
      n_distinct(
        asv_id[
          significant &
          posterior_mean > 0
        ]
      ),
    n_lower_HC_all =
      n_distinct(
        asv_id[
          significant &
          posterior_mean < 0
        ]
      ),
    n_significant_ASVs_display_eligible =
      n_distinct(
        asv_id[
          significant &
          display_eligible
        ]
      ),
    n_higher_HC_display_eligible =
      n_distinct(
        asv_id[
          significant &
          display_eligible &
          posterior_mean > 0
        ]
      ),
    n_lower_HC_display_eligible =
      n_distinct(
        asv_id[
          significant &
          display_eligible &
          posterior_mean < 0
        ]
      ),
    .groups = "drop"
  )

summary_by_day <- sample_counts %>%
  left_join(
    asv_counts_by_day,
    by = "day"
  ) %>%
  arrange(
    match(
      day,
      DAY_LEVELS
    )
  )

# Unique ASVs across all days
overall_summary <- data.frame(
  tissue = TISSUE,
  n_samples_total =
    n_distinct(
      metadata_used$sample_id
    ),
  n_HC_total =
    sum(
      metadata_used$diet == "HC"
    ),
  n_NC_total =
    sum(
      metadata_used$diet == "NC"
    ),
  n_ASVs_tested =
    n_distinct(
      contrast_annotated$asv_id
    ),
  n_ASVs_diagnostic_pass =
    n_distinct(
      contrast_annotated$asv_id[
        contrast_annotated$diagnostic_pass
      ]
    ),
  n_ASVs_display_eligible =
    n_distinct(
      contrast_annotated$asv_id[
        contrast_annotated$display_eligible
      ]
    ),
  n_unique_significant_ASVs_all =
    n_distinct(
      contrast_annotated$asv_id[
        contrast_annotated$significant
      ]
    ),
  n_unique_significant_ASVs_displayed =
    n_distinct(
      contrast_annotated$asv_id[
        contrast_annotated$significant &
        contrast_annotated$display_eligible
      ]
    )
)

# Export
write_csv(
  summary_by_day,
  file.path(
    OUTPUT_DIR,
    "microbiome_DA_sample_size_and_ASV_counts_by_day.csv"
  )
)

write_csv(
  overall_summary,
  file.path(
    OUTPUT_DIR,
    "microbiome_DA_overall_summary.csv"
  )
)

write_tsv(
  contrast_annotated %>%
    filter(significant) %>%
    arrange(
      factor(
        day,
        levels = DAY_LEVELS
      ),
      desc(abs(posterior_mean))
    ),
  file.path(
    OUTPUT_DIR,
    "microbiome_DA_all_significant_ASVs_by_day.tsv"
  )
)

write_tsv(
  contrast_annotated %>%
    filter(
      significant,
      display_eligible
    ) %>%
    arrange(
      factor(
        day,
        levels = DAY_LEVELS
      ),
      desc(abs(posterior_mean))
    ),
  file.path(
    OUTPUT_DIR,
    "microbiome_DA_significant_ASVs_display_eligible.tsv"
  )
)

print(summary_by_day)
print(overall_summary)

message(
  "Saved outputs in: ",
  OUTPUT_DIR
)

