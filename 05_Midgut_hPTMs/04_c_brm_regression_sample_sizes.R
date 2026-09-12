suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
})

# Paths

base_dir <- paste0(
  "04_brms_regression/",
  "brms_best_family_final_batch_fixed_20260628_162843"
)

summary_dir <- file.path(base_dir, "summary_tables")

# Load original model inputs

histone_data <- readRDS(
  "02_process_matrix/H3_H4_5PTM_midgut.rds"
) %>%
  as.data.frame() %>%
  rownames_to_column("sample_id")

metadata <- data.table::fread(
  "STPN2309_metadata_all.tsv"
) %>%
  as.data.frame()

names(metadata)[1] <- "sample_id"

# Retain only features present in the final contrast table

features_used <- read_csv(
  file.path(
    summary_dir,
    "all_features_contrast_diet_by_day_best_family_posterior.csv"
  ),
  show_col_types = FALSE
) %>%
  filter(contrast == "HC - NC") %>%
  distinct(feature) %>%
  pull(feature)

# Reproduce the feature-specific filtering used before each brms model

sample_sizes <- histone_data %>%
  pivot_longer(
    cols = -sample_id,
    names_to = "feature",
    values_to = "response"
  ) %>%
  filter(feature %in% features_used) %>%
  left_join(
    metadata %>%
      select(
        sample_id,
        diet,
        day_categ,
        gut_hptm_batch
      ),
    by = "sample_id"
  ) %>%
  filter(
    diet %in% c("HC", "NC"),
    !is.na(response),
    !is.na(day_categ),
    !is.na(gut_hptm_batch)
  ) %>%
  count(
    feature,
    day_categ,
    diet,
    name = "n"
  ) %>%
  arrange(
    feature,
    day_categ,
    diet
  )

# Long output

write_csv(
  sample_sizes,
  file.path(
    summary_dir,
    "midgut_hPTM_sample_sizes_long.csv"
  )
)

# HC and NC in separate columns

sample_sizes_wide <- sample_sizes %>%
  pivot_wider(
    names_from = diet,
    values_from = n,
    values_fill = 0
  ) %>%
  arrange(
    feature,
    day_categ
  )

write_csv(
  sample_sizes_wide,
  file.path(
    summary_dir,
    "midgut_hPTM_sample_sizes_HC_NC.csv"
  )
)

print(sample_sizes_wide)