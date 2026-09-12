library(dplyr)
library(stringr)
library(tidyr)

# Settings
load("01_MSDAP_DEA/2026-04-29_16-53-21/dataset.RData")

tissue_name <- "Liver"
out_dir <- "03_DEA_statistic"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Sample numbers after exclusion
sample_summary <- dataset$samples %>%
  filter(exclude == FALSE, diet %in% c("HC", "NC")) %>%
  mutate(
    diet_clean = diet,
    day = as.integer(str_extract(as.character(day_categorical), "\\d+"))
  ) %>%
  count(day, diet_clean, name = "n_samples") %>%
  pivot_wider(
    names_from = diet_clean,
    values_from = n_samples,
    values_fill = 0
  ) %>%
  mutate(
    n_total = HC + NC,
    tissue = tissue_name
  ) %>%
  select(tissue, day, HC, NC, n_total) %>%
  arrange(day)

sample_total <- sample_summary %>%
  summarise(
    tissue = tissue_name,
    total_HC = sum(HC),
    total_NC = sum(NC),
    total_samples = sum(n_total)
  )
sample_exclusion_summary <- dataset$samples %>%
  mutate(
    diet_clean = diet,
    day = as.integer(str_extract(as.character(day_categorical), "\\d+"))
  ) %>%
  count(exclude, diet_clean, day, name = "n_samples") %>%
  mutate(tissue = tissue_name) %>%
  select(tissue, exclude, diet = diet_clean, day, n_samples) %>%
  arrange(exclude, day, diet)

# Protein numbers per contrast and DEA algorithm
dea_contrast_summary <- dataset$de_proteins %>%
  mutate(
    tissue = tissue_name,
    day = as.integer(str_match(contrast, "day_(\\d+)")[, 2]),
    is_contaminant = str_detect(protein_id, fixed("Cont_")),
    sig = is.finite(qvalue) & qvalue < 0.05
  ) %>%
  filter(!is.na(day)) %>%
  group_by(tissue, day, contrast, dea_algorithm) %>%
  summarise(
    n_proteins_raw = n_distinct(protein_id),
    n_contaminants = n_distinct(protein_id[is_contaminant]),
    n_proteins_tested_no_contaminants = n_distinct(protein_id[!is_contaminant]),
    n_valid_qvalue = n_distinct(protein_id[!is_contaminant & is.finite(qvalue)]),
    n_sig_total = n_distinct(protein_id[!is_contaminant & sig]),
    n_sig_HC_gt_NC = n_distinct(
      protein_id[!is_contaminant & sig & foldchange.log2 > 0]
    ),
    n_sig_HC_lt_NC = n_distinct(
      protein_id[!is_contaminant & sig & foldchange.log2 < 0]
    ),
    pct_sig_of_tested = round(
      100 * n_sig_total / n_proteins_tested_no_contaminants,
      2
    ),
    .groups = "drop"
  ) %>%
  arrange(day, dea_algorithm)

#  tables
reviewer_contrast_table <- dea_contrast_summary %>%
  left_join(sample_summary, by = c("tissue", "day")) %>%
  select(
    tissue, day, contrast, dea_algorithm,
    HC, NC, n_total,
    n_proteins_raw, n_contaminants,
    n_proteins_tested_no_contaminants,
    n_valid_qvalue,
    n_sig_total, n_sig_HC_gt_NC, n_sig_HC_lt_NC,
    pct_sig_of_tested
  )

dea_day_summary <- reviewer_contrast_table %>%
  group_by(tissue, day, HC, NC, n_total) %>%
  summarise(
    n_algorithms = n_distinct(dea_algorithm),
    min_proteins_tested = min(n_proteins_tested_no_contaminants),
    median_proteins_tested = median(n_proteins_tested_no_contaminants),
    max_proteins_tested = max(n_proteins_tested_no_contaminants),
    min_sig = min(n_sig_total),
    median_sig = median(n_sig_total),
    max_sig = max(n_sig_total),
    min_pct_sig = min(pct_sig_of_tested),
    median_pct_sig = median(pct_sig_of_tested),
    max_pct_sig = max(pct_sig_of_tested),
    .groups = "drop"
  ) %>%
  arrange(day)

# 4. Global protein and sample summary
protein_global_summary <- dataset$de_proteins %>%
  summarise(
    tissue = tissue_name,
    total_unique_proteins_raw = n_distinct(protein_id),
    contaminant_proteins = n_distinct(
      protein_id[str_detect(protein_id, fixed("Cont_"))]
    ),
    total_unique_proteins_no_contaminants = n_distinct(
      protein_id[!str_detect(protein_id, fixed("Cont_"))]
    )
  )

global_summary <- protein_global_summary %>%
  bind_cols(sample_total %>% select(total_HC, total_NC, total_samples))

write.csv(
  sample_summary,
  file.path(out_dir, "Liver_sample_n_by_day_diet.csv"),
  row.names = FALSE
)

write.csv(
  sample_total,
  file.path(out_dir, "Liver_sample_n_total.csv"),
  row.names = FALSE
)

write.csv(
  sample_exclusion_summary,
  file.path(out_dir, "Liver_sample_exclusion_summary.csv"),
  row.names = FALSE
)

write.csv(
  reviewer_contrast_table,
  file.path(out_dir, "Liver_reviewer_contrast_table.csv"),
  row.names = FALSE
)

write.csv(
  dea_day_summary,
  file.path(out_dir, "Liver_DEA_day_summary.csv"),
  row.names = FALSE
)

write.csv(
  global_summary,
  file.path(out_dir, "Liver_global_summary.csv"),
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(out_dir, "sessionInfo.txt")
)

print(sample_summary, n = Inf)
print(sample_total)
print(global_summary)
print(dea_day_summary, n = Inf)

