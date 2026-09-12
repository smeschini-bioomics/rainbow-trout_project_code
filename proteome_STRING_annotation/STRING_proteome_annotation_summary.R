suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})


# 1) Paths and settings
string_terms_file <- "110079946.protein.enrichment.terms.v12.0.txt"
reference_fasta <- 
  "24_10_02_uniprot-ref_proteome_trout_UP000193380.fasta"
out_dir <- "annotation_statistics_reference_proteome"
TIFF_DPI <- 1200

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(string_terms_file))
stopifnot(file.exists(reference_fasta))



# 2) Helper functions
clean_character <- function(x) {
  x %>%
    as.character() %>%
    str_trim() %>%
    na_if("")
}


# Standard UniProt headers:
# >sp|A0A123ABC4|PROTEIN_NAME ...
# >tr|A0A123ABC4|PROTEIN_NAME ...
parse_uniprot_id <- function(headers) {
  
  headers_clean <- str_remove(headers, "^>")
  
  uniprot_id <- str_match(
    headers_clean,
    "^(?:sp|tr)\\|([^|]+)\\|"
  )[, 2]
  
  fallback_id <- str_extract(headers_clean, "^[^\\s|]+")
  
  dplyr::coalesce(uniprot_id, fallback_id)
}


save_plot <- function(plot_object, filename_stub, width, height) {
  
  ggsave(
    filename = file.path(out_dir, paste0(filename_stub, ".pdf")),
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    device = "pdf",
    bg = "white"
  )
  
  ggsave(
    filename = file.path(out_dir, paste0(filename_stub, ".tiff")),
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = TIFF_DPI,
    compression = "lzw",
    device = "tiff",
    bg = "white"
  )
}


theme_annotation <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, color = "grey35"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "grey25"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "none",
    plot.margin = margin(8, 32, 8, 8)
  )



# 3) Import STRING annotation file


ann_raw <- read_tsv(
  string_terms_file,
  comment = "#",
  col_names = c(
    "string_protein_id",
    "category",
    "term",
    "description"
  ),
  col_types = cols(.default = col_character()),
  progress = FALSE
)

ann <- ann_raw %>%
  mutate(
    protein_id = string_protein_id %>%
      str_remove("^[^.]+\\.") %>%
      clean_character(),
    
    category = clean_character(category),
    term = clean_character(term),
    description = clean_character(description)
  ) %>%
  filter(
    !is.na(protein_id),
    !is.na(category),
    !is.na(term),
    !str_detect(
      protein_id,
      regex("^(Cont_|REV_|DECOY_)", ignore_case = TRUE)
    )
  ) %>%
  distinct(
    protein_id,
    category,
    term,
    .keep_all = TRUE
  )

if (nrow(ann) == 0) {
  stop("No valid annotation rows were retained.", call. = FALSE)
}



# 4) Import UniProt reference proteome FASTA
fasta_lines <- readLines(reference_fasta, warn = FALSE)

fasta_headers <- fasta_lines[
  str_starts(fasta_lines, ">")
]

if (length(fasta_headers) == 0) {
  stop("No FASTA headers found.", call. = FALSE)
}

reference_proteome_raw <- tibble(
  fasta_header = str_remove(fasta_headers, "^>"),
  protein_id = parse_uniprot_id(fasta_headers)
) %>%
  mutate(
    protein_id = clean_character(protein_id)
  ) %>%
  filter(
    !is.na(protein_id),
    !str_detect(
      protein_id,
      regex("^(Cont_|REV_|DECOY_)", ignore_case = TRUE)
    )
  )

reference_proteome <- reference_proteome_raw %>%
  distinct(protein_id)

n_fasta_headers <- length(fasta_headers)
n_reference_proteins <- nrow(reference_proteome)
n_duplicate_fasta_ids <- nrow(reference_proteome_raw) - n_reference_proteins

if (n_reference_proteins == 0) {
  stop("No protein IDs could be parsed from the FASTA file.", call. = FALSE)
}



# 5) Restrict STRING annotations to reference FASTA proteins
ann_reference <- ann %>%
  semi_join(reference_proteome, by = "protein_id")

annotated_reference_proteins <- ann_reference %>%
  distinct(protein_id)

unannotated_reference_proteins <- reference_proteome %>%
  anti_join(annotated_reference_proteins, by = "protein_id")

string_proteins_not_in_reference <- ann %>%
  distinct(protein_id) %>%
  anti_join(reference_proteome, by = "protein_id")

n_annotated_reference_proteins <- nrow(annotated_reference_proteins)
n_unannotated_reference_proteins <- nrow(unannotated_reference_proteins)

pct_reference_annotated <-
  100 * n_annotated_reference_proteins / n_reference_proteins

if (n_annotated_reference_proteins == 0) {
  stop(
    paste0(
      "No protein IDs overlap between FASTA and STRING.\n\n",
      "Example FASTA IDs: ",
      paste(head(reference_proteome$protein_id, 5), collapse = ", "),
      "\nExample STRING IDs: ",
      paste(head(ann$protein_id, 5), collapse = ", ")
    ),
    call. = FALSE
  )
}



# 6) Per-protein statistics within annotation category


protein_annotation_by_category <- ann_reference %>%
  group_by(category, protein_id) %>%
  summarise(
    n_terms = n_distinct(term),
    terms = paste(sort(unique(term)), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(category, desc(n_terms), protein_id)



# 7) Per-term statistics within annotation category

term_annotation_by_category <- ann_reference %>%
  group_by(category, term) %>%
  summarise(
    description = first(description),
    n_proteins = n_distinct(protein_id),
    proteins = paste(sort(unique(protein_id)), collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(
    pct_reference_proteome =
      100 * n_proteins / n_reference_proteins
  ) %>%
  arrange(category, desc(n_proteins), term)



# 8) Summary statistics by category

category_coverage <- ann_reference %>%
  group_by(category) %>%
  summarise(
    annotation_pairs = n(),
    annotated_reference_proteins = n_distinct(protein_id),
    unique_terms = n_distinct(term),
    .groups = "drop"
  ) %>%
  mutate(
    pct_reference_proteome_annotated =
      100 * annotated_reference_proteins / n_reference_proteins
  )


protein_summary_by_category <- protein_annotation_by_category %>%
  group_by(category) %>%
  summarise(
    min_terms_per_protein = min(n_terms),
    q1_terms_per_protein = unname(quantile(n_terms, 0.25)),
    median_terms_per_protein = median(n_terms),
    mean_terms_per_protein = mean(n_terms),
    q3_terms_per_protein = unname(quantile(n_terms, 0.75)),
    max_terms_per_protein = max(n_terms),
    sd_terms_per_protein = sd(n_terms),
    .groups = "drop"
  )


term_summary_by_category <- term_annotation_by_category %>%
  group_by(category) %>%
  summarise(
    min_proteins_per_term = min(n_proteins),
    q1_proteins_per_term = unname(quantile(n_proteins, 0.25)),
    median_proteins_per_term = median(n_proteins),
    mean_proteins_per_term = mean(n_proteins),
    q3_proteins_per_term = unname(quantile(n_proteins, 0.75)),
    max_proteins_per_term = max(n_proteins),
    sd_proteins_per_term = sd(n_proteins),
    .groups = "drop"
  )


annotation_statistics_by_category <- category_coverage %>%
  left_join(protein_summary_by_category, by = "category") %>%
  left_join(term_summary_by_category, by = "category") %>%
  arrange(desc(pct_reference_proteome_annotated))



# 9) Global statistics and identifier matching

global_annotation_statistics <- tibble(
  metric = c(
    "FASTA header entries",
    "Unique reference FASTA proteins",
    "Duplicate FASTA protein IDs",
    "STRING annotation pairs in full term file",
    "STRING annotated proteins in full term file",
    "Reference proteins with >=1 STRING annotation",
    "Reference proteins without STRING annotation",
    "Percent reference proteome annotated",
    "STRING annotation pairs within reference proteome",
    "STRING categories within reference proteome",
    "STRING terms within reference proteome"
  ),
  value = c(
    n_fasta_headers,
    n_reference_proteins,
    n_duplicate_fasta_ids,
    nrow(ann),
    n_distinct(ann$protein_id),
    n_annotated_reference_proteins,
    n_unannotated_reference_proteins,
    pct_reference_annotated,
    nrow(ann_reference),
    n_distinct(ann_reference$category),
    n_distinct(ann_reference$term)
  )
)

identifier_overlap_diagnostic <- tibble(
  metric = c(
    "Protein IDs in reference FASTA",
    "Protein IDs in STRING term file",
    "Shared FASTA and STRING protein IDs",
    "FASTA proteins without STRING annotation",
    "STRING proteins absent from FASTA"
  ),
  value = c(
    n_reference_proteins,
    n_distinct(ann$protein_id),
    n_annotated_reference_proteins,
    n_unannotated_reference_proteins,
    nrow(string_proteins_not_in_reference)
  )
)



# 10) Export tables

write_tsv(
  global_annotation_statistics,
  file.path(out_dir, "global_annotation_statistics.tsv")
)

write_tsv(
  identifier_overlap_diagnostic,
  file.path(out_dir, "identifier_overlap_diagnostic.tsv")
)

write_tsv(
  annotation_statistics_by_category,
  file.path(out_dir, "annotation_statistics_by_category.tsv")
)

write_tsv(
  protein_annotation_by_category,
  file.path(out_dir, "protein_annotation_statistics_by_category.tsv")
)

write_tsv(
  term_annotation_by_category,
  file.path(out_dir, "term_annotation_statistics_by_category.tsv")
)

write_tsv(
  unannotated_reference_proteins,
  file.path(out_dir, "unannotated_reference_proteins.tsv")
)

write_tsv(
  string_proteins_not_in_reference,
  file.path(out_dir, "string_proteins_not_in_reference_fasta.tsv")
)



# 11) Plot data
category_plot <- annotation_statistics_by_category %>%
  mutate(
    category_label = str_wrap(category, width = 30)
  ) %>%
  arrange(pct_reference_proteome_annotated) %>%
  mutate(
    category_label = factor(
      category_label,
      levels = category_label
    )
  )

category_levels <- levels(category_plot$category_label)

protein_plot_data <- protein_annotation_by_category %>%
  mutate(
    category_label = str_wrap(category, width = 30),
    category_label = factor(
      category_label,
      levels = category_levels
    )
  )

coverage_plot_data <- tibble(
  group = c(
    "Reference FASTA proteins",
    "With STRING annotation",
    "Without STRING annotation"
  ),
  n_proteins = c(
    n_reference_proteins,
    n_annotated_reference_proteins,
    n_unannotated_reference_proteins
  )
) %>%
  mutate(
    group = factor(group, levels = rev(group))
  )


# 12) Plot A: Overall coverage
p_overall_coverage <- ggplot(
  coverage_plot_data,
  aes(x = n_proteins, y = group, fill = group)
) +
  geom_col(width = 0.68) +
  geom_text(
    aes(label = scales::comma(n_proteins)),
    hjust = -0.10,
    fontface = "bold",
    size = 4
  ) +
  scale_fill_manual(
    values = c(
      "Reference FASTA proteins" = "#2C7FB8",
      "With STRING annotation" = "#41AB5D",
      "Without STRING annotation" = "#BDBDBD"
    )
  ) +
  scale_x_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.18))
  ) +
  labs(
    title = "Reference proteome coverage",
    subtitle = "FASTA proteins with >=1 STRING annotation",
    x = "Proteins",
    y = NULL
  ) +
  theme_annotation



# 13) Plot B: Coverage by category
p_category_coverage <- ggplot(
  category_plot,
  aes(
    x = category_label,
    y = pct_reference_proteome_annotated
  )
) +
  geom_col(
    fill = "#6BAED6",
    width = 0.70
  ) +
  geom_text(
    aes(
      label = paste0(
        round(pct_reference_proteome_annotated, 1),
        "%"
      )
    ),
    hjust = -0.12,
    size = 3.1
  ) +
  coord_flip(clip = "off") +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.20))
  ) +
  labs(
    title = "Coverage by annotation category",
    subtitle = "Reference proteins with >=1 term",
    x = NULL,
    y = "Reference proteome annotated"
  ) +
  theme_annotation



# 14) Plot C: Number of terms by category
p_unique_terms <- ggplot(
  category_plot,
  aes(
    x = category_label,
    y = unique_terms
  )
) +
  geom_col(
    fill = "#9ECAE1",
    width = 0.70
  ) +
  geom_text(
    aes(label = scales::comma(unique_terms)),
    hjust = -0.12,
    size = 3.1
  ) +
  coord_flip(clip = "off") +
  scale_y_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.20))
  ) +
  labs(
    title = "Unique terms by annotation category",
    subtitle = "Terms represented in the reference proteome",
    x = NULL,
    y = "Terms"
  ) +
  theme_annotation



# 15) Plot D: Terms per protein by category

p_terms_per_protein <- ggplot(
  protein_plot_data,
  aes(
    x = n_terms,
    y = category_label
  )
) +
  geom_boxplot(
    fill = "#9ECAE1",
    width = 0.65,
    outlier.alpha = 0.18,
    outlier.size = 0.55
  ) +
  scale_x_log10(
    breaks = c(1, 2, 5, 10, 20, 50, 100, 200, 500),
    labels = scales::comma
  ) +
  labs(
    title = "Terms per annotated protein",
    subtitle = "Within annotation category",
    x = "Distinct terms per protein (log10 scale)",
    y = NULL
  ) +
  theme_annotation



# 16) Combined figure

combined_annotation_figure <- (
  (p_overall_coverage | p_category_coverage) /
    (p_unique_terms | p_terms_per_protein)
) +
  plot_annotation(
    title = "STRING annotation statistics",
    subtitle = paste0(
      "UniProt rainbow trout reference proteome: n = ",
      scales::comma(n_reference_proteins),
      " proteins"
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 10.5, color = "grey35")
    )
  )



# 17) Save plots
save_plot(
  p_overall_coverage,
  "01_reference_proteome_coverage",
  width = 7.2,
  height = 4.5
)

save_plot(
  p_category_coverage,
  "02_coverage_by_category",
  width = 8.0,
  height = 5.8
)

save_plot(
  p_unique_terms,
  "03_unique_terms_by_category",
  width = 8.0,
  height = 5.8
)

save_plot(
  p_terms_per_protein,
  "04_terms_per_protein_by_category",
  width = 8.0,
  height = 5.8
)

save_plot(
  combined_annotation_figure,
  "05_combined_annotation_statistics",
  width = 14,
  height = 10.5
)



# 18) Console output
cat("\n============================================================\n")
cat("STRING ANNOTATION STATISTICS\n")
cat("============================================================\n\n")

cat(
  "Reference FASTA proteins: ",
  scales::comma(n_reference_proteins),
  "\n",
  sep = ""
)

cat(
  "Reference proteins with >=1 STRING annotation: ",
  scales::comma(n_annotated_reference_proteins),
  " (",
  round(pct_reference_annotated, 1),
  "%)\n\n",
  sep = ""
)

print(annotation_statistics_by_category)

cat(
  "\nOutput directory:\n",
  normalizePath(out_dir),
  "\n",
  sep = ""
)

print(combined_annotation_figure)