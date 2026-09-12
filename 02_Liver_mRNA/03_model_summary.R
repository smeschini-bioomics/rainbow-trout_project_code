
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

genes <- list.files("input_files", pattern = "\\.tsv$") |>
  tools::file_path_sans_ext()

summary_all <- lapply(genes, function(gene) {

  model_file <- file.path("results", gene, paste0(gene, "_model_summary.txt"))

  # No model-summary file = no inferential analysis was completed
  if (!file.exists(model_file)) {
    return(tibble(
      gene = gene,
      statistical_analysis = NA_character_,
      model_family = NA_character_,
      model_link = NA_character_
    ))
  }

  model_text <- trimws(readLines(model_file, warn = FALSE))
  family_line <- model_text[str_detect(model_text, "Family:")][1]

  # Typical glmmTMB summary line: Family: t  ( identity )
  model_family <- str_trim(str_remove(
    str_split(family_line, "\\(")[[1]][1],
    "^Family:"
  ))

  model_link <- str_extract(family_line, "\\([^)]*\\)") |>
    str_remove_all("[()]") |>
    str_trim()

  tibble(
    gene = gene,
    statistical_analysis = "completed",
    model_family = model_family,
    model_link = model_link
  )
}) |>
  bind_rows() |>
  arrange(gene)

dir.create("results/summary_all_genes", recursive = TRUE, showWarnings = FALSE)

write_csv(
  summary_all,
  "results/summary_all_genes/all_genes_model_summary.csv"
)

print(summary_all)

