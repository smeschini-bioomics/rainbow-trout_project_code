library(mixOmics)
library(readr)
library(tidyverse)
library(data.table)

# Load metadata
metadata <- fread("STPN2309_metadata_all.tsv")
metadata <- as.data.frame(metadata)

rownames(metadata) <- metadata[[1]]
metadata <- metadata[, -1, drop = FALSE]

# Load protein abundance matrix
proteome <- fread(
  "01_MSDAP_DEA/2026-04-25_07-16-26/protein_abundance__global data filter.tsv"
)
proteome <- as.data.frame(proteome)

# Use protein IDs as row names
rownames(proteome) <- proteome[[1]]
proteome <- proteome[, -1, drop = FALSE]

# Remove annotation columns
proteome <- proteome %>%
  dplyr::select(-fasta_headers, -gene_symbols_or_id)

# Remove contaminant proteins
proteome <- proteome[
  !stringr::str_detect(rownames(proteome), "Cont_"),
  ,
  drop = FALSE
]

# Transpose to samples × proteins
proteome <- t(as.matrix(proteome))

# Match proteome sample IDs to metadata row names
stopifnot("sampleID_proteomics_gut" %in% colnames(metadata))

sample_id_map <- setNames(
  rownames(metadata),
  metadata$sampleID_proteomics_gut
)

old_sample_ids <- rownames(proteome)
new_sample_ids <- unname(sample_id_map[old_sample_ids])

unmatched_sample_ids <- old_sample_ids[is.na(new_sample_ids)]

if (length(unmatched_sample_ids) > 0) {
  stop(
    "Proteome sample IDs not found in metadata$sampleID_proteomics_gut: ",
    paste(head(unmatched_sample_ids, 10), collapse = ", "),
    ifelse(length(unmatched_sample_ids) > 10, " ...", "")
  )
}

if (anyDuplicated(new_sample_ids)) {
  stop(
    "Duplicated metadata sample IDs after proteome renaming: ",
    paste(unique(new_sample_ids[duplicated(new_sample_ids)]), collapse = ", ")
  )
}

rownames(proteome) <- new_sample_ids

stopifnot(all(rownames(proteome) %in% rownames(metadata)))

# Remove proteins with more than 30% missing values
missing_rate <- colMeans(is.na(proteome))

plot(
  missing_rate,
  type = "h",
  xlab = "Protein index",
  ylab = "Missing-value rate",
  main = "Protein missing-value rate"
)

proteome_filtered <- proteome[
  ,
  missing_rate <= 0.30,
  drop = FALSE
]

# Check final matrix
cat(
  "Proteome matrix before filtering:",
  nrow(proteome), "samples ×", ncol(proteome), "proteins\n"
)

cat(
  "Proteome matrix after filtering:",
  nrow(proteome_filtered), "samples ×", ncol(proteome_filtered), "proteins\n"
)

cat(
  "Remaining missing values:",
  sum(is.na(proteome_filtered)),
  "out of",
  length(proteome_filtered),
  "(",
  round(100 * mean(is.na(proteome_filtered)), 2),
  "%)\n"
)

# Save filtered proteome matrix
RDS_DIR <- "rds"
dir.create(RDS_DIR, showWarnings = FALSE, recursive = TRUE)

saveRDS(
  proteome_filtered,
  file = file.path(RDS_DIR, "midgut_proteome_filtered.rds")
)
