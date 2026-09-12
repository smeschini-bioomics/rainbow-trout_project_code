# Recover the sample size used in the Python Spearman analysis

suppressPackageStartupMessages({
  library(biomformat)
  library(readr)
  library(dplyr)
})

biom_file <- "output_files/digesta_filtered_table_with_taxonomy.biom"
metadata_file <- "output_files/metadata_model_used_for_BIRDMAn.tsv"

physio_cols <- c(
  "Body_Weight",
  "Plasma_Glucose",
  "Plasma_L-Lactate",
  "Plasma_D-Lactate",
  "HSI"
)

# Read sample IDs from an HDF5 BIOM file

if (!requireNamespace("rhdf5", quietly = TRUE)) {
  stop(
    "Package 'rhdf5' is required. Install it with:\n",
    "BiocManager::install('rhdf5')"
  )
}

# Optional: inspect the HDF5 structure
print(
  rhdf5::h5ls(biom_file)
)

# BIOM 2.x sample IDs are normally stored here
biom_sample_ids <- rhdf5::h5read(
  biom_file,
  "sample/ids"
)

biom_sample_ids <- as.character(biom_sample_ids)

cat(
  "Number of samples found in BIOM:",
  length(biom_sample_ids),
  "\n"
)
# Load metadata, using the first column as sample ID
metadata <- read.delim(
  metadata_file,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

# Apply the same sample filtering as the Python script
metadata_model <- metadata %>%
  filter(
    sample_type == "Digesta",
    Diet %in% c("NC", "HC")
  )

# Same intersection performed in Python
common_samples <- intersect(
  biom_sample_ids,
  rownames(metadata_model)
)

meta_common <- metadata_model[
  common_samples,
  physio_cols,
  drop = FALSE
]

# Same complete-case filter used in Python
complete_mask <- complete.cases(meta_common)

n_common <- length(common_samples)
n_used <- sum(complete_mask)
n_removed <- n_common - n_used

cat(
  "Samples shared between BIOM and metadata:", n_common, "\n",
  "Samples used for all Spearman correlations:", n_used, "\n",
  "Samples removed because at least one trait was missing:", n_removed, "\n"
)

# Export summary
sample_size_summary <- tibble(
  tissue = "Digesta",
  n_shared_samples = n_common,
  n_complete_samples_used = n_used,
  n_removed_missing_physiology = n_removed
)

write_tsv(
  sample_size_summary,
  file.path(
    out_dir,
    "digesta_spearman_CLR_sample_size_summary.tsv"
  )
)

# Export exact sample IDs used
samples_used <- tibble(
  sample_id = common_samples[complete_mask]
)

write_tsv(
  samples_used,
  file.path(
    out_dir,
    "digesta_spearman_CLR_samples_used.tsv"
  )
)
