# =============================================================================
# Liver histone PTM matrix processing
# =============================================================================

library(data.table)
library(stringr)

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------

METADATA_FILE <- "STPN2309_metadata_all.tsv"
INPUT_FILE    <- "01_output_files/TroutLiver_hPTM_H3H4.csv"
OUTDIR        <- "02_process_matrix"

MISSINGNESS_THRESHOLD <- 0.30
COLLAPSE_TOLERANCE    <- 6


dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Helper: collapse columns with identical signal profiles
# -----------------------------------------------------------------------------

collapse_identical_cols <- function(X, sep = "|", tol = NULL) {
  
  X <- as.matrix(X)
  X_key <- if (is.null(tol)) X else round(X, tol)
  
  column_key <- apply(
    X_key,
    2,
    function(x) {
      paste(
        ifelse(
          is.na(x),
          "NA",
          format(x, scientific = FALSE, digits = 15, trim = TRUE)
        ),
        collapse = "\u241F"
      )
    }
  )
  
  keep_columns <- !duplicated(column_key)
  X_collapsed <- X[, keep_columns, drop = FALSE]
  
  groups <- split(colnames(X), column_key)
  
  new_names <- vapply(
    colnames(X_collapsed),
    function(column_name) {
      key <- column_key[match(column_name, colnames(X))]
      paste(groups[[key]], collapse = sep)
    },
    character(1)
  )
  
  colnames(X_collapsed) <- new_names
  
  X_collapsed
}

# -----------------------------------------------------------------------------
# Import metadata
# -----------------------------------------------------------------------------

metadata <- fread(METADATA_FILE)
metadata <- as.data.frame(metadata)

rownames(metadata) <- metadata[[1]]
metadata <- metadata[, -1, drop = FALSE]

stopifnot("sampleID_hPTMs_liver" %in% colnames(metadata))

cat(
  "Metadata:",
  nrow(metadata),
  "samples x",
  ncol(metadata),
  "variables\n"
)

# -----------------------------------------------------------------------------
# Import liver hPTM matrix
# -----------------------------------------------------------------------------

histones_PTM <- read.csv(
  INPUT_FILE,
  sep = ",",
  check.names = FALSE
)

histones_PTM <- as.data.frame(histones_PTM)

rownames(histones_PTM) <- histones_PTM[[1]]
histones_PTM <- histones_PTM[, -1, drop = FALSE]

histones_PTM <- as.matrix(histones_PTM)
histones_PTM <- t(histones_PTM)

cat(
  "Raw hPTM matrix:",
  nrow(histones_PTM),
  "samples x",
  ncol(histones_PTM),
  "features\n"
)

# -----------------------------------------------------------------------------
# Standardize PTM feature names
# Example: H3 [9] (K) Ac -> H3K9Ac
# -----------------------------------------------------------------------------

original_names <- colnames(histones_PTM)

name_match <- str_match(
  original_names,
  "^(H[0-9AB]+)\\s+\\[(\\d+)\\]\\s+\\(([A-Z])\\)\\s+([A-Za-z0-9]+)$"
)

standardized_names <- ifelse(
  is.na(name_match[, 1]),
  original_names,
  paste0(
    name_match[, 2],
    name_match[, 4],
    name_match[, 3],
    name_match[, 5]
  )
)

if (anyDuplicated(standardized_names)) {
  standardized_names <- make.unique(standardized_names)
}

colnames(histones_PTM) <- standardized_names

# -----------------------------------------------------------------------------
# Select H3 and H4 PTMs
# -----------------------------------------------------------------------------

is_H3_H4 <- grepl(
  "^H[34]",
  colnames(histones_PTM),
  ignore.case = TRUE
)

H3_H4_PTM <- histones_PTM[, is_H3_H4, drop = FALSE]

cat("H3/H4 PTM columns:", ncol(H3_H4_PTM), "\n")

# -----------------------------------------------------------------------------
# Select PTM types of interest
# -----------------------------------------------------------------------------

selected_modifications <- c(
  "Ac",
  "Me",
  "Cr",
  "La",
  "Bu",
  "Prop"
)

modification_pattern <- paste0(
  "(",
  paste(selected_modifications, collapse = "|"),
  ")"
)

is_selected_modification <- grepl(
  modification_pattern,
  colnames(H3_H4_PTM),
  ignore.case = TRUE
)

H3_H4_selected <- H3_H4_PTM[
  ,
  is_selected_modification,
  drop = FALSE
]

cat(
  "Selected H3/H4 PTM columns before collapse:",
  ncol(H3_H4_selected),
  "\n"
)

# -----------------------------------------------------------------------------
# Collapse PTMs with identical signal profiles
# -----------------------------------------------------------------------------

H3_H4_selected_collapsed <- collapse_identical_cols(
  H3_H4_selected,
  sep = "|",
  tol = COLLAPSE_TOLERANCE
)

cat(
  "Selected H3/H4 PTM columns after collapse:",
  ncol(H3_H4_selected_collapsed),
  "\n"
)

# -----------------------------------------------------------------------------
# Filter PTMs with more than 30% missing values
# -----------------------------------------------------------------------------

na_rate <- colMeans(is.na(H3_H4_selected_collapsed))

removed_for_missingness <- names(na_rate)[
  na_rate > MISSINGNESS_THRESHOLD
]

H3_H4_selected_filtered <- H3_H4_selected_collapsed[
  ,
  na_rate <= MISSINGNESS_THRESHOLD,
  drop = FALSE
]

cat(
  "PTM columns removed for missingness:",
  length(removed_for_missingness),
  "\n"
)

cat(
  "PTM columns retained after missingness filter:",
  ncol(H3_H4_selected_filtered),
  "\n"
)

cat(
  "Overall remaining missing-value proportion:",
  round(mean(is.na(H3_H4_selected_filtered)), 4),
  "\n"
)

# -----------------------------------------------------------------------------
# Rename samples using metadata IDs
# -----------------------------------------------------------------------------

sample_id_map <- setNames(
  rownames(metadata),
  metadata$sampleID_hPTMs_liver
)

old_sample_ids <- rownames(H3_H4_selected_filtered)

new_sample_ids <- unname(
  sample_id_map[old_sample_ids]
)

unmatched_ids <- old_sample_ids[is.na(new_sample_ids)]

if (length(unmatched_ids) > 0) {
  stop(
    "hPTM sample IDs not found in metadata$sampleID_hPTMs_liver: ",
    paste(head(unmatched_ids, 10), collapse = ", "),
    ifelse(length(unmatched_ids) > 10, " ...", "")
  )
}

if (anyDuplicated(new_sample_ids)) {
  
  duplicated_ids <- unique(
    new_sample_ids[duplicated(new_sample_ids)]
  )
  
  stop(
    "Duplicate sample IDs after metadata renaming: ",
    paste(duplicated_ids, collapse = ", ")
  )
}

H3_H4_8PTM_liver <- H3_H4_selected_filtered
rownames(H3_H4_8PTM_liver) <- new_sample_ids

cat(
  "Renamed samples:",
  length(new_sample_ids),
  "/",
  length(old_sample_ids),
  "\n"
)

# -----------------------------------------------------------------------------
# Final checks
# -----------------------------------------------------------------------------

stopifnot(
  nrow(H3_H4_8PTM_liver) == nrow(H3_H4_selected_filtered),
  !anyDuplicated(rownames(H3_H4_8PTM_liver)),
  ncol(H3_H4_8PTM_liver) > 0
)

cat(
  "Final liver hPTM matrix:",
  nrow(H3_H4_8PTM_liver),
  "samples x",
  ncol(H3_H4_8PTM_liver),
  "PTM features\n"
)

print(colnames(H3_H4_8PTM_liver))

# -----------------------------------------------------------------------------
# Export
# -----------------------------------------------------------------------------

saveRDS(
  H3_H4_8PTM_liver,
  file = file.path(OUTDIR, "H3_H4_8PTM_liver.rds")
)

write.table(
  H3_H4_8PTM_liver,
  file = file.path(OUTDIR, "H3_H4_8PTM_liver.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(OUTDIR, "sessionInfo.txt")
)
