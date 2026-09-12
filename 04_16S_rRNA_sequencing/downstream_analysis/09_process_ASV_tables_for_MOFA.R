library(mixOmics)
library(readr)
library(tidyverse)
library(data.table)

# Output directory
out_dir <- "09_process_ASV_tables_for_MOFA"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# import metadata 
metadata <- fread("STPN2309_metadata_all.tsv") # metadata
dim(metadata)

# Set the first column as rownames
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata[[1]]
metadata <- metadata[, -1]
# check
class(metadata)
dim(metadata)
rownames(metadata)
str(metadata)

# Convert factors to the correct format
metadata$day_categ <- as.factor(metadata$day_categ)
metadata$rearing_tank <- as.factor(metadata$rearing_tank)
metadata$diet <- as.factor(metadata$diet)
metadata$group <- as.factor(metadata$group)
metadata$diet_phase <- as.factor(metadata$diet_phase)
# check the data structure
str(metadata)


# import genus table for multiomics (MOFA)
digesta_asv_table_mofa <- fread("08_ASVs_tables_tsv_for_mofa/digesta_ASV_table_with_ASVbest_names.tsv")
dim(digesta_asv_table_mofa)
# Set the first column as rownames
digesta_asv_table_mofa <- as.data.frame(digesta_asv_table_mofa)
rownames(digesta_asv_table_mofa) <- digesta_asv_table_mofa[[1]]
digesta_asv_table_mofa <- digesta_asv_table_mofa[, -1]
class(digesta_asv_table_mofa)
dim(digesta_asv_table_mofa)
rownames(digesta_asv_table_mofa)
digesta_asv_table_mofa <- as.matrix(digesta_asv_table_mofa)
dim(digesta_asv_table_mofa)
class(digesta_asv_table_mofa)

# rename sampleID
stopifnot("SampleID_digesta_microbiota_PGTB" %in% colnames(metadata))

# 1) Build mapping: SampleID_digesta_microbiota_PGTB -> desired new name (metadata rowname)
map_old_to_new <- setNames(rownames(metadata), metadata$SampleID_digesta_microbiota_PGTB)

# 2) Create new rownames for asv_mat_t using the mapping
old <- rownames(digesta_asv_table_mofa)
new <- unname(map_old_to_new[old])  # will be NA if not found

# 3) Report any unmatched or duplicate targets (diagnostics)
unmatched <- old[is.na(new)]
if (length(unmatched)) {
  message("IDs in ASV_table_mefisto_labeled not found in data$metadata$SampleID_digesta_microbiota_PGTB: ",
          paste(head(unmatched, 10), collapse = ", "),
          ifelse(length(unmatched) > 10, " ...", ""))
}

# if you want to **keep original names** when not matched:
new_fallback <- ifelse(is.na(new), old, new)

# detect potential duplicates after renaming
dups <- new_fallback[duplicated(new_fallback)]
if (length(dups)) {
  warning("Duplicate target rownames after renaming: ",
          paste(unique(dups), collapse = ", "))
}

# 4) Apply renaming (no subsetting, metadata remains untouched)
digesta_asv_table_mofa_renamed <- digesta_asv_table_mofa
rownames(digesta_asv_table_mofa_renamed) <- new_fallback

# (optional) quick check
cat("Renamed rows:", sum(!is.na(new)), " / ", length(old), "\n")
rownames(digesta_asv_table_mofa_renamed)



# import genus table for multiomics (MOFA)
mucus_asv_table_mofa <- fread("08_ASVs_tables_tsv_for_mofa/mucus_ASV_table_with_ASVbest_names.tsv")
dim(mucus_asv_table_mofa)
# Set the first column as rownames
mucus_asv_table_mofa <- as.data.frame(mucus_asv_table_mofa)
rownames(mucus_asv_table_mofa) <- mucus_asv_table_mofa[[1]]
mucus_asv_table_mofa <- mucus_asv_table_mofa[, -1]
class(mucus_asv_table_mofa)
dim(mucus_asv_table_mofa)
rownames(mucus_asv_table_mofa)
mucus_asv_table_mofa <- as.matrix(mucus_asv_table_mofa)
dim(mucus_asv_table_mofa)
class(mucus_asv_table_mofa)

# rename sampleID
stopifnot("SampleID_mucus_microbiota_PGTB" %in% colnames(metadata))

# 1) Build mapping: SampleID_digesta_microbiota_PGTB -> desired new name (metadata rowname)
map_old_to_new <- setNames(rownames(metadata), metadata$SampleID_mucus_microbiota_PGTB)

# 2) Create new rownames for asv_mat_t using the mapping
old <- rownames(mucus_asv_table_mofa)
new <- unname(map_old_to_new[old])  # will be NA if not found

# 3) Report any unmatched or duplicate targets (diagnostics)
unmatched <- old[is.na(new)]
if (length(unmatched)) {
  message("IDs in ASV_table_mefisto_labeled not found in data$metadata$SampleID_digesta_microbiota_PGTB: ",
          paste(head(unmatched, 10), collapse = ", "),
          ifelse(length(unmatched) > 10, " ...", ""))
}

# if you want to **keep original names** when not matched:
new_fallback <- ifelse(is.na(new), old, new)

# detect potential duplicates after renaming
dups <- new_fallback[duplicated(new_fallback)]
if (length(dups)) {
  warning("Duplicate target rownames after renaming: ",
          paste(unique(dups), collapse = ", "))
}

# 4) Apply renaming (no subsetting, metadata remains untouched)
mucus_asv_table_mofa_renamed <- mucus_asv_table_mofa
rownames(mucus_asv_table_mofa_renamed) <- new_fallback

# (optional) quick check
cat("Renamed rows:", sum(!is.na(new)), " / ", length(old), "\n")
rownames(mucus_asv_table_mofa_renamed)
class(mucus_asv_table_mofa_renamed)


# Save ASV-best microbiome tables as RDS
# Sanity checks
stopifnot(is.matrix(digesta_asv_table_mofa_renamed))
stopifnot(is.matrix(mucus_asv_table_mofa_renamed))

stopifnot(all(rownames(digesta_asv_table_mofa_renamed) %in% rownames(metadata)))
stopifnot(all(rownames(mucus_asv_table_mofa_renamed) %in% rownames(metadata)))

stopifnot(sum(duplicated(rownames(digesta_asv_table_mofa_renamed))) == 0)
stopifnot(sum(duplicated(rownames(mucus_asv_table_mofa_renamed))) == 0)

stopifnot(sum(duplicated(colnames(digesta_asv_table_mofa_renamed))) == 0)
stopifnot(sum(duplicated(colnames(mucus_asv_table_mofa_renamed))) == 0)


# Save sample x feature matrices
saveRDS(
  digesta_asv_table_mofa_renamed,
  file = file.path(
    out_dir,
    "digesta_ASVbest_sample_by_feature.rds"
  )
)

saveRDS(
  mucus_asv_table_mofa_renamed,
  file = file.path(
    out_dir,
    "mucus_ASVbest_sample_by_feature.rds"
  )
)


writeLines(
  capture.output(sessionInfo()),
  con = file.path(out_dir, "sessionInfo.txt")
)

cat("Saved MOFA ASV inputs in:", out_dir, "\n")

