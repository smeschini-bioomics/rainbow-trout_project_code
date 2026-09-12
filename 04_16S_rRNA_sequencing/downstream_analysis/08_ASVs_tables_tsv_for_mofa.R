library(microeco)
library(magrittr)
library(dplyr)

# SETTINGS
outdir <- "08_ASVs_tables_tsv_for_mofa"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

freq_filter <- 0.20
rel_abund_filter <- 0

bad_tax_values <- c(
  "", NA,
  "k__", "p__", "c__", "o__", "f__", "g__", "s__",
  "unclassified", "uncultured", "unknown", "metagenome",
  "NA"
)

# DIGESTA - ASV LEVEL

load("RData/16s_dataset_digesta.RData")

# Work on a copy, not original object
mt_dig_asv <- mt_dig$clone(deep = TRUE)
mt_dig_asv$tidy_dataset()

cat("\n================ DIGESTA: ORIGINAL ================\n")
print(mt_dig_asv)
print(range(mt_dig_asv$sample_sums()))
cat("Taxa abundance range:\n")
print(range(mt_dig_asv$taxa_sums()))

# Sanity check: ASV IDs must match taxonomy rows
cat("\nDIGESTA: ASV-taxonomy matching before filtering\n")
print(all(rownames(mt_dig_asv$otu_table) %in% rownames(mt_dig_asv$tax_table)))
stopifnot(all(rownames(mt_dig_asv$otu_table) %in% rownames(mt_dig_asv$tax_table)))

# Filter ASVs
mt_dig_asv$filter_taxa(
  rel_abund = rel_abund_filter,
  freq = freq_filter,
  include_lowest = TRUE,
  for_taxa_abund = FALSE
)

mt_dig_asv$tidy_dataset()

cat("\n================ DIGESTA: AFTER FILTERING ================\n")
print(mt_dig_asv)
cat("Sample depth range:\n")
print(range(mt_dig_asv$sample_sums()))
cat("Taxa abundance range:\n")
print(range(mt_dig_asv$taxa_sums()))

# CLR normalization
tmp <- trans_norm$new(dataset = mt_dig_asv)
mt_dig_asv_clr <- tmp$norm(method = "clr")
mt_dig_asv_clr$tidy_dataset()

cat("\n================ DIGESTA: AFTER CLR ================\n")
print(mt_dig_asv_clr)

# ASV abundance table for your current TSV pipeline
# rows = samples, columns = ASVs
otu_dig_asv <- t(as.matrix(mt_dig_asv_clr$otu_table))

cat("\nDIGESTA: exported abundance orientation\n")
cat("Rows should be samples, columns should be ASVs\n")
print(dim(otu_dig_asv))
print(head(rownames(otu_dig_asv)))
print(head(colnames(otu_dig_asv)))

# Taxonomy table
tax_dig_asv <- as.data.frame(mt_dig_asv_clr$tax_table)
tax_dig_asv$ASV_ID <- rownames(tax_dig_asv)

# Sanity check: abundance ASVs and taxonomy ASVs
cat("\nDIGESTA: ASV-taxonomy matching after CLR\n")
print(all(colnames(otu_dig_asv) %in% tax_dig_asv$ASV_ID))
print(all(tax_dig_asv$ASV_ID %in% colnames(otu_dig_asv)))

stopifnot(all(colnames(otu_dig_asv) %in% tax_dig_asv$ASV_ID))
stopifnot(all(tax_dig_asv$ASV_ID %in% colnames(otu_dig_asv)))

# Reorder taxonomy to match abundance table columns
tax_dig_asv <- tax_dig_asv[colnames(otu_dig_asv), , drop = FALSE]

cat("\nDIGESTA: after reordering taxonomy\n")
print(identical(tax_dig_asv$ASV_ID, colnames(otu_dig_asv)))

stopifnot(identical(tax_dig_asv$ASV_ID, colnames(otu_dig_asv)))

# Create best available taxonomy label
tax_ranks_dig <- c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom")
tax_ranks_dig <- tax_ranks_dig[tax_ranks_dig %in% colnames(tax_dig_asv)]

cat("\nDIGESTA: taxonomy ranks used for best_tax\n")
print(tax_ranks_dig)

tax_dig_asv$best_tax <- apply(
  tax_dig_asv[, tax_ranks_dig, drop = FALSE],
  1,
  function(x) {
    x <- as.character(x)
    x <- x[!is.na(x)]
    x <- trimws(x)
    x <- x[x != ""]
    x <- x[!x %in% c("k__", "p__", "c__", "o__", "f__", "g__", "s__")]
    x <- x[!grepl("unclassified|uncultured|unknown|metagenome", x, ignore.case = TRUE)]
    
    if (length(x) == 0) {
      return("Unclassified")
    } else {
      return(x[1])
    }
  }
)

# Create ASV_best label
tax_dig_asv$ASV_best <- paste(tax_dig_asv$best_tax, tax_dig_asv$ASV_ID, sep = "|")

# Final ASV-taxonomy association matrix
cols_keep_dig <- c(
  "ASV_ID",
  "ASV_best",
  "best_tax",
  "Kingdom",
  "Phylum",
  "Class",
  "Order",
  "Family",
  "Genus",
  "Species"
)

cols_keep_dig <- cols_keep_dig[cols_keep_dig %in% colnames(tax_dig_asv)]

digesta_ASV_best_map <- tax_dig_asv[, cols_keep_dig, drop = FALSE]

cat("\n================ DIGESTA: FINAL ASV MAP ================\n")
print(head(digesta_ASV_best_map))
cat("ASV map dim:\n")
print(dim(digesta_ASV_best_map))

# Sanity checks for final map
cat("\nDIGESTA: final sanity checks\n")
cat("nrow map == ncol abundance:\n")
print(nrow(digesta_ASV_best_map) == ncol(otu_dig_asv))

cat("ASV_ID identical to abundance columns:\n")
print(identical(digesta_ASV_best_map$ASV_ID, colnames(otu_dig_asv)))

cat("Duplicated ASV_ID:\n")
print(sum(duplicated(digesta_ASV_best_map$ASV_ID)))

cat("Duplicated ASV_best:\n")
print(sum(duplicated(digesta_ASV_best_map$ASV_best)))

stopifnot(nrow(digesta_ASV_best_map) == ncol(otu_dig_asv))
stopifnot(identical(digesta_ASV_best_map$ASV_ID, colnames(otu_dig_asv)))
stopifnot(sum(duplicated(digesta_ASV_best_map$ASV_ID)) == 0)
stopifnot(sum(duplicated(digesta_ASV_best_map$ASV_best)) == 0)

# Human-readable abundance table with ASV_best column names
# Do NOT use this one for training unless you really want long feature names.
otu_dig_asv_best <- otu_dig_asv
colnames(otu_dig_asv_best) <- digesta_ASV_best_map$ASV_best

cat("\nDIGESTA: ASV_best abundance sanity check\n")
print(identical(colnames(otu_dig_asv_best), digesta_ASV_best_map$ASV_best))

# Export files
write.table(
  otu_dig_asv,
  file = file.path(outdir, "digesta_ASV_table_mofa_f20_ra0.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

write.table(
  digesta_ASV_best_map,
  file = file.path(outdir, "digesta_ASV_best_taxonomy_map.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

write.table(
  otu_dig_asv_best,
  file = file.path(outdir, "digesta_ASV_table_with_ASVbest_names.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

cat("\nDIGESTA files saved.\n")


# MUCUS - ASV LEVEL

load("RData/16s_dataset_mucus.RData")

# Work on a copy, not original object
mt_mucus_asv <- mt_mucus$clone(deep = TRUE)
mt_mucus_asv$tidy_dataset()

cat("\n================ MUCUS: ORIGINAL ================\n")
print(mt_mucus_asv)
cat("OTU table dim:\n")
print(dim(mt_mucus_asv$otu_table))
cat("Tax table dim:\n")
print(dim(mt_mucus_asv$tax_table))
cat("Sample depth range:\n")
print(range(mt_mucus_asv$sample_sums()))
cat("Taxa abundance range:\n")
print(range(mt_mucus_asv$taxa_sums()))

# Sanity check: ASV IDs must match taxonomy rows
cat("\nMUCUS: ASV-taxonomy matching before filtering\n")
print(all(rownames(mt_mucus_asv$otu_table) %in% rownames(mt_mucus_asv$tax_table)))
print(all(rownames(mt_mucus_asv$tax_table) %in% rownames(mt_mucus_asv$otu_table)))

stopifnot(all(rownames(mt_mucus_asv$otu_table) %in% rownames(mt_mucus_asv$tax_table)))

# Filter ASVs
mt_mucus_asv$filter_taxa(
  rel_abund = rel_abund_filter,
  freq = freq_filter,
  include_lowest = TRUE,
  for_taxa_abund = FALSE
)

mt_mucus_asv$tidy_dataset()

cat("\n================ MUCUS: AFTER FILTERING ================\n")
print(mt_mucus_asv)
cat("OTU table dim:\n")
print(dim(mt_mucus_asv$otu_table))
cat("Tax table dim:\n")
print(dim(mt_mucus_asv$tax_table))
cat("Sample depth range:\n")
print(range(mt_mucus_asv$sample_sums()))
cat("Taxa abundance range:\n")
print(range(mt_mucus_asv$taxa_sums()))

# CLR normalization
tmp <- trans_norm$new(dataset = mt_mucus_asv)
mt_mucus_asv_clr <- tmp$norm(method = "clr")
mt_mucus_asv_clr$tidy_dataset()

cat("\n================ MUCUS: AFTER CLR ================\n")
cat("CLR OTU table dim:\n")
print(dim(mt_mucus_asv_clr$otu_table))
cat("CLR tax table dim:\n")
print(dim(mt_mucus_asv_clr$tax_table))

# ASV abundance table for your current TSV pipeline
# rows = samples, columns = ASVs
otu_mucus_asv <- t(as.matrix(mt_mucus_asv_clr$otu_table))

cat("\nMUCUS: exported abundance orientation\n")
cat("Rows should be samples, columns should be ASVs\n")
print(dim(otu_mucus_asv))
print(head(rownames(otu_mucus_asv)))
print(head(colnames(otu_mucus_asv)))

# Taxonomy table
tax_mucus_asv <- as.data.frame(mt_mucus_asv_clr$tax_table)
tax_mucus_asv$ASV_ID <- rownames(tax_mucus_asv)

# Sanity check: abundance ASVs and taxonomy ASVs
cat("\nMUCUS: ASV-taxonomy matching after CLR\n")
print(all(colnames(otu_mucus_asv) %in% tax_mucus_asv$ASV_ID))
print(all(tax_mucus_asv$ASV_ID %in% colnames(otu_mucus_asv)))

stopifnot(all(colnames(otu_mucus_asv) %in% tax_mucus_asv$ASV_ID))
stopifnot(all(tax_mucus_asv$ASV_ID %in% colnames(otu_mucus_asv)))

# Reorder taxonomy to match abundance table columns
tax_mucus_asv <- tax_mucus_asv[colnames(otu_mucus_asv), , drop = FALSE]

cat("\nMUCUS: after reordering taxonomy\n")
print(identical(tax_mucus_asv$ASV_ID, colnames(otu_mucus_asv)))

stopifnot(identical(tax_mucus_asv$ASV_ID, colnames(otu_mucus_asv)))

# Create best available taxonomy label
tax_ranks_mucus <- c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom")
tax_ranks_mucus <- tax_ranks_mucus[tax_ranks_mucus %in% colnames(tax_mucus_asv)]

cat("\nMUCUS: taxonomy ranks used for best_tax\n")
print(tax_ranks_mucus)

tax_mucus_asv$best_tax <- apply(
  tax_mucus_asv[, tax_ranks_mucus, drop = FALSE],
  1,
  function(x) {
    x <- as.character(x)
    x <- x[!is.na(x)]
    x <- trimws(x)
    x <- x[x != ""]
    x <- x[!x %in% c("k__", "p__", "c__", "o__", "f__", "g__", "s__")]
    x <- x[!grepl("unclassified|uncultured|unknown|metagenome", x, ignore.case = TRUE)]
    
    if (length(x) == 0) {
      return("Unclassified")
    } else {
      return(x[1])
    }
  }
)

# Create ASV_best label
tax_mucus_asv$ASV_best <- paste(tax_mucus_asv$best_tax, tax_mucus_asv$ASV_ID, sep = "|")

# Final ASV-taxonomy association matrix
cols_keep_mucus <- c(
  "ASV_ID",
  "ASV_best",
  "best_tax",
  "Kingdom",
  "Phylum",
  "Class",
  "Order",
  "Family",
  "Genus",
  "Species"
)

cols_keep_mucus <- cols_keep_mucus[cols_keep_mucus %in% colnames(tax_mucus_asv)]

mucus_ASV_best_map <- tax_mucus_asv[, cols_keep_mucus, drop = FALSE]

cat("\n================ MUCUS: FINAL ASV MAP ================\n")
print(head(mucus_ASV_best_map))
cat("ASV map dim:\n")
print(dim(mucus_ASV_best_map))

# Sanity checks for final map
cat("\nMUCUS: final sanity checks\n")
cat("nrow map == ncol abundance:\n")
print(nrow(mucus_ASV_best_map) == ncol(otu_mucus_asv))

cat("ASV_ID identical to abundance columns:\n")
print(identical(mucus_ASV_best_map$ASV_ID, colnames(otu_mucus_asv)))

cat("Duplicated ASV_ID:\n")
print(sum(duplicated(mucus_ASV_best_map$ASV_ID)))

cat("Duplicated ASV_best:\n")
print(sum(duplicated(mucus_ASV_best_map$ASV_best)))

stopifnot(nrow(mucus_ASV_best_map) == ncol(otu_mucus_asv))
stopifnot(identical(mucus_ASV_best_map$ASV_ID, colnames(otu_mucus_asv)))
stopifnot(sum(duplicated(mucus_ASV_best_map$ASV_ID)) == 0)
stopifnot(sum(duplicated(mucus_ASV_best_map$ASV_best)) == 0)

otu_mucus_asv_best <- otu_mucus_asv
colnames(otu_mucus_asv_best) <- mucus_ASV_best_map$ASV_best

cat("\nMUCUS: ASV_best abundance sanity check\n")
print(identical(colnames(otu_mucus_asv_best), mucus_ASV_best_map$ASV_best))

# Export files
write.table(
  otu_mucus_asv,
  file = file.path(outdir, "mucus_ASV_table_mofa_f20_ra0.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

write.table(
  mucus_ASV_best_map,
  file = file.path(outdir, "mucus_ASV_best_taxonomy_map.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

write.table(
  otu_mucus_asv_best,
  file = file.path(outdir, "mucus_ASV_table_with_ASVbest_names.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

cat("\nMUCUS files saved.\n")

# Record package and R versions
writeLines(
  capture.output(sessionInfo()),
  con = file.path(outdir, "sessionInfo.txt")
)

