# Load libraries
library(microeco)
library(ape)
library(sessioninfo)
library(biomformat)

# Input and output paths
input_file <- file.path("RData", "16s_dataset.RData")
output_dir <- "05_export_files_for_BIRDMAn_DA_analysis"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load the data
load(input_file)

# Extract ASV table
otu_df <- as.data.frame(mt_raw$otu_table)

write.table(
  otu_df,
  file = file.path(output_dir, "birdman_otu_table.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# Extract taxonomy table
tax_df <- as.data.frame(mt_raw$tax_table)

tax_df$ASV <- rownames(tax_df)
tax_df <- tax_df[, c("ASV", setdiff(colnames(tax_df), "ASV"))]

write.table(
  tax_df,
  file = file.path(output_dir, "birdman_taxonomy.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Export phylogenetic tree
tree <- mt_raw$phylo_tree

ape::write.tree(
  tree,
  file = file.path(output_dir, "tree.nwk")
)

# Record package and R versions
writeLines(
  capture.output(sessionInfo()),
  con = file.path(output_dir, "sessionInfo.txt")
)

