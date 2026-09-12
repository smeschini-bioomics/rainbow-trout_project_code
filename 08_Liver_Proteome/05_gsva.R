# libraries
library(GSVA)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(patchwork)
library(ggplot2)
library(tibble)

# Input and output paths
PROTEOME_RDS <- file.path("rds", "liver_proteome_filtered.rds")
METADATA_FILE <- "STPN2309_metadata_all.tsv"
ANNOTATION_FILE <- "110079946.protein.enrichment.terms.v12.0.txt"

OUTPUT_DIR <- "05_gsva"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Load filtered proteome matrix
# Expected orientation in the RDS: samples × proteins
proteome_samples_x_proteins <- readRDS(PROTEOME_RDS)

if (!is.matrix(proteome_samples_x_proteins)) {
  proteome_samples_x_proteins <- as.matrix(proteome_samples_x_proteins)
}

# GSVA requires proteins × samples
proteome <- t(proteome_samples_x_proteins)

storage.mode(proteome) <- "numeric"

cat(
  "GSVA input matrix:",
  nrow(proteome), "proteins ×",
  ncol(proteome), "samples\n"
)

# Keep only the first protein ID before ";"
rownames(proteome) <- sub(";.*$", "", rownames(proteome))

if (anyDuplicated(rownames(proteome))) {
  duplicated_ids <- unique(rownames(proteome)[duplicated(rownames(proteome))])
  
  stop(
    "Duplicated protein IDs after removing text after ';': ",
    paste(head(duplicated_ids, 10), collapse = ", "),
    ifelse(length(duplicated_ids) > 10, " ...", "")
  )
}

# Load metadata
metadata <- read_tsv(METADATA_FILE, show_col_types = FALSE) %>%
  as.data.frame()

rownames(metadata) <- metadata[[1]]
metadata <- metadata[, -1, drop = FALSE]

if (!all(colnames(proteome) %in% rownames(metadata))) {
  unmatched_samples <- setdiff(colnames(proteome), rownames(metadata))
  
  stop(
    "Proteome samples not found in metadata: ",
    paste(head(unmatched_samples, 10), collapse = ", "),
    ifelse(length(unmatched_samples) > 10, " ...", "")
  )
}

metadata <- metadata[colnames(proteome), , drop = FALSE]

stopifnot(identical(colnames(proteome), rownames(metadata)))

# Load STRING pathway annotations
ann <- read_tsv(
  ANNOTATION_FILE,
  comment = "#",
  col_names = c(
    "string_protein_id",
    "category",
    "term",
    "description"
  ),
  show_col_types = FALSE
) %>%
  mutate(
    protein_id = sub("^[^.]+\\.", "", string_protein_id)
  )

ann_bp <- ann %>%
  filter(category == "Biological Process (Gene Ontology)")

ann_mf <- ann %>%
  filter(category == "Molecular Function (Gene Ontology)")

ann_cc <- ann %>%
  filter(category == "Cellular Component (Gene Ontology)")

ann_kegg <- ann %>%
  filter(str_detect(category, "KEGG"))

ann_rp <- ann %>%
  filter(category == "Reactome Pathways")

t2g_bp <- ann_bp %>%
  transmute(term = term, gene = protein_id) %>%
  distinct()

t2n_bp <- ann_bp %>%
  dplyr::select(term, name = description) %>%
  distinct()

t2g_mf <- ann_mf %>%
  transmute(term = term, gene = protein_id) %>%
  distinct()

t2n_mf <- ann_mf %>%
  dplyr::select(term, name = description) %>%
  distinct()

t2g_cc <- ann_cc %>%
  transmute(term = term, gene = protein_id) %>%
  distinct()

t2n_cc <- ann_cc %>%
  dplyr::select(term, name = description) %>%
  distinct()

t2g_kegg <- ann_kegg %>%
  transmute(term = term, gene = protein_id) %>%
  distinct()

t2n_kegg <- ann_kegg %>%
  dplyr::select(term, name = description) %>%
  distinct()

t2g_rp <- ann_rp %>%
  transmute(term = term, gene = protein_id) %>%
  distinct()

t2n_rp <- ann_rp %>%
  dplyr::select(term, name = description) %>%
  distinct()

keep_ids <- unique(c(
  t2g_bp$gene,
  t2g_mf$gene,
  t2g_cc$gene,
  t2g_kegg$gene,
  t2g_rp$gene
))

# Retain proteins present in at least one annotation category
common_ids <- intersect(rownames(proteome), keep_ids)

if (length(common_ids) == 0) {
  stop("No protein IDs overlap between the proteome matrix and STRING annotations.")
}

prot_gsva <- proteome[common_ids, , drop = FALSE]

cat(
  "Annotated proteins retained for GSVA:",
  nrow(prot_gsva), "\n"
)

# Create gene-set lists
make_gs_list <- function(t2g, genes_in_matrix, min.size = 5, max.size = 500) {
  gene_sets <- t2g %>%
    filter(gene %in% genes_in_matrix) %>%
    group_by(term) %>%
    summarise(genes = list(unique(gene)), .groups = "drop") %>%
    deframe()
  
  gene_sets[lengths(gene_sets) >= min.size & lengths(gene_sets) <= max.size]
}

gs_bp <- make_gs_list(t2g_bp, rownames(prot_gsva))
gs_mf <- make_gs_list(t2g_mf, rownames(prot_gsva))
gs_cc <- make_gs_list(t2g_cc, rownames(prot_gsva))
gs_kegg <- make_gs_list(t2g_kegg, rownames(prot_gsva))
gs_rp <- make_gs_list(t2g_rp, rownames(prot_gsva))

names(gs_bp) <- paste0("BP:", names(gs_bp))
names(gs_mf) <- paste0("MF:", names(gs_mf))
names(gs_cc) <- paste0("CC:", names(gs_cc))
names(gs_kegg) <- paste0("KEGG:", names(gs_kegg))
names(gs_rp) <- paste0("RP:", names(gs_rp))

gs_all <- c(gs_bp, gs_mf, gs_cc, gs_kegg, gs_rp)

cat("Gene sets used for GSVA:", length(gs_all), "\n")

# Run GSVA
gpar <- gsvaParam(
  exprData = as.matrix(prot_gsva),
  geneSets = gs_all,
  minSize = 5,
  maxSize = 500,
  use = "na.rm",
  sparse = FALSE,
  checkNA = "yes",
  kcdf = "Gaussian",
  maxDiff = TRUE
)

es_all <- gsva(gpar, verbose = TRUE)
# Convert GSVA result to a plain matrix
es_mat <- es_all
class(es_mat) <- "matrix"

cat(
  "GSVA score matrix:",
  nrow(es_mat), "pathways ×",
  ncol(es_mat), "samples\n"
)

# Create a regular data frame before converting wide to long
es_all_df <- as.data.frame(
  es_mat,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

es_all_df$term_full <- rownames(es_mat)

es_all_df <- es_all_df %>%
  dplyr::select(term_full, dplyr::everything())

t2n_all <- bind_rows(
  t2n_bp %>% mutate(term_full = paste0("BP:", term)) %>% select(term_full, name),
  t2n_mf %>% mutate(term_full = paste0("MF:", term)) %>% select(term_full, name),
  t2n_cc %>% mutate(term_full = paste0("CC:", term)) %>% select(term_full, name),
  t2n_kegg %>% mutate(term_full = paste0("KEGG:", term)) %>% select(term_full, name),
  t2n_rp %>% mutate(term_full = paste0("RP:", term)) %>% select(term_full, name)
) %>%
  distinct(term_full, .keep_all = TRUE)
# Convert GSVA scores to long format and add pathway annotation
es_all_long <- tidyr::pivot_longer(
  es_all_df,
  cols = -term_full,
  names_to = "sample",
  values_to = "gsva_score"
) %>%
  left_join(t2n_all, by = "term_full") %>%
  mutate(
    ontology = case_when(
      str_detect(term_full, "^BP:") ~ "BP",
      str_detect(term_full, "^MF:") ~ "MF",
      str_detect(term_full, "^CC:") ~ "CC",
      str_detect(term_full, "^KEGG:") ~ "KEGG",
      str_detect(term_full, "^RP:") ~ "Reactome",
      TRUE ~ "Other"
    )
  )
cols_needed <- c(
  "diet",
  "day_cont",
  "day_categ",
  "rearing_tank",
  "group",
  "phase",
  "diet_phase"
)

missing_metadata_columns <- setdiff(cols_needed, colnames(metadata))

if (length(missing_metadata_columns) > 0) {
  stop(
    "Missing metadata columns: ",
    paste(missing_metadata_columns, collapse = ", ")
  )
}

sample_meta_clean <- metadata %>%
  rownames_to_column("sample") %>%
  dplyr::select(sample, all_of(cols_needed))

es_all_long <- es_all_long %>%
  left_join(sample_meta_clean, by = "sample")



# Split GSVA matrices by annotation category
term_full <- rownames(es_mat)

ontology <- case_when(
  str_detect(term_full, "^BP:") ~ "BP",
  str_detect(term_full, "^MF:") ~ "MF",
  str_detect(term_full, "^CC:") ~ "CC",
  str_detect(term_full, "^KEGG:") ~ "KEGG",
  str_detect(term_full, "^RP:") ~ "Reactome",
  TRUE ~ "Other"
)

es_bp <- es_mat[ontology == "BP", , drop = FALSE]
es_mf <- es_mat[ontology == "MF", , drop = FALSE]
es_cc <- es_mat[ontology == "CC", , drop = FALSE]
es_kegg <- es_mat[ontology == "KEGG", , drop = FALSE]
es_rp <- es_mat[ontology == "Reactome", , drop = FALSE]

# Save GSVA outputs using term IDs
gsva_list <- list(
  es_all = es_mat,
  es_bp = es_bp,
  es_mf = es_mf,
  es_cc = es_cc,
  es_kegg = es_kegg,
  es_rp = es_rp,
  t2n_all = t2n_all,
  samples = colnames(es_mat),
  metadata = sample_meta_clean
)

saveRDS(
  gsva_list,
  file.path(OUTPUT_DIR, "sample_geneterms_matrix.rds")
)

# Replace pathway IDs by descriptive names
rename_rows_to_full_name <- function(mat, t2n_all) {
  term_ids <- rownames(mat)
  
  full_names <- t2n_all$name[
    match(term_ids, t2n_all$term_full)
  ]
  
  full_names[is.na(full_names)] <- term_ids[is.na(full_names)]
  rownames(mat) <- make.unique(full_names)
  
  mat
}

es_all_named <- rename_rows_to_full_name(es_mat, t2n_all)
es_bp_named <- rename_rows_to_full_name(es_bp, t2n_all)
es_mf_named <- rename_rows_to_full_name(es_mf, t2n_all)
es_cc_named <- rename_rows_to_full_name(es_cc, t2n_all)
es_kegg_named <- rename_rows_to_full_name(es_kegg, t2n_all)
es_rp_named <- rename_rows_to_full_name(es_rp, t2n_all)

gsva_list_named <- list(
  es_all = es_all_named,
  es_bp = es_bp_named,
  es_mf = es_mf_named,
  es_cc = es_cc_named,
  es_kegg = es_kegg_named,
  es_rp = es_rp_named,
  t2n_all = t2n_all,
  samples = colnames(es_all_named),
  metadata = sample_meta_clean
)

saveRDS(
  gsva_list_named,
  file.path(OUTPUT_DIR, "sample_geneterms_matrix_fullnames.rds")
)

write_tsv(
  es_all_long,
  file.path(OUTPUT_DIR, "GSVA_scores_long.tsv")
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(OUTPUT_DIR, "sessionInfo.txt")
)

cat("GSVA analysis complete.\n")
cat("Outputs saved in:", OUTPUT_DIR, "\n")
