# libraries
library(phyloseq)
library(microViz)
library(corncob)
library(ape)
library(dplyr)
library(stringr)
library(patchwork)
library(ggplot2)
library(ggrepel)

set.seed(123)

# Set Output directory
out_dir <- "07_DA_analysis_mucus_microViz"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Import filtered BIOM table used for DA analysis in BIRDMAn
ps <- suppressWarnings(
  import_biom("DA_analysis_BIRDMAn/mucus/output_files/mucus_filtered_table_with_sample_metadata.biom")
)


# Check imported object
ps
sample_variables(ps)
rank_names(ps)
sample_data(ps)[1:3, 1:5]
tax_table(ps)[1:3, ]


# Rename taxonomy ranks
colnames(tax_table(ps)) <- c(
  "Kingdom", "Phylum", "Class", "Order",
  "Family", "Genus", "Species"
)

rank_names(ps)
tax_table(ps)[1:5, ]

# Remove Species rank
# Species is mostly empty and we are modelling only Phylum to Genus.
tax <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)
tax$Species <- NULL
tax_table(ps) <- as.matrix(tax)

rank_names(ps)
tax_table(ps)[1:5, ]

# Check sample metadata
sample_data(ps)$Diet <- factor(
  sample_data(ps)$Diet,
  levels = c("NC", "HC")
)

sample_data(ps)$Day <- factor(
  sample_data(ps)$Day,
  levels = c(
    "day_01", "day_02", "day_03", "day_04",
    "day_10", "day_15", "day_22"
  )
)

# Check metadata
table(sample_data(ps)$Diet)
table(sample_data(ps)$Day)
sample_data(ps)[1:5, c("Diet", "Day")]

# Optional: import and match phylogenetic tree
# taxatree_plots() uses taxonomy, not phylogenetic branch lengths,
# but keeping the tree in the phyloseq object is useful.
tree <- ape::read.tree("05_export_files_for_BIRDMAn_DA_analysis/tree.nwk")

cat("Taxa in BIOM:", ntaxa(ps), "\n")
cat("Tips in tree:", length(tree$tip.label), "\n")
cat("Shared taxa:", length(intersect(taxa_names(ps), tree$tip.label)), "\n")

common_taxa <- intersect(taxa_names(ps), tree$tip.label)

ps <- prune_taxa(common_taxa, ps)
tree <- ape::keep.tip(tree, common_taxa)

ps <- merge_phyloseq(ps, tree)

stopifnot(setequal(taxa_names(ps), phy_tree(ps)$tip.label))

ps

# Prepare object exactly in tutorial style
ps_bb <- ps %>%
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(
    min_prevalence = 0.10,
    undetected = 0,
    use_counts = TRUE
  )

# Check after fixing/prepending/filtering
ps_bb
rank_names(ps_bb)
tax_table(ps_bb)[1:10, ]
table(sample_data(ps_bb)$Diet)

cat("Taxa after filtering:", ntaxa(ps_bb), "\n")
cat("Samples after filtering:", nsamples(ps_bb), "\n")

# Save prepared object
saveRDS(
  ps_bb,
  file.path(out_dir, "ps_microViz_prepared_taxfix_prepend_filtered.rds")
)

# Test model on a few taxa first
bb_test <- ps_bb %>%
  tax_model(
    type = corncob::bbdml,
    rank = "Genus",
    taxa = 1:3,
    variables = "Diet",
    return_psx = FALSE
  )

bb_test
summary(bb_test[[1]])


# Run beta-binomial models across ranks
bb_models <- ps_bb %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
    variables = "Diet"
  )

bb_models

saveRDS(
  bb_models,
  file.path(out_dir, "beta_binomial_taxatree_models_Diet.rds")
)

# Extract abundance statistics: mu
bb_stats <- taxatree_models2stats(
  bb_models,
  param = "mu"
)

bb_stats

bb_df <- bb_stats %>%
  taxatree_stats_get()

head(bb_df)
unique(bb_df$term)
unique(bb_df$rank)


# Adjust p-values by rank, as in tutorial
bb_stats <- taxatree_stats_p_adjust(
  data = bb_stats,
  method = "BH",
  grouping = "rank"
)

bb_results <- bb_stats %>%
  taxatree_stats_get() %>%
  filter(term == "DietHC") %>%
  arrange(p.adj.BH.rank)

bb_results %>%
  select(
    term, taxon, rank, estimate, std.error,
    t.statistic, p.value, p.adj.BH.rank
  ) %>%
  head(30)

# Save results
write.csv(
  bb_results,
  file.path(out_dir, "beta_binomial_taxatree_DietHC_vs_NC_results.csv"),
  row.names = FALSE
)

saveRDS(
  bb_stats,
  file.path(out_dir, "beta_binomial_taxatree_stats_Diet_BH.rds")
)

# Final checkpoint
cat("\nFinal checkpoint\n")
cat("Samples:", nsamples(ps_bb), "\n")
cat("Taxa in prepared object:", ntaxa(ps_bb), "\n")
cat("Ranks modelled:", paste(unique(bb_results$rank), collapse = ", "), "\n")
cat("Terms:", paste(unique(bb_results$term), collapse = ", "), "\n")
cat("Significant taxa BH < 0.05:", sum(bb_results$p.adj.BH.rank < 0.05, na.rm = TRUE), "\n")
cat("Output directory:", out_dir, "\n")


# Taxatree plot showing both adjusted and raw significance

tree_labelled_sig <- bb_stats %>%
  taxatree_label(
    p.adj.BH.rank < 0.05
    #rank == "Genus"
  ) %>%
  taxatree_plots(
    sig_stat = c("p.adj.BH.rank", "p.value"),
    sig_threshold = 0.05,
    sig_shape = c("cross", "circle filled"),
    sig_colour = "white",
    sig_size = c(1.5, 1),
    sig_stroke = c(1, 0.25),
    title_size = 14,
    drop_ranks = TRUE,
    node_size_range = c(1, 5),
    colour_trans = "identity",
    colour_lims = c(-2, 2),
    colour_oob = scales::oob_squish,
    var_renamer = function(x) "Mucus-Associated Microbiota HC vs NC"
  ) %>%
  .[[1]] %>%
  taxatree_plot_labels(
    taxon_renamer = function(x) {
      x %>%
        stringr::str_remove("^[KPCOFG]:\\s*")
    },
    fun = ggrepel::geom_label_repel,
    x_nudge = 0.15,
    hjust = 0.5,
    size = 2.5
  )

tree_labelled_sig
tree_labelled_sig <- tree_labelled_sig +
  theme(
    plot.title = element_text(family = "mono")
  ) + theme(
  legend.position = c(1.25, 0.5),
  legend.justification = c(0, 0.5)
)

tree_labelled_sig
ggsave(
  file.path(out_dir, "beta_binomial_taxatree_HC_vs_NC_labelled_all_BH_significant.pdf"),
  tree_labelled_sig,
  width = 9,
  height = 4.5,
  dpi = 1200
)

ggsave(
  file.path(out_dir, "beta_binomial_taxatree_HC_vs_NC_labelled_all_BH_significant.tiff"),
  tree_labelled_sig,
  width = 9,
  height = 4.5,
  dpi = 1200,
  compression = "lzw"
)

# Record package and R versions
writeLines(
  capture.output(sessionInfo()),
  con = file.path(out_dir, "sessionInfo.txt")
)

