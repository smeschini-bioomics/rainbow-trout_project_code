#Load libraries
library(microeco) 
library(file2meco)
library(mecodev)
library(magrittr)
library(stringr)

set.seed(123)

## Import data from qiime2 output

# ASVs count table
abund_file_path <- "STPN2309_decontam-filtered-merged-table-freq1-with-phyla.qza"
# metadata
sample_file_path <- "metadata.tsv"
# taxonomy file
taxonomy_file_path <- "STPN2309-taxonomy-silva-138.2.qza"
# phylogenetic tree
tree_data <- "STPN2309-silva_nb-rooted-tree.qza"
# fasta seqs
rep_data <- "STPN2309_decontam-filtered-merged-rep-seqs-freq1-with-phyla.qza"


#Create microtable object
mt <- qiime2meco(abund_file_path,
                 sample_table = sample_file_path,
                 taxonomy_table = taxonomy_file_path,
                 rep_fasta = rep_data,
                 phylo_tree = tree_data,
                 auto_tidy = TRUE
)
mt

# check the sequence numbers in each sample
mt$sample_sums() %>% range

# sum the abundance for each taxon
mt$taxa_sums() %>% range

# show sample name
mt$sample_names()

# show taxa name
mt$taxa_names()

# Remove ASVs which are not assigned in the Kingdom “k__Archaea” or “k__Bacteria”
mt$tax_table %<>% base::subset(Kingdom == "k__Bacteria")
mt$tidy_dataset()
mt

# filter taxa
mt$filter_pollution(taxa = c("Mitochondria"))
mt$tidy_dataset()
mt
mt$filter_pollution(taxa = c("Chloroplast"))
mt$tidy_dataset()
mt

# check the sequence numbers in each sample
mt$sample_sums() %>% range

# sum the abundance for each taxon
mt$taxa_sums() %>% range

# show sample name
mt$sample_names()

# show taxa name
mt$taxa_names()


# export tables
output_dir <- "01_output_files"
dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
mt$save_table(dirpath = "01_16s-basic_files", sep = ",")

# Rename featuresIDs to "ASV", including the rownames of otu_table, rownames of tax_table, tip labels of phylo_tree and rep_fasta
mt$rename_taxa(newname_prefix = "ASV_")

# Add the rownames of microtable$tax_table as its last column
mt$add_rownames2taxonomy(use_name = "ASV")

# order the factors level
mt$sample_table$Diet %<>% factor(., levels = c("Fasted",
                                               "NC",
                                               "HC"))
mt$sample_table$sample_type %<>% factor(., levels = c("Feed",
                                                      "Digesta",
                                                      "Mucus"))

mt$sample_table$sample_type_diet %<>% factor(., levels = c("Digesta_Fasted",
                                                           "Mucus_Fasted",
                                                           "Feed_NC",
                                                           "Feed_HC",
                                                           "Digesta_NC",
                                                           "Digesta_HC",
                                                           "Mucus_NC",
                                                           "Mucus_HC"))


mt$sample_table$sample_type_diet_day %<>% factor(., levels = c("Feed_NC",
                                                               "Feed_HC",
                                                               "Digesta_Fasted",
                                                               "Digesta_NC_day1",
                                                               "Digesta_NC_day2",
                                                               "Digesta_NC_day3",
                                                               "Digesta_NC_day4",
                                                               "Digesta_NC_day10",
                                                               "Digesta_NC_day15",
                                                               "Digesta_NC_day22",
                                                               "Digesta_HC_day1",
                                                               "Digesta_HC_day2",
                                                               "Digesta_HC_day3",
                                                               "Digesta_HC_day4",
                                                               "Digesta_HC_day10",
                                                               "Digesta_HC_day15",
                                                               "Digesta_HC_day22",
                                                               "Mucus_Fasted",
                                                               "Mucus_NC_day1",
                                                               "Mucus_NC_day2",
                                                               "Mucus_NC_day3",
                                                               "Mucus_NC_day4",
                                                               "Mucus_NC_day10",
                                                               "Mucus_NC_day15",
                                                               "Mucus_NC_day22",
                                                               "Mucus_HC_day1",
                                                               "Mucus_HC_day2",
                                                               "Mucus_HC_day3",
                                                               "Mucus_HC_day4",
                                                               "Mucus_HC_day10",
                                                               "Mucus_HC_day15",
                                                               "Mucus_HC_day22"))


# ---- RData output directory ----
rdata_dir <- "RData"

dir.create(
  rdata_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Save dataset
mt_raw <- clone(mt)

save(
  mt_raw,
  file = file.path(rdata_dir, "16s_dataset.RData")
)

# load the data
#load("RData/16s_dataset.RData")

# Dataset Rarefaction
# rarefaction curve: "Observed"
t1 <- trans_rarefy$new(dataset = mt_raw, method = "SRS",
                       alphadiv = c("Observed"),
                       depth = c(0, 10000, 20000,
                                 50000, 80000, 100000,
                                 150000, 200000))
t1$plot_rarefy(color = "sample_type_diet_day",
               show_legend = FALSE,
               show_point = FALSE,
               show_samplename = FALSE,
               add_fitting = FALSE)


# rarefaction curve: "Shannon"
t1 <- trans_rarefy$new(dataset = mt_raw,
                       method = "SRS",
                       alphadiv = c("Shannon"),
                       depth = c(0, 10000, 20000,
                                 50000, 80000, 100000,
                                 150000, 200000))
t1$plot_rarefy(color = "sample_type_diet_day",
               show_legend = FALSE,
               show_point = FALSE,
               show_samplename = FALSE,
               add_fitting = FALSE)


# rarefaction curve: "Pielou"
t1 <- trans_rarefy$new(dataset = mt_raw,
                       method = "SRS",
                       alphadiv = c("Pielou"),
                       depth = c(0, 10000, 20000,
                                 50000, 80000, 100000,
                                 150000, 200000))
t1$plot_rarefy(color = "sample_type_diet_day",
               show_legend = FALSE,
               show_point = FALSE,
               show_samplename = FALSE,
               add_fitting = FALSE)


# clone the data for rarefaction
mt_rarefied <- clone(mt_raw)
mt_rarefied
# filter
mt_rarefied$sample_table <- subset(mt_rarefied$sample_table,
                                   Diet %in% c("NC", "HC")) #select your day of interest
mt_rarefied$sample_table <- subset(mt_rarefied$sample_table,
                                   sample_type %in% c("Digesta", "Mucus")) #select your day of interest

mt_rarefied$tidy_dataset()
mt_rarefied

# Convert Day to a factor, then drop unused levels from all factor columns
mt_rarefied$sample_table$Day <- as.factor(mt_rarefied$sample_table$Day)
is_factor <- vapply(mt_rarefied$sample_table, is.factor, logical(1))
mt_rarefied$sample_table[is_factor] <- lapply(
  mt_rarefied$sample_table[is_factor],
  droplevels
)
mt_rarefied$tidy_dataset()
mt_rarefied



# check the sequence numbers in each sample
mt_rarefied$sample_sums() %>% range
# sum the abundance for each taxon
mt_rarefied$taxa_sums() %>% range

# "SRS": scaling with ranked subsampling method based on the SRS package provided by Lukas Beule and Petr Karlovsky (2020) <doi:10.7717/peerj.9593>
mt_rarefied$rarefy_samples(method = "SRS",
                           sample.size = 50000)
# check the sequence numbers in each sample
mt_rarefied$sample_sums() %>% range
# sum the abundance for each taxon
mt_rarefied$taxa_sums() %>% range

# we call alpha diversity
# SRS dataset
mt_rarefied$cal_alphadiv(measures = c("Observed",
                                          "Shannon",
                                          "Pielou"), PD = TRUE)
# save alpha diversity table to the computer
diversity_dir <- "01_diversity_tables"

dir.create(
  diversity_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Save alpha-diversity table
mt_rarefied$save_alphadiv(
  dirpath = file.path(diversity_dir, "alpha_div_rarefied_50k")
)

## we call beta diversity
mt_rarefied$cal_betadiv(method = c("bray", "jaccard"),
                            unifrac = TRUE,
                            binary = FALSE)
# save beta diversity matrix to the computer
mt_rarefied$save_betadiv(
  dirpath = file.path(diversity_dir, "beta_div_rarefied_50k")
)

# order the factors level
mt_rarefied$sample_table$Diet %<>% factor(., levels = c("NC",
                                                        "HC"))
mt_rarefied$sample_table$sample_type %<>% factor(., levels = c("Digesta",
                                                               "Mucus"))

mt_rarefied$sample_table$sample_type_diet_day %<>% factor(., levels = c(
  "Digesta_NC_day1",
  "Digesta_NC_day2",
  "Digesta_NC_day3",
  "Digesta_NC_day4",
  "Digesta_NC_day10",
  "Digesta_NC_day15",
  "Digesta_NC_day22",
  "Digesta_HC_day1",
  "Digesta_HC_day2",
  "Digesta_HC_day3",
  "Digesta_HC_day4",
  "Digesta_HC_day10",
  "Digesta_HC_day15",
  "Digesta_HC_day22",
  "Mucus_NC_day1",
  "Mucus_NC_day2",
  "Mucus_NC_day3",
  "Mucus_NC_day4",
  "Mucus_NC_day10",
  "Mucus_NC_day15",
  "Mucus_NC_day22",
  "Mucus_HC_day1",
  "Mucus_HC_day2",
  "Mucus_HC_day3",
  "Mucus_HC_day4",
  "Mucus_HC_day10",
  "Mucus_HC_day15",
  "Mucus_HC_day22"))

mt_rarefied$sample_table$sample_type_diet %<>% factor(., levels = c(
  "Digesta_NC",
  "Digesta_HC",
  "Mucus_NC",
  "Mucus_HC"))

mt_rarefied$tidy_dataset()
class(mt_rarefied)

save(mt_rarefied, file = "Rdata/16s_dataset_rarefied.RData")



# Subset the dataset to include only "Digesta" samples
mt_dig <- clone(mt_raw)
# check the sequence numbers in each sample
mt_dig$sample_sums() %>% range
# sum the abundance for each taxon
mt_dig$taxa_sums() %>% range

mt_dig$sample_table <- subset(mt_dig$sample_table,
                              sample_type %in% c("Digesta"))
mt_dig$sample_table <- subset(mt_dig$sample_table,
                              Diet %in% c("NC", "HC"))

# trim all the data
mt_dig$tidy_dataset()
# check the sequence numbers in each sample
mt_dig$sample_sums() %>% range
# sum the abundance for each taxon
mt_dig$taxa_sums() %>% range

# save taxonomic abundance as local file
out_dir <- file.path("01_ASVs_relative_abundance_table", "digesta")

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
mt_dig$cal_abund()
mt_dig$save_abund(dirpath = out_dir)

# Convert Day to a factor, then drop unused levels from all factor columns
mt_dig$sample_table$Day <- as.factor(mt_dig$sample_table$Day)
is_factor <- vapply(mt_dig$sample_table, is.factor, logical(1))
mt_dig$sample_table[is_factor] <- lapply(
  mt_dig$sample_table[is_factor],
  droplevels
)
mt_dig$tidy_dataset()

# we order factor levels
mt_dig$sample_table$Diet %<>% factor(., levels = c("NC", "HC"))


save(mt_dig, file = "RData/16s_dataset_digesta.RData")



# digesta dataset rarefied
mt_dig_rarefied <- clone(mt_raw)
mt_dig_rarefied$taxa_sums() %>% range
mt_dig_rarefied$sample_table <- subset(mt_dig_rarefied$sample_table,
                                       sample_type %in% c("Digesta"))
mt_dig_rarefied$sample_table <- subset(mt_dig_rarefied$sample_table,
                                       Diet %in% c("NC", "HC"))
mt_dig_rarefied$tidy_dataset()
# check the sequence numbers in each sample
mt_dig_rarefied$sample_sums() %>% range
# sum the abundance for each taxon
mt_dig_rarefied$taxa_sums() %>% range

# Convert Day to a factor, then drop unused levels from all factor columns
mt_dig_rarefied$sample_table$Day <- as.factor(mt_dig_rarefied$sample_table$Day)
is_factor <- vapply(mt_dig_rarefied$sample_table, is.factor, logical(1))
mt_dig_rarefied$sample_table[is_factor] <- lapply(
  mt_dig_rarefied$sample_table[is_factor],
  droplevels
)
mt_dig_rarefied$tidy_dataset()

# "SRS": scaling with ranked subsampling method based on the SRS package provided by Lukas Beule and Petr Karlovsky (2020) <doi:10.7717/peerj.9593>
mt_dig_rarefied$rarefy_samples(method = "SRS",
                               sample.size = 50000)
mt_dig_rarefied$tidy_dataset()
# check the sequence numbers in each sample
mt_dig_rarefied$sample_sums() %>% range
# sum the abundance for each taxon
mt_dig_rarefied$taxa_sums() %>% range


# we call alpha diversity
# SRS dataset
mt_dig_rarefied$cal_alphadiv(measures = c("Observed",
                                          "Shannon",
                                          "Pielou"), PD = TRUE)
# save alpha diversity table to the computer
diversity_dir <- "01_diversity_tables"

dir.create(
  diversity_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Save alpha-diversity table
mt_dig_rarefied$save_alphadiv(
  dirpath = file.path(diversity_dir, "alpha_div_digesta_rarefied_50k")
)

## we call beta diversity
mt_dig_rarefied$cal_betadiv(method = c("bray", "jaccard"),
                            unifrac = TRUE,
                            binary = FALSE)
# save beta diversity matrix to the computer
mt_dig_rarefied$save_betadiv(
  dirpath = file.path(diversity_dir, "beta_div_digesta_rarefied_50k")
)

# we order factor levels
mt_dig_rarefied$sample_table$Diet %<>% factor(., levels = c("NC", "HC"))
mt_dig_rarefied$sample_table$Day %<>% factor(., levels = c("day_01",
                                                           "day_02",
                                                           "day_03",
                                                           "day_04",
                                                           "day_10",
                                                           "day_15",
                                                           "day_22"))
mt_dig_rarefied$sample_table$sample_type_diet_day %<>% factor(., levels = c("Digesta_NC_day1", "Digesta_HC_day1",
                                                                            "Digesta_NC_day2", "Digesta_HC_day2",
                                                                            "Digesta_NC_day3", "Digesta_HC_day3",
                                                                            "Digesta_NC_day4", "Digesta_HC_day4",
                                                                            "Digesta_NC_day10", "Digesta_HC_day10",
                                                                            "Digesta_NC_day15", "Digesta_HC_day15",
                                                                            "Digesta_NC_day22", "Digesta_HC_day22"))

save(mt_dig_rarefied, file = "RData/16s_dataset_dig_rarefied.RData")


# Mucus
mt_mucus <- clone(mt_raw)

# check the sequence numbers in each sample
mt_mucus$sample_sums() %>% range
# sum the abundance for each taxon
mt_mucus$taxa_sums() %>% range

mt_mucus$sample_table <- subset(mt_mucus$sample_table,
                                sample_type %in% c("Mucus"))
mt_mucus$sample_table <- subset(mt_mucus$sample_table,
                                Diet %in% c("NC", "HC"))

# trim all the data
mt_mucus$tidy_dataset()
# check the sequence numbers in each sample
mt_mucus$sample_sums() %>% range
# sum the abundance for each taxon
mt_mucus$taxa_sums() %>% range

# save taxonomic abundance as local file
out_dir <- file.path("01_ASVs_relative_abundance_table", "mucus")

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
mt_mucus$cal_abund()
mt_mucus$save_abund(dirpath = out_dir)

# Convert Day to a factor, then drop unused levels from all factor columns
mt_mucus$sample_table$Day <- as.factor(mt_mucus$sample_table$Day)
is_factor <- vapply(mt_mucus$sample_table, is.factor, logical(1))
mt_mucus$sample_table[is_factor] <- lapply(
  mt_mucus$sample_table[is_factor],
  droplevels
)
mt_mucus$tidy_dataset()

# we order factor levels
mt_mucus$sample_table$Diet %<>% factor(., levels = c("NC", "HC"))
mt_mucus$sample_table$sample_type %<>% factor(., levels = c("Mucus"))
mt_mucus$sample_table$Day %<>% factor(., levels = c("day_01",
                                                    "day_02",
                                                    "day_03",
                                                    "day_04",
                                                    "day_10",
                                                    "day_15",
                                                    "day_22"))
mt_mucus$sample_table$sample_type_diet_day %<>% factor(., levels = c(
  "Mucus_NC_day1", "Mucus_NC_day2",
  "Mucus_NC_day3", "Mucus_NC_day4",
  "Mucus_NC_day10", "Mucus_NC_day15",
  "Mucus_NC_day22",
  "Mucus_HC_day1", "Mucus_HC_day2",
  "Mucus_HC_day3", "Mucus_HC_day4",
  "Mucus_HC_day10",
  "Mucus_HC_day15", "Mucus_HC_day22"))

mt_mucus

save(mt_mucus, file = "RData/16s_dataset_mucus.RData")


# mucus rarefied
mt_mucus_rarefied <- clone(mt_raw)
# check the sequence numbers in each sample
mt_mucus_rarefied$sample_sums() %>% range
# sum the abundance for each taxon
mt_mucus_rarefied$taxa_sums() %>% range
mt_mucus_rarefied$sample_table <- subset(mt_mucus_rarefied$sample_table,
                                         sample_type %in% c("Mucus"))
mt_mucus_rarefied$sample_table <- subset(mt_mucus_rarefied$sample_table,
                                         Diet %in% c("NC", "HC"))
mt_mucus_rarefied$tidy_dataset()
# check the sequence numbers in each sample
mt_mucus_rarefied$sample_sums() %>% range
# sum the abundance for each taxon
mt_mucus_rarefied$taxa_sums() %>% range

# "SRS": scaling with ranked subsampling method based on the SRS package provided by Lukas Beule and Petr Karlovsky (2020) <doi:10.7717/peerj.9593>
mt_mucus_rarefied$rarefy_samples(method = "SRS",
                                 sample.size = 50000)
mt_mucus_rarefied$tidy_dataset()
# check the sequence numbers in each sample
mt_mucus_rarefied$sample_sums() %>% range
# sum the abundance for each taxon
mt_mucus_rarefied$taxa_sums() %>% range

# Convert Day to a factor, then drop unused levels from all factor columns
mt_mucus_rarefied$sample_table$Day <- as.factor(mt_mucus_rarefied$sample_table$Day)
is_factor <- vapply(mt_mucus_rarefied$sample_table, is.factor, logical(1))
mt_mucus_rarefied$sample_table[is_factor] <- lapply(
  mt_mucus_rarefied$sample_table[is_factor],
  droplevels
)
mt_mucus_rarefied$tidy_dataset()

# we call alpha diversity
# SRS dataset
mt_mucus_rarefied$cal_alphadiv(measures = c("Observed",
                                          "Shannon",
                                          "Pielou"), PD = TRUE)
# save alpha diversity table to the computer
diversity_dir <- "01_diversity_tables"

dir.create(
  diversity_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Save alpha-diversity table
mt_mucus_rarefied$save_alphadiv(
  dirpath = file.path(diversity_dir, "alpha_div_mucus_rarefied_50k")
)

## we call beta diversity
mt_mucus_rarefied$cal_betadiv(method = c("bray", "jaccard"),
                            unifrac = TRUE,
                            binary = FALSE)
# save beta diversity matrix to the computer
mt_mucus_rarefied$save_betadiv(
  dirpath = file.path(diversity_dir, "beta_div_mucus_rarefied_50k")
)

# we order factor levels
mt_mucus_rarefied$sample_table$Diet %<>% factor(., levels = c("NC", "HC"))
mt_mucus_rarefied$sample_table$Day %<>% factor(., levels = c("day_01",
                                                             "day_02",
                                                             "day_03",
                                                             "day_04",
                                                             "day_10",
                                                             "day_15",
                                                             "day_22"))
mt_mucus_rarefied$sample_table$sample_type_diet_day %<>% factor(., levels = c("Mucus_NC_day1", "Mucus_HC_day1",
                                                                              "Mucus_NC_day2", "Mucus_HC_day2",
                                                                              "Mucus_NC_day3", "Mucus_HC_day3",
                                                                              "Mucus_NC_day4", "Mucus_HC_day4",
                                                                              "Mucus_NC_day10", "Mucus_HC_day10",
                                                                              "Mucus_NC_day15", "Mucus_HC_day15",
                                                                              "Mucus_NC_day22", "Mucus_HC_day22"))
save(mt_mucus_rarefied, file = "RData/16s_dataset_mucus_rarefied.RData")





# load dataset and rename factor for MicroViz
load("RData/16s_dataset.RData")
mt_abund <- clone(mt_raw)
mt_abund$sample_table <- subset(mt_abund$sample_table,
                                sample_type %in% c("Digesta",
                                                   "Mucus"))
mt_abund$sample_table <- subset(mt_abund$sample_table,
                                Diet %in% c("NC", "HC"))
mt_abund$tidy_dataset()

# Convert Day to a factor, then drop unused levels from all factor columns
mt_abund$sample_table$Day <- as.factor(mt_abund$sample_table$Day)
is_factor <- vapply(mt_abund$sample_table, is.factor, logical(1))
mt_abund$sample_table[is_factor] <- lapply(
  mt_abund$sample_table[is_factor],
  droplevels
)
mt_abund$tidy_dataset()
mt_abund

mt_abund$cal_abund()
mt_abund

# save taxonomic abundance as local file
mt_abund$save_abund(dirpath = "01_ASVs_relative_abundance_table")

# rename "Day"
mt_abund$sample_table$Day_clean <- sapply(
  mt_abund$sample_table$Day,
  function(v) {
    num <- as.integer(str_match(v, "day_(\\d+)")[,2])
    paste("Day", num)
  }
)
# order "Day"
mt_abund$sample_table$Day_clean %<>% factor(., levels = c("Day 1","Day 2",
                                                          "Day 3", "Day 4",
                                                          "Day 10", "Day 15",
                                                          "Day 22"))

# rename "sample_type_diet": replace "_" with space
mt_abund$sample_table$sample_type_diet_clean <- gsub(
  "_", 
  " ", 
  mt_abund$sample_table$sample_type_diet
)

# order "sample_type_diet"
mt_abund$sample_table$sample_type_diet_clean %<>% factor(
  .,
  levels = c(
    "Digesta NC",
    "Digesta HC",
    "Mucus NC",
    "Mucus HC"
  )
)
save(mt_abund, file = "RData/16s_dataset_microViz.RData")
