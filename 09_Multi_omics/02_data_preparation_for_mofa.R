# libraries
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(MOFA2)
  library(matrixStats)
  library(ggplot2)
})

# Load data
liver_hPTM <- readRDS("input_files/H3_H4_8PTM_liver_batch_corrected_limma.rds")
liver_hPTM <- liver_hPTM$ptm_batch_corrected_limma
midgut_hPTM <- readRDS("input_files/H3_H4_5PTM_midgut_batch_corrected_limma.rds") 
midgut_hPTM <- midgut_hPTM$ptm_batch_corrected_limma
digesta_microbiota <- readRDS("input_files/digesta_ASVbest_sample_by_feature.rds")
mucus_microbiota <- readRDS("input_files/mucus_ASVbest_sample_by_feature.rds")
midgut_proteome <- readRDS("input_files/midgut_proteome_filtered_batch_corrected_limma.rds")
midgut_proteome <- midgut_proteome$proteome_batch_corrected_limma
liver_proteome <- readRDS("input_files/liver_proteome_filtered_batch_corrected_limma.rds")
liver_proteome <- liver_proteome$proteome_batch_corrected_limma

data <- readRDS("metadata_for_omics_data.rds")
names(data)
metadata <- data$metadata
host_traits_all <- data$host_traits

methylome <- as.data.frame(data$Liver_DNA_methylation)
methylome_clean <- methylome[, c("dC_rel", "mdC_rel", "hmdC_rel")] %>%
  rename(
    C = dC_rel,
    `5-mC`  = mdC_rel,
    `5-hmC` = hmdC_rel
  )
stopifnot(!is.null(rownames(methylome_clean)))

# Keep only core experimental metadata
meta_core_cols <- c(
  "diet",
  "day_categ",
  "phase",
  "diet_phase"
)

metadata_core <- metadata[, meta_core_cols, drop = FALSE]

stopifnot(identical(rownames(metadata_core), rownames(metadata)))

# Remove fasted condition
meta_filtered <- metadata_core[metadata_core$diet != "Fasted", , drop = FALSE]

table(meta_filtered$diet)
dim(meta_filtered)

# Align host traits to non-fasted metadata
stopifnot(!is.null(rownames(meta_filtered)))
stopifnot(!is.null(rownames(host_traits_all)))

host_traits_for_meta <- host_traits_all[rownames(meta_filtered), , drop = FALSE]
stopifnot(identical(rownames(meta_filtered), rownames(host_traits_for_meta)))

# Bind design metadata + host traits
meta_filtered_with_traits <- cbind(meta_filtered, host_traits_for_meta)

# Add methylome variables to metadata
meta_traits_tbl <- meta_filtered_with_traits %>%
  rownames_to_column("sample_id")

methylome_tbl <- methylome_clean %>%
  rownames_to_column("sample_id")

meta_full_tbl <- meta_traits_tbl %>%
  left_join(methylome_tbl, by = "sample_id")

meta_filtered_with_traits <- meta_full_tbl %>%
  column_to_rownames("sample_id")

dim(meta_filtered_with_traits)
colnames(meta_filtered_with_traits)


# Convert all omics to features × samples

HPTM_mat              <- t(as.matrix(midgut_hPTM))
Liver_HPTM_mat        <- t(as.matrix(liver_hPTM))
midgut_proteome_mat   <- t(as.matrix(midgut_proteome))
liver_proteome_mat    <- t(as.matrix(liver_proteome))
MicrobiotaDigesta_mat <- t(as.matrix(digesta_microbiota))
MicrobiotaMucus_mat   <- t(as.matrix(mucus_microbiota))
Host_Traits_mat       <- t(as.matrix(host_traits_all))

class(HPTM_mat)              <- "matrix"
class(Liver_HPTM_mat)        <- "matrix"
class(midgut_proteome_mat)   <- "matrix"
class(liver_proteome_mat)    <- "matrix"
class(MicrobiotaDigesta_mat) <- "matrix"
class(MicrobiotaMucus_mat)   <- "matrix"
class(Host_Traits_mat)       <- "matrix"


# Keep only non-fasted samples available in each view
get_nonfasted_view <- function(M, meta_filtered_with_traits) {
  common_samples <- intersect(colnames(M), rownames(meta_filtered_with_traits))
  M[, common_samples, drop = FALSE]
}

HPTM_common              <- get_nonfasted_view(HPTM_mat, meta_filtered_with_traits)
Liver_HPTM_common        <- get_nonfasted_view(Liver_HPTM_mat, meta_filtered_with_traits)
midgut_proteome_common   <- get_nonfasted_view(midgut_proteome_mat, meta_filtered_with_traits)
liver_proteome_common    <- get_nonfasted_view(liver_proteome_mat, meta_filtered_with_traits)
MicrobiotaDigesta_common <- get_nonfasted_view(MicrobiotaDigesta_mat, meta_filtered_with_traits)
MicrobiotaMucus_common   <- get_nonfasted_view(MicrobiotaMucus_mat, meta_filtered_with_traits)
Host_Traits_common       <- get_nonfasted_view(Host_Traits_mat, meta_filtered_with_traits)

## Sanity check: number of samples per view
sapply(
  list(
    HPTM = HPTM_common,
    Liver_HPTM = Liver_HPTM_common,
    Midgut_Proteome = midgut_proteome_common,
    Liver_Proteome = liver_proteome_common,
    Microbiota_Digesta = MicrobiotaDigesta_common,
    Microbiota_Mucus = MicrobiotaMucus_common,
    Host_Traits = Host_Traits_common
  ),
  ncol
)

# Feature variance / SD diagnostics
vars <- list(
  HPTM               = rowVars(HPTM_common, na.rm = TRUE),
  Liver_HPTM         = rowVars(Liver_HPTM_common, na.rm = TRUE),
  Midgut_Proteome    = rowVars(midgut_proteome_common, na.rm = TRUE),
  Liver_Proteome     = rowVars(liver_proteome_common, na.rm = TRUE),
  Microbiota_Digesta = rowVars(MicrobiotaDigesta_common, na.rm = TRUE),
  Microbiota_Mucus   = rowVars(MicrobiotaMucus_common, na.rm = TRUE)
)

lapply(vars, summary)

hist(vars$HPTM, breaks = 50, main = "Midgut hPTM feature variances")
hist(vars$Liver_HPTM, breaks = 50, main = "Liver hPTM feature variances")
hist(vars$Midgut_Proteome, breaks = 50, main = "Midgut proteome feature variances")
hist(vars$Liver_Proteome, breaks = 50, main = "Liver proteome feature variances")
hist(vars$Microbiota_Digesta, breaks = 50, main = "Microbiota DIGESTA feature variances")
hist(vars$Microbiota_Mucus, breaks = 50, main = "Microbiota MUCUS feature variances")

sd_ptm        <- sqrt(rowVars(HPTM_common, na.rm = TRUE))
sd_liver_ptm  <- sqrt(rowVars(Liver_HPTM_common, na.rm = TRUE))
sd_prot       <- sqrt(rowVars(midgut_proteome_common, na.rm = TRUE))
sd_liver_prot <- sqrt(rowVars(liver_proteome_common, na.rm = TRUE))
sd_micro_d    <- sqrt(rowVars(MicrobiotaDigesta_common, na.rm = TRUE))
sd_micro_m    <- sqrt(rowVars(MicrobiotaMucus_common, na.rm = TRUE))

summary(sd_ptm)
summary(sd_liver_ptm)
summary(sd_prot)
summary(sd_liver_prot)
summary(sd_micro_d)
summary(sd_micro_m)

# Microbiota prevalence diagnostics
prev_digesta <- rowMeans(MicrobiotaDigesta_common > 0, na.rm = TRUE)
summary(prev_digesta)

sum(prev_digesta == 0)
sum(prev_digesta > 0 & prev_digesta < 0.05)
sum(prev_digesta >= 0.05 & prev_digesta < 0.20)
sum(prev_digesta >= 0.20 & prev_digesta < 0.80)
sum(prev_digesta >= 0.50)
sum(prev_digesta >= 0.80)

keep_ids_micro_digesta <- which(prev_digesta >= 0.00)

prev_mucus <- rowMeans(MicrobiotaMucus_common > 0, na.rm = TRUE)
summary(prev_mucus)

sum(prev_mucus == 0)
sum(prev_mucus > 0 & prev_mucus < 0.05)
sum(prev_mucus >= 0.05 & prev_mucus < 0.20)
sum(prev_mucus >= 0.20 & prev_mucus < 0.80)
sum(prev_mucus >= 0.50)
sum(prev_mucus >= 0.80)

keep_ids_micro_mucus <- which(prev_mucus >= 0.00)


# Proteome SD filtering
cutoffs <- seq(0, 1, by = 0.05)

n_kept_midgut <- sapply(cutoffs, function(th) sum(sd_prot > th, na.rm = TRUE))
n_kept_liver  <- sapply(cutoffs, function(th) sum(sd_liver_prot > th, na.rm = TRUE))

data.frame(
  cutoff = cutoffs,
  midgut_proteins_kept = n_kept_midgut,
  midgut_percent_kept = round(100 * n_kept_midgut / length(sd_prot), 1),
  liver_proteins_kept = n_kept_liver,
  liver_percent_kept = round(100 * n_kept_liver / length(sd_liver_prot), 1)
)

keep_prot <- sd_prot >= 0.30 & !is.na(sd_prot)
keep_liver_prot <- sd_liver_prot >= 0.30 & !is.na(sd_liver_prot)

table(keep_prot)
sum(keep_prot)

table(keep_liver_prot)
sum(keep_liver_prot)

# Create filtered views

HPTM_filtered              <- HPTM_common
Liver_HPTM_filtered        <- Liver_HPTM_common
midgut_proteome_filtered   <- midgut_proteome_common[keep_prot, , drop = FALSE]
liver_proteome_filtered    <- liver_proteome_common[keep_liver_prot, , drop = FALSE]
MicrobiotaDigesta_filtered <- MicrobiotaDigesta_common[keep_ids_micro_digesta, , drop = FALSE]
MicrobiotaMucus_filtered   <- MicrobiotaMucus_common[keep_ids_micro_mucus, , drop = FALSE]
Host_Traits_filtered       <- Host_Traits_common

# Z-score filtered views before alignment
zscore_rows <- function(M) {
  M <- as.matrix(M)
  M[!is.finite(M)] <- NA_real_
  
  mu <- rowMeans(M, na.rm = TRUE)
  sd <- rowSds(M, na.rm = TRUE)
  
  keep <- is.finite(sd) & sd > 0
  
  Z <- sweep(M[keep, , drop = FALSE], 1, mu[keep], "-")
  Z <- sweep(Z, 1, sd[keep], "/")
  
  Z[!is.finite(Z)] <- NA_real_
  Z
}

HPTM_filtered_z              <- zscore_rows(HPTM_filtered)
Liver_HPTM_filtered_z        <- zscore_rows(Liver_HPTM_filtered)
midgut_proteome_filtered_z   <- zscore_rows(midgut_proteome_filtered)
liver_proteome_filtered_z    <- zscore_rows(liver_proteome_filtered)
MicrobiotaDigesta_filtered_z <- zscore_rows(MicrobiotaDigesta_filtered)
MicrobiotaMucus_filtered_z   <- zscore_rows(MicrobiotaMucus_filtered)

## Sanity check: after z-scoring, feature means should be ~0 and SDs ~1
zscore_check <- data.frame(
  view = c(
    "Mid-Gut hPTMs",
    "Liver hPTMs",
    "Mid-Gut Proteome",
    "Liver Proteome",
    "Digesta Microbiota",
    "Mucus Microbiota"
  ),
  median_feature_mean = c(
    median(rowMeans(HPTM_filtered_z, na.rm = TRUE), na.rm = TRUE),
    median(rowMeans(Liver_HPTM_filtered_z, na.rm = TRUE), na.rm = TRUE),
    median(rowMeans(midgut_proteome_filtered_z, na.rm = TRUE), na.rm = TRUE),
    median(rowMeans(liver_proteome_filtered_z, na.rm = TRUE), na.rm = TRUE),
    median(rowMeans(MicrobiotaDigesta_filtered_z, na.rm = TRUE), na.rm = TRUE),
    median(rowMeans(MicrobiotaMucus_filtered_z, na.rm = TRUE), na.rm = TRUE)
  ),
  median_feature_sd = c(
    median(rowSds(HPTM_filtered_z, na.rm = TRUE), na.rm = TRUE),
    median(rowSds(Liver_HPTM_filtered_z, na.rm = TRUE), na.rm = TRUE),
    median(rowSds(midgut_proteome_filtered_z, na.rm = TRUE), na.rm = TRUE),
    median(rowSds(liver_proteome_filtered_z, na.rm = TRUE), na.rm = TRUE),
    median(rowSds(MicrobiotaDigesta_filtered_z, na.rm = TRUE), na.rm = TRUE),
    median(rowSds(MicrobiotaMucus_filtered_z, na.rm = TRUE), na.rm = TRUE)
  ),
  n_features = c(
    nrow(HPTM_filtered_z),
    nrow(Liver_HPTM_filtered_z),
    nrow(midgut_proteome_filtered_z),
    nrow(liver_proteome_filtered_z),
    nrow(MicrobiotaDigesta_filtered_z),
    nrow(MicrobiotaMucus_filtered_z)
  )
)

zscore_check

# Align all views to the same non-fasted sample universe

all_samples <- rownames(meta_filtered_with_traits)

align_to_samples <- function(M, all_samples) {
  current <- colnames(M)
  common  <- intersect(current, all_samples)
  
  M_sub <- M[, common, drop = FALSE]
  
  M_aligned <- matrix(
    NA_real_,
    nrow = nrow(M_sub),
    ncol = length(all_samples),
    dimnames = list(rownames(M_sub), all_samples)
  )
  
  M_aligned[, common] <- M_sub
  M_aligned
}

HPTM_aligned              <- align_to_samples(HPTM_filtered, all_samples)
Liver_HPTM_aligned        <- align_to_samples(Liver_HPTM_filtered, all_samples)

midgut_proteome_aligned   <- align_to_samples(midgut_proteome_filtered, all_samples)
liver_proteome_aligned    <- align_to_samples(liver_proteome_filtered, all_samples)

MicrobiotaDigesta_aligned <- align_to_samples(MicrobiotaDigesta_filtered, all_samples)
MicrobiotaMucus_aligned   <- align_to_samples(MicrobiotaMucus_filtered, all_samples)

Host_Traits_aligned       <- align_to_samples(Host_Traits_filtered, all_samples)

meta_common <- meta_filtered_with_traits[all_samples, , drop = FALSE]

stopifnot(identical(colnames(HPTM_aligned), all_samples))
stopifnot(identical(colnames(Liver_HPTM_aligned), all_samples))
stopifnot(identical(colnames(midgut_proteome_aligned), all_samples))
stopifnot(identical(colnames(liver_proteome_aligned), all_samples))
stopifnot(identical(colnames(MicrobiotaDigesta_aligned), all_samples))
stopifnot(identical(colnames(MicrobiotaMucus_aligned), all_samples))
stopifnot(identical(colnames(Host_Traits_aligned), all_samples))
stopifnot(identical(rownames(meta_common), all_samples))


# Z-scored filtered data, aligned after z-scoring
HPTM_aligned_z              <- align_to_samples(HPTM_filtered_z, all_samples)
Liver_HPTM_aligned_z        <- align_to_samples(Liver_HPTM_filtered_z, all_samples)

midgut_proteome_aligned_z   <- align_to_samples(midgut_proteome_filtered_z, all_samples)
liver_proteome_aligned_z    <- align_to_samples(liver_proteome_filtered_z, all_samples)

MicrobiotaDigesta_aligned_z <- align_to_samples(MicrobiotaDigesta_filtered_z, all_samples)
MicrobiotaMucus_aligned_z   <- align_to_samples(MicrobiotaMucus_filtered_z, all_samples)


# Missingness diagnostics

check_missing_samples <- function(mat, name) {
  sample_ids <- colnames(mat)
  
  all_na     <- colSums(!is.na(mat)) == 0
  partial_na <- colSums(!is.na(mat)) > 0 & colSums(is.na(mat)) > 0
  no_na      <- colSums(is.na(mat)) == 0
  
  cat("====", name, "====\n")
  cat("Total samples:", length(sample_ids), "\n")
  cat("Fully NA samples:     ", sum(all_na), "\n")
  cat("Partially NA samples: ", sum(partial_na), "\n")
  cat("Fully observed:       ", sum(no_na), "\n\n")
  
  if (sum(all_na) > 0) {
    cat("Samples fully NA in", name, ":\n")
    print(sample_ids[all_na])
    cat("\n")
  }
}

check_missing_samples(HPTM_aligned, "Midgut hPTM")
check_missing_samples(Liver_HPTM_aligned, "Liver hPTM")
check_missing_samples(midgut_proteome_aligned, "Midgut proteome")
check_missing_samples(liver_proteome_aligned, "Liver proteome")
check_missing_samples(MicrobiotaDigesta_aligned, "Microbiota Digesta")
check_missing_samples(MicrobiotaMucus_aligned, "Microbiota Mucus")
check_missing_samples(Host_Traits_aligned, "Host traits")


# Build MOFA input object

rownames(midgut_proteome_aligned) <-
  paste0("Midgut_", rownames(midgut_proteome_aligned))

rownames(liver_proteome_aligned) <-
  paste0("Liver_", rownames(liver_proteome_aligned))

rownames(MicrobiotaDigesta_aligned) <-
  paste0("Digesta_", rownames(MicrobiotaDigesta_aligned))

rownames(MicrobiotaMucus_aligned) <-
  paste0("Mucus_", rownames(MicrobiotaMucus_aligned))

rownames(HPTM_aligned) <-
  paste0("Midgut_", rownames(HPTM_aligned))

rownames(Liver_HPTM_aligned) <-
  paste0("Liver_", rownames(Liver_HPTM_aligned))

mofa_data_all <- list(
  "Midgut hPTMs"      = HPTM_aligned,
  "Liver hPTMs"        = Liver_HPTM_aligned,
  "Midgut Proteome"   = midgut_proteome_aligned,
  "Liver Proteome"     = liver_proteome_aligned,
  "Digesta Microbiota" = MicrobiotaDigesta_aligned,
  "Mucus Microbiota"   = MicrobiotaMucus_aligned
)

rownames(midgut_proteome_aligned_z) <-
  paste0("Midgut_", rownames(midgut_proteome_aligned_z))

rownames(liver_proteome_aligned_z) <-
  paste0("Liver_", rownames(liver_proteome_aligned_z))

rownames(MicrobiotaDigesta_aligned_z) <-
  paste0("Digesta_", rownames(MicrobiotaDigesta_aligned_z))

rownames(MicrobiotaMucus_aligned_z) <-
  paste0("Mucus_", rownames(MicrobiotaMucus_aligned_z))

rownames(HPTM_aligned_z) <-
  paste0("Midgut_", rownames(HPTM_aligned_z))

rownames(Liver_HPTM_aligned_z) <-
  paste0("Liver_", rownames(Liver_HPTM_aligned_z))
mofa_data_all_zscore <- list(
  "Midgut hPTMs"      = HPTM_aligned_z,
  "Liver hPTMs"        = Liver_HPTM_aligned_z,
  "Midgut Proteome"   = midgut_proteome_aligned_z,
  "Liver Proteome"     = liver_proteome_aligned_z,
  "Digesta Microbiota" = MicrobiotaDigesta_aligned_z,
  "Mucus Microbiota"   = MicrobiotaMucus_aligned_z
)


# Save mofa views
saveRDS(
  mofa_data_all_zscore,
  file = "views.rds"
)

# Save common metadata
saveRDS(
  meta_common,
  file = "covariates.rds"
)


# MultiPower input: same filtered/aligned matrices, before feature-wise z-scoring
saveRDS(
  mofa_data_all,
  file = "views_not_zscored_for_MultiPower.rds"
)