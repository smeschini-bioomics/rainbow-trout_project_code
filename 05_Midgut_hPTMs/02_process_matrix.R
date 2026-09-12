library(readr)
library(tidyverse)
library(data.table)
library(stringr)

## import mid-gut metadata 
metadata <- fread("STPN2309_metadata_all.tsv") # metadata
dim(metadata)
# Set the first column as rownames
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata[[1]]
metadata <- metadata[, -1]
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

## histone gut PTMs dataset import
histones_PTM <- read.csv("01_output_files/TroutGut_hPTM_H3H4.csv", sep = ",", check.names = FALSE)
# Set the first column as rownames
histones_PTM <- as.data.frame(histones_PTM)
rownames(histones_PTM) <- histones_PTM[[1]]
histones_PTM <- histones_PTM[, -1]
class(histones_PTM)
dim(histones_PTM)
rownames(histones_PTM)
histones_PTM <- as.matrix(histones_PTM)
dim(histones_PTM)
histones_PTM <- t(histones_PTM)
class(histones_PTM)
dim(histones_PTM)
class(histones_PTM)
dim(histones_PTM)
rownames(histones_PTM)


# Apply to your target matrix, e.g. H34_mod or histones_PTM
cn <- colnames(histones_PTM)  # or histones_PTM

# Capture: Histone, position, residue, modification
m <- str_match(cn, "^(H[0-9AB]+)\\s+\\[(\\d+)\\]\\s+\\(([A-Z])\\)\\s+([A-Za-z0-9]+)$")
# m[,2]=Histone, m[,3]=Position, m[,4]=Residue, m[,5]=Modification

# Build compact names: H3K9Ac
new_names <- ifelse(
  is.na(m[,1]),
  cn,                                    # keep original if no match
  paste0(m[,2], m[,4], m[,3], m[,5])     # e.g., H3 + K + 9 + Ac
)

# Optional: ensure uniqueness if collisions occur
if (any(duplicated(new_names))) new_names <- make.unique(new_names)

colnames(histones_PTM) <- new_names


# Collapse identical columns (histone PTMS that share same intensities); merge their original names with "|"
collapse_identical_cols <- function(X, sep = "|", tol = NULL) {
  X <- as.matrix(X)
  Xk <- if (is.null(tol)) X else round(X, tol)   # optional rounding
  
  # build a stable key per column (explicit NA token)
  key <- apply(Xk, 2, function(v)
    paste(ifelse(is.na(v), "NA",
                 format(v, scientific = FALSE, digits = 15, trim = TRUE)),
          collapse = "\u241F"))
  
  keep <- !duplicated(key)
  X_out <- X[, keep, drop = FALSE]
  
  # map: key -> all original column names
  groups <- split(colnames(X), key)
  
  # merged names for kept columns
  new_names <- vapply(colnames(X_out), function(nm) {
    k <- key[match(nm, colnames(X))]
    paste(groups[[k]], collapse = sep)
  }, character(1))
  
  colnames(X_out) <- new_names
  X_out
}


histones_PTM_collapsed <- collapse_identical_cols(histones_PTM, sep = "|", tol = 6)




# quick checks
cat("Before:", ncol(histones_PTM), " After:", ncol(histones_PTM_collapsed), "\n")
head(colnames(histones_PTM_collapsed))

stopifnot("sampleID_hPTMs_gut" %in% colnames(metadata))

# 1) Build mapping: SampleID_proteomics -> desired new name (metadata rowname)
map_old_to_new <- setNames(rownames(metadata), metadata$sampleID_hPTMs_gut)

# 2) Create new rownames for asv_mat_t using the mapping
old <- rownames(histones_PTM_collapsed)
new <- unname(map_old_to_new[old])  # will be NA if not found

# 3) Report any unmatched or duplicate targets (diagnostics)
unmatched <- old[is.na(new)]
if (length(unmatched)) {
  message("IDs in histones_PTM_collapsed not found in metadata$sampleID_hPTMs: ",
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
histones_PTM_renamed_collapsed <- histones_PTM_collapsed
rownames(histones_PTM_renamed_collapsed) <- new_fallback
# (optional) quick check
cat("Renamed rows:", sum(!is.na(new)), " / ", length(old), "\n")


# starting point: histones_PTM_renamed  (samples x features matrix)

# 1) Keep only columns whose names contain H3 or H4
cn <- colnames(histones_PTM)
sel_H34 <- grepl("^H[34]", cn, ignore.case = TRUE)
H3_H4_PTM <- histones_PTM[, sel_H34, drop = FALSE]

# 2) From that, keep only columns containing Ac, Me2, Me3, Cr, or La
mods <- c("Ac", "Me2","Me3","Cr","La")
rx_mod <- paste0("(", paste(mods, collapse = "|"), ")")
sel_mod <- grepl(rx_mod, colnames(H3_H4_PTM), ignore.case = TRUE)
H3_H4_5PTM_midgut <- H3_H4_PTM[, sel_mod, drop = FALSE]

# Quick checks
cat("Cols H3/H4:", sum(sel_H34), "\n")
cat("Cols H3/H4 with mods:", sum(sel_mod), "\n")

# collapse same signal
H3_H4_PTM_collapsed <- collapse_identical_cols(H3_H4_PTM, sep = "|", tol = 6)
# quick checks
cat("Before:", ncol(H3_H4_PTM), " After:", ncol(H3_H4_PTM_collapsed), "\n")
head(colnames(H3_H4_PTM_collapsed))

H3_H4_5PTM_midgut_collapsed <- collapse_identical_cols(H3_H4_5PTM_midgut, sep = "|", tol = 6)
# Remove excluded PTM before downstream filtering and export
FEATURES_TO_DROP <- c("H4K12Cr|H4K16La|H4R17Me2")

features_found <- intersect(
  FEATURES_TO_DROP,
  colnames(H3_H4_5PTM_midgut_collapsed)
)

if (length(features_found) > 0) {
  message(
    "Dropping selected PTM(s): ",
    paste(features_found, collapse = ", ")
  )
  
  H3_H4_5PTM_midgut_collapsed <-
    H3_H4_5PTM_midgut_collapsed[
      ,
      !colnames(H3_H4_5PTM_midgut_collapsed) %in% features_found,
      drop = FALSE
    ]
} else {
  warning(
    "FEATURES_TO_DROP was not found: ",
    paste(FEATURES_TO_DROP, collapse = ", ")
  )
}

# quick checks
cat("Before:", ncol(H3_H4_5PTM_midgut), " After:", ncol(H3_H4_5PTM_midgut_collapsed), "\n")
head(colnames(H3_H4_5PTM_midgut_collapsed))


# first, sum the number of missing values par varaible
sum.na.per.var <- apply(H3_H4_5PTM_midgut_collapsed, 2, function(x){sum(is.na(x))})
# a simple plot to show the NA rate par variable
plot(sum.na.per.var/nrow(H3_H4_5PTM_midgut_collapsed), type ='h',
     xlab = 'variable index', ylab = 'NA rate', main = 'NA rate  per variable')
# these variables could be removed according to a threshold, e.g. 20%
H3_H4_5PTM_midgut_collapsed_filtered <- H3_H4_5PTM_midgut_collapsed
remove.var <- which(sum.na.per.var/nrow(H3_H4_5PTM_midgut_collapsed_filtered) > 0.30)
H3_H4_5PTM_midgut_collapsed_filtered <- H3_H4_5PTM_midgut_collapsed_filtered[, -c(remove.var)]
# chack dimension of the data
dim(H3_H4_5PTM_midgut_collapsed_filtered)
#calculate the proportion of missing values remaining
sum(is.na(H3_H4_5PTM_midgut_collapsed_filtered)) / (length(H3_H4_5PTM_midgut_collapsed_filtered)) # proportion of missing value
# number of cells with NA
sum(is.na(H3_H4_5PTM_midgut_collapsed_filtered))
# rename sample 
stopifnot("sampleID_hPTMs_gut" %in% colnames(metadata))

# 1) Build mapping: SampleID_proteomics -> desired new name (metadata rowname)
map_old_to_new <- setNames(rownames(metadata), metadata$sampleID_hPTMs_gut)

# 2) Create new rownames for asv_mat_t using the mapping
old <- rownames(H3_H4_5PTM_midgut_collapsed_filtered)
new <- unname(map_old_to_new[old])  # will be NA if not found

# 3) Report any unmatched or duplicate targets (diagnostics)
unmatched <- old[is.na(new)]
if (length(unmatched)) {
  message("IDs in histones_PTM_collapsed not found in metadata$sampleID_hPTMs: ",
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
H3_H4_5PTM_midgut_collapsed_filtered_renamed <- H3_H4_5PTM_midgut_collapsed_filtered
rownames(H3_H4_5PTM_midgut_collapsed_filtered_renamed) <- new_fallback
# (optional) quick check
cat("Renamed rows:", sum(!is.na(new)), " / ", length(old), "\n")
# number of cells with NA
sum(is.na(H3_H4_5PTM_midgut_collapsed_filtered_renamed))


dim(H3_H4_5PTM_midgut_collapsed_filtered_renamed)
H3_H4_5PTM_midgut_collapsed_filtered_renamed <- H3_H4_5PTM_midgut_collapsed_filtered_renamed[, -1]
dim(H3_H4_5PTM_midgut_collapsed_filtered_renamed)

# Output directory
OUTDIR <- "02_process_matrix"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

saveRDS(
  H3_H4_5PTM_midgut_collapsed_filtered_renamed,
  file = file.path(OUTDIR, "H3_H4_5PTM_midgut.rds")
)

