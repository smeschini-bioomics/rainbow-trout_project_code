# Midgut hPTM covariance analysis by diet

# Libraries
library(Hmisc)
library(corrplot)

# Settings
INPUT_RDS <- "05_batch_correction_limma/H3_H4_5PTM_midgut_batch_corrected_limma.rds"
METADATA_FILE <- "STPN2309_metadata_all.tsv"

OUTPUT_DIR <- "08_a_covariance analysis_corr"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

SAMPLE_ID_COL <- "sampleID_STPN2309"
DIET_COL <- "diet"
DIET_LEVELS <- c("NC", "HC")

P_CUTOFF <- 0.05

TITLE_FAMILY <- "mono"
TITLE_SIZE <- 16
HCLUST_METHOD <- "complete"

# Check input files
if (!file.exists(INPUT_RDS)) {
  stop("Corrected hPTM RDS not found: ", INPUT_RDS)
}

if (!file.exists(METADATA_FILE)) {
  stop("Metadata TSV not found: ", METADATA_FILE)
}

# Load corrected hPTM matrix and metadata
rds_obj <- readRDS(INPUT_RDS)

if (!"ptm_batch_corrected_limma" %in% names(rds_obj)) {
  stop(
    "Object 'ptm_batch_corrected_limma' not found in: ",
    INPUT_RDS
  )
}

ptm_df <- as.data.frame(
  rds_obj$ptm_batch_corrected_limma,
  stringsAsFactors = FALSE
)

meta_df <- read.delim(
  METADATA_FILE,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (is.null(rownames(ptm_df))) {
  stop("Corrected hPTM matrix must contain sample IDs as rownames.")
}

required_cols <- c(SAMPLE_ID_COL, DIET_COL)

missing_cols <- setdiff(required_cols, names(meta_df))

if (length(missing_cols) > 0) {
  stop(
    "Missing metadata column(s): ",
    paste(missing_cols, collapse = ", ")
  )
}

# Align hPTM matrix and metadata
rownames(ptm_df) <- trimws(rownames(ptm_df))

meta_df[[SAMPLE_ID_COL]] <- trimws(
  as.character(meta_df[[SAMPLE_ID_COL]])
)

meta_df[[DIET_COL]] <- trimws(
  as.character(meta_df[[DIET_COL]])
)

if (anyDuplicated(rownames(ptm_df)) > 0) {
  stop("Duplicated sample IDs found in corrected hPTM matrix.")
}

if (anyDuplicated(meta_df[[SAMPLE_ID_COL]]) > 0) {
  stop("Duplicated sample IDs found in metadata.")
}

meta_df <- meta_df[
  meta_df[[DIET_COL]] %in% DIET_LEVELS,
  ,
  drop = FALSE
]

common_ids <- rownames(ptm_df)[
  rownames(ptm_df) %in% meta_df[[SAMPLE_ID_COL]]
]

if (length(common_ids) == 0) {
  stop("No matching sample IDs found between hPTM matrix and metadata.")
}

ptm_df <- ptm_df[
  common_ids,
  ,
  drop = FALSE
]

meta_df <- meta_df[
  match(common_ids, meta_df[[SAMPLE_ID_COL]]),
  ,
  drop = FALSE
]

rownames(meta_df) <- meta_df[[SAMPLE_ID_COL]]

stopifnot(
  identical(
    rownames(ptm_df),
    rownames(meta_df)
  )
)

meta_df[[DIET_COL]] <- factor(
  meta_df[[DIET_COL]],
  levels = DIET_LEVELS
)

ptm_df[] <- lapply(
  ptm_df,
  function(x) suppressWarnings(as.numeric(x))
)

if (anyDuplicated(colnames(ptm_df)) > 0) {
  stop("Duplicated PTM names found in corrected hPTM matrix.")
}

normalize_ptm <- function(x) {
  vapply(x, function(entry) {
    parts <- strsplit(entry, "\\|")[[1]]
    norm_parts <- vapply(parts, function(p) {
      g <- regmatches(p, regexec("^(H[34])(K|R)([0-9]{1,3})([A-Za-z0-9]+)$", p))[[1]]
      if (length(g) < 5 || is.na(g[1])) return(NA_character_)
      paste0(g[2], g[3], g[4], tolower(g[5]))
    }, character(1))
    if (any(is.na(norm_parts))) return(NA_character_)
    paste(norm_parts, collapse = "|")
  }, character(1), USE.NAMES = FALSE)
}

colnames(ptm_df) <- normalize_ptm(colnames(ptm_df))

# QC check — make sure nothing became NA
print(colnames(ptm_df))
if (anyDuplicated(colnames(ptm_df)) > 0) {
  stop("Duplicated PTM names found AFTER normalization — check regex/parsing.")
}


cat("Aligned matrix:", nrow(ptm_df), "samples x", ncol(ptm_df), "PTMs\n")
print(table(meta_df[[DIET_COL]]))

# Define one common PTM order using all samples
global_cor <- rcorr(
  as.matrix(ptm_df),
  type = "spearman"
)$r

global_cor_for_clust <- global_cor
global_cor_for_clust[is.na(global_cor_for_clust)] <- 0
diag(global_cor_for_clust) <- 1

global_hc <- hclust(
  as.dist(1 - global_cor_for_clust),
  method = HCLUST_METHOD
)

ordered_features <- colnames(global_cor_for_clust)[global_hc$order]

# Run correlation analysis for one diet
run_corrplot_by_group <- function(
    group_name,
    ptm_sub,
    ordered_features,
    output_dir
) {
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  cat("\nProcessing:", group_name, "\n")
  cat(nrow(ptm_sub), "samples x", ncol(ptm_sub), "PTMs\n")
  
  if (nrow(ptm_sub) < 3) {
    warning("Fewer than 3 samples available for: ", group_name)
    return(invisible(NULL))
  }
  
  cor_res <- rcorr(
    as.matrix(ptm_sub),
    type = "spearman"
  )
  
  cor_mat <- cor_res$r
  p_mat <- cor_res$P
  n_mat <- cor_res$n
  
  common_features <- ordered_features[
    ordered_features %in% colnames(cor_mat)
  ]
  
  cor_mat <- cor_mat[
    common_features,
    common_features,
    drop = FALSE
  ]
  
  p_mat <- p_mat[
    common_features,
    common_features,
    drop = FALSE
  ]
  
  n_mat <- n_mat[
    common_features,
    common_features,
    drop = FALSE
  ]
  
  keep_mat <- !is.na(cor_mat) &
    !is.na(p_mat) &
    p_mat < P_CUTOFF
  
  diag(keep_mat) <- FALSE
  
  p_plot_mat <- p_mat
  p_plot_mat[!keep_mat] <- 1
  
  n_displayed <- sum(
    keep_mat[upper.tri(keep_mat)],
    na.rm = TRUE
  )
  
  cat(
    "Displayed correlations with nominal p <",
    P_CUTOFF,
    ":",
    n_displayed,
    "\n"
  )
  
  write.table(
    cor_mat,
    file = file.path(
      output_dir,
      paste0("spearman_correlation_matrix_", group_name, ".tsv")
    ),
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
  
  write.table(
    p_mat,
    file = file.path(
      output_dir,
      paste0("spearman_pvalue_matrix_", group_name, ".tsv")
    ),
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
  
  write.table(
    n_mat,
    file = file.path(
      output_dir,
      paste0("spearman_pairwise_n_matrix_", group_name, ".tsv")
    ),
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
  
  plot_corr <- function() {
    
    corrplot(
      cor_mat,
      method = "pie",
      order = "original",
      col = colorRampPalette(
        c("#08519c", "white", "#8B0000")
      )(200),
      tl.col = "black",
      tl.cex = 1,
      p.mat = p_plot_mat,
      sig.level = P_CUTOFF,
      insig = "blank",
      bg = "grey95",
      diag = FALSE,
      mar = c(0, 0, 4, 0),
      title = ""
    )
  }
  
  pdf(
    file = file.path(
      output_dir,
      paste0("corrplot_ptm_", group_name, ".pdf")
    ),
    width = 6,
    height = 6
  )
  
  plot_corr()
  dev.off()
  
  tiff(
    filename = file.path(
      output_dir,
      paste0("corrplot_ptm_", group_name, ".tiff")
    ),
    width = 6,
    height = 6,
    units = "in",
    res = 1200,
    compression = "lzw"
  )
  
  plot_corr()
  dev.off()
  
  invisible(
    list(
      correlation = cor_mat,
      p_value = p_mat,
      pairwise_n = n_mat
    )
  )
}

# Run correlations by diet
results_list <- list()

for (grp in DIET_LEVELS) {
  
  idx <- which(meta_df[[DIET_COL]] == grp)
  
  if (length(idx) == 0) {
    warning("No samples found for: ", grp)
    next
  }
  
  results_list[[grp]] <- run_corrplot_by_group(
    group_name = grp,
    ptm_sub = ptm_df[idx, , drop = FALSE],
    ordered_features = ordered_features,
    output_dir = file.path(OUTPUT_DIR, grp)
  )
}

writeLines(
  capture.output(sessionInfo()),
  con = file.path(OUTPUT_DIR, "sessionInfo.txt")
)

cat("\nDone. Outputs saved in:", OUTPUT_DIR, "\n")
