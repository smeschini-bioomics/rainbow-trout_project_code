library(readr)
library(dplyr)
library(tidyr)

# change to digesta or simply""01_ASVs_relative_abundance_table" for global stat
input_dir <- "01_ASVs_relative_abundance_table/mucus"

summarize_abundance_file <- function(file_path, output_file = NULL) {
  
  dat <- read.csv(
    file_path,
    header = TRUE,
    sep = ",",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  if (colnames(dat)[1] != "Taxa") {
    stop(paste("First column must be named 'Taxa' in file:", file_path))
  }
  
  taxa <- dat$Taxa
  abund <- dat[, -1, drop = FALSE]
  abund[] <- lapply(abund, function(x) as.numeric(as.character(x)))
  abund_mat <- as.matrix(abund)
  
  stats_df <- data.frame(
    Taxa       = taxa,
    N_samples  = rowSums(!is.na(abund_mat)),
    Mean       = rowMeans(abund_mat, na.rm = TRUE),
    Median     = apply(abund_mat, 1, median, na.rm = TRUE),
    SD         = apply(abund_mat, 1, sd, na.rm = TRUE),
    Min        = apply(abund_mat, 1, min, na.rm = TRUE),
    Q1         = apply(abund_mat, 1, quantile, probs = 0.25, na.rm = TRUE),
    Q3         = apply(abund_mat, 1, quantile, probs = 0.75, na.rm = TRUE),
    Max        = apply(abund_mat, 1, max, na.rm = TRUE),
    IQR        = apply(abund_mat, 1, IQR, na.rm = TRUE),
    Sum        = rowSums(abund_mat, na.rm = TRUE),
    Nonzero_n  = rowSums(abund_mat > 0, na.rm = TRUE),
    Zero_n     = rowSums(abund_mat == 0, na.rm = TRUE),
    Prevalence = rowSums(abund_mat > 0, na.rm = TRUE) / rowSums(!is.na(abund_mat)),
    CV_percent = ifelse(
      rowMeans(abund_mat, na.rm = TRUE) == 0,
      NA,
      apply(abund_mat, 1, sd, na.rm = TRUE) / rowMeans(abund_mat, na.rm = TRUE) * 100
    ),
    stringsAsFactors = FALSE
  )
  
  numeric_cols <- sapply(stats_df, is.numeric)
  stats_df[numeric_cols] <- lapply(stats_df[numeric_cols], function(x) round(x, 6))
  
  if (!is.null(output_file)) {
    write.csv(stats_df, output_file, row.names = FALSE, quote = TRUE)
  }
  
  return(stats_df)
}

files <- list.files(input_dir, pattern = "_abund\\.csv$", full.names = TRUE)

for (f in files) {
  out_file <- file.path(
    input_dir,
    paste0(tools::file_path_sans_ext(basename(f)), "_descriptive_stats.csv")
  )
  
  cat("\nProcessing:", basename(f), "\n")
  res <- summarize_abundance_file(f, out_file)
  print(head(res, 5))
}

