# =============================================================================
# MultiGroupPower: HC vs NC within each sampling day
# =============================================================================

suppressPackageStartupMessages({
  library(MultiPower)
  library(FDRsampsize)
  library(lpSolve)
})

set.seed(42)

# -----------------------------------------------------------------------------
# Input and output
# -----------------------------------------------------------------------------

views_file <- "views_not_zscored_for_MultiPower.rds"
covariates_file <- "covariates.rds"

results_root <- "MultiGroupPower"
figures_dir <- file.path(results_root, "figures")
tables_dir <- file.path(results_root, "tables")
rds_dir <- file.path(results_root, "rds")

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(views_file))
stopifnot(file.exists(covariates_file))

# -----------------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------------

type1 <- rep(2, 6)

omic_power_target <- 0.60
average_power_target <- 0.80
fdr_target <- 0.05
max_sample_size <- 2000
omic_cost <- 1

# -----------------------------------------------------------------------------
# Colours
# -----------------------------------------------------------------------------

omic_names <- c(
  "Midgut hPTMs",
  "Liver hPTMs",
  "Midgut Proteome",
  "Liver Proteome",
  "Digesta Microbiota",
  "Mucus Microbiota"
)

omic_colors <- c(
  "#fb8072",
  "#e66101",
  "#1f78b4",
  "#ffd92f",
  "#2ca25f",
  "#7570b3"
)

names(omic_colors) <- omic_names

# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------

statdata <- readRDS(views_file)
metadata <- readRDS(covariates_file)

if (!is.list(statdata) || is.null(names(statdata))) {
  stop("views.rds must contain a named list of omics matrices.")
}

if (is.null(rownames(metadata))) {
  stop("covariates.rds must have sample IDs as row names.")
}

if (!all(c("diet", "day_categ") %in% colnames(metadata))) {
  stop("covariates.rds must contain columns: diet and day_categ.")
}

day_order <- c(
  "day_01", "day_02", "day_03", "day_04",
  "day_10", "day_15", "day_22"
)

metadata <- metadata[
  metadata$diet %in% c("NC", "HC") &
    metadata$day_categ %in% day_order,
  ,
  drop = FALSE
]

# -----------------------------------------------------------------------------
# Align each omics matrix with metadata
#
# No filtering, imputation, transformation, or manual parameter estimation.
# The matrices are passed to MultiGroupPower as provided.
# -----------------------------------------------------------------------------

statdesign <- vector("list", length(statdata))
names(statdesign) <- names(statdata)

for (omic in names(statdata)) {
  
  statdata[[omic]] <- as.matrix(statdata[[omic]])
  storage.mode(statdata[[omic]]) <- "numeric"
  
  sample_ids <- colnames(statdata[[omic]])
  sample_ids <- sample_ids[sample_ids %in% rownames(metadata)]
  
  statdata[[omic]] <- statdata[[omic]][, sample_ids, drop = FALSE]
  
  statdesign[[omic]] <- paste(
    metadata[sample_ids, "diet"],
    metadata[sample_ids, "day_categ"],
    sep = "_"
  )
}

# Keep only views present in the colour vector
keep_omics <- names(statdata) %in% names(omic_colors)

statdata <- statdata[keep_omics]
statdesign <- statdesign[keep_omics]

type1 <- rep(2, length(statdata))
names(type1) <- names(statdata)

omic_colors <- omic_colors[names(statdata)]

# -----------------------------------------------------------------------------
# Requested within-day HC vs NC comparisons
# -----------------------------------------------------------------------------

comparisons <- rbind(
  paste("NC", day_order, sep = "_"),
  paste("HC", day_order, sep = "_")
)
colnames(comparisons) <- c(
  "Day 1",
  "Day 2",
  "Day 3",
  "Day 4",
  "Day 10",
  "Day 15",
  "Day 22"
)

write.csv(
  as.data.frame(comparisons),
  file.path(tables_dir, "requested_day_specific_comparisons.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    data = statdata,
    groups = statdesign,
    comparisons = comparisons
  ),
  file.path(rds_dir, "MultiGroupPower_inputs.rds")
)

# -----------------------------------------------------------------------------
# Equal-size MultiGroupPower
#
# Balanced experimental design:
# HC vs NC at each sampling day with equal replication.
# -----------------------------------------------------------------------------

pdf(
  file.path(figures_dir, "MultiGroupPower_equal_size_package_plots.pdf"),
  width = 10,
  height = 7
)

results_equal <- MultiGroupPower(
  data = statdata,
  groups = statdesign,
  type = type1,
  comparisons = comparisons,
  omicPower = omic_power_target,
  averagePower = average_power_target,
  fdr = fdr_target,
  cost = omic_cost,
  equalSize = TRUE,
  max.size = max_sample_size,
  omicCol = omic_colors,
  powerPlots = TRUE,
  summaryPlot = TRUE
)

dev.off()

save(
  results_equal,
  file = file.path(
    rds_dir,
    "MultiGroupPower_equal_size_results.RData"
  )
)

write.csv(
  results_equal$GlobalSummary,
  file.path(
    tables_dir,
    "MultiGroupPower_equal_size_global_summary.csv"
  ),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Post MultiGroupPower analysis
#
# Evaluate the detectable effect sizes using the actual experimental
# replication (12 fish per diet per sampling day).
# -----------------------------------------------------------------------------

pdf(
  file.path(figures_dir, "PostMultiPower_equal_size_plots.pdf"),
  width = 10,
  height = 7
)

post_results <- lapply(
  names(results_equal)[1:7],
  function(day_comparison) {
    
    postMultiPower(
      optResults = results_equal[[day_comparison]],
      max.size = 12,
      omicCol = omic_colors
    )
    
  }
)

dev.off()

names(post_results) <- names(results_equal)[1:7]


# -----------------------------------------------------------------------------
# Save results
# -----------------------------------------------------------------------------

save(
  post_results,
  file = file.path(
    rds_dir,
    "PostMultiPower_equal_size_results.RData"
  )
)


# -----------------------------------------------------------------------------
# Export tables for each day
# -----------------------------------------------------------------------------

for (day_name in names(post_results)) {
  
  write.csv(
    post_results[[day_name]]$Power,
    file.path(
      tables_dir,
      paste0("PostMultiPower_", day_name, "_Power.csv")
    ),
    row.names = TRUE
  )
  
  write.csv(
    post_results[[day_name]]$SampleSize,
    file.path(
      tables_dir,
      paste0("PostMultiPower_", day_name, "_SampleSize.csv")
    ),
    row.names = TRUE
  )
  
  write.csv(
    post_results[[day_name]]$NumFeat,
    file.path(
      tables_dir,
      paste0("PostMultiPower_", day_name, "_NumFeat.csv")
    ),
    row.names = TRUE
  )
  
}


# -----------------------------------------------------------------------------
# Session information
# -----------------------------------------------------------------------------

writeLines(
  capture.output(sessionInfo()),
  file.path(results_root, "sessionInfo.txt")
)

message("MultiGroupPower completed successfully.")
message("PostMultiPower completed successfully.")
message("Results written to: ", normalizePath(results_root))