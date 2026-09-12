library(reticulate)

py <- "C:/mofa_venv/Scripts/python.exe"
stopifnot(file.exists(py))
use_python(py, required = TRUE)
py_config()

suppressPackageStartupMessages({
  library(MOFA2)
  library(tibble)
  library(dplyr)
})

# Load data
mofa_data   <- readRDS("views.rds")
meta_common <- readRDS("covariates.rds")
# Check metadata column names
colnames(meta_common)
samples <- colnames(mofa_data[[1]])
stopifnot(identical(rownames(meta_common), samples))

meta_common$diet <- as.factor(meta_common$diet)

meta_common2 <- meta_common
meta_common2$sample <- rownames(meta_common2)

# ---- Create MOFA object ----
mofa <- create_mofa(mofa_data)

samples_metadata(mofa) <- meta_common2

# Data options
data_opts <- get_default_data_options(mofa)
data_opts$scale_views <- TRUE

# Model options
model_opts <- get_default_model_options(mofa)
model_opts$num_factors <- 20
model_opts$spikeslab_factors <- FALSE
model_opts$spikeslab_weights <- TRUE

# Training options
train_opts <- get_default_training_options(mofa)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 85000
train_opts$maxiter <- 5000
train_opts$weight_views <- FALSE
train_opts$verbose <- TRUE

# Prepare MOFA model (no MEFISTO)
MOFAobject <- prepare_mofa(
  object = mofa,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

# Output directory
model_dir <- "mofa_models"
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

# Train
MOFAobject.trained <- run_mofa(
  MOFAobject,
  outfile = file.path(model_dir, "K20_setseed_85000.hdf5"),
  use_basilisk = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(model_dir, "sessionInfo.txt")
)
