suppressPackageStartupMessages({
  library(MOFA2)
  library(ggplot2)
})

# -------------------------------------------------------------------------
# Input / output
# -------------------------------------------------------------------------

model_dir  <- "mofa_models"
output_dir <- "elbo_mofa_models"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# Load models
# -------------------------------------------------------------------------

model_files <- sort(list.files(
  model_dir,
  pattern = "\\.hdf5$",
  full.names = TRUE
))

if (length(model_files) < 2) {
  stop("At least two .hdf5 models are required in: ", model_dir)
}

models <- lapply(model_files, MOFA2::load_model)
names(models) <- tools::file_path_sans_ext(basename(model_files))

# -------------------------------------------------------------------------
# Raw ELBO table
# -------------------------------------------------------------------------

elbo_data <- data.frame(
  model = names(models),
  raw_ELBO = vapply(models, MOFA2::get_elbo, numeric(1)),
  stringsAsFactors = FALSE
)

# Extract seed from names such as K20_setseed_42
elbo_data$seed <- suppressWarnings(
  as.integer(sub(".*setseed_([0-9]+).*", "\\1", elbo_data$model))
)

if (anyNA(elbo_data$seed)) {
  stop(
    "Could not extract seeds from model names. ",
    "Expected names such as: K20_setseed_42"
  )
}

elbo_data <- elbo_data[order(elbo_data$seed), ]
elbo_data$model <- factor(elbo_data$model, levels = elbo_data$model)

best_row <- which.max(elbo_data$raw_ELBO)
best_model <- as.character(elbo_data$model[best_row])

write.csv(
  elbo_data,
  file.path(output_dir, "ELBO_comparison_data.csv"),
  row.names = FALSE
)

writeLines(
  c(
    paste("Best model:", best_model),
    paste(
      "Raw ELBO:",
      format(elbo_data$raw_ELBO[best_row], big.mark = ",")
    ),
    "Selection criterion: highest raw ELBO."
  ),
  file.path(output_dir, "best_model.txt")
)

# -------------------------------------------------------------------------
# Raw ELBO plot
# -------------------------------------------------------------------------

elbo_range <- range(elbo_data$raw_ELBO)
elbo_pad <- max(diff(elbo_range) * 0.10, 1)

p_elbo <- ggplot(
  elbo_data,
  aes(x = model, y = raw_ELBO)
) +
  geom_point(size = 3) +
  geom_point(
    data = elbo_data[best_row, , drop = FALSE],
    shape = 18,
    size = 4
  ) +
  geom_text(
    aes(
      label = format(raw_ELBO, big.mark = ",", trim = TRUE)
    ),
    vjust = -0.8,
    size = 3
  ) +
  coord_cartesian(
    ylim = elbo_range + c(-elbo_pad, elbo_pad)
  ) +
  scale_y_continuous(
    labels = scales::label_number(big.mark = ",")
  ) +
  labs(
    x = "MOFA model / random initialisation seed",
    y = "Final raw ELBO\n(higher / less negative is better)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    plot.margin = margin(12, 12, 12, 20)
  )

ggsave(
  file.path(output_dir, "ELBO_raw_by_seed.pdf"),
  p_elbo,
  width = 8,
  height = 5
)

ggsave(
  file.path(output_dir, "ELBO_raw_by_seed.tiff"),
  p_elbo,
  width = 8,
  height = 5,
  dpi = 1200,
  compression = "lzw"
)

# -------------------------------------------------------------------------
# Factor reproducibility heatmap
# -------------------------------------------------------------------------

# `silent = TRUE` returns the pheatmap object without drawing it first.
factor_heatmap <- MOFA2::compare_factors(
  models,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize = 10,
  show_rownames = FALSE,
  show_colnames = FALSE,
  silent = TRUE
)

save_factor_heatmap <- function(filename, type = c("pdf", "tiff")) {
  
  type <- match.arg(type)
  
  if (type == "pdf") {
    grDevices::pdf(
      filename,
      width = 8,
      height = 8,
      useDingbats = FALSE
    )
  } else {
    grDevices::tiff(
      filename,
      width = 8,
      height = 8,
      units = "in",
      res = 1200,
      compression = "lzw"
    )
  }
  
  grid::grid.newpage()
  grid::grid.draw(factor_heatmap$gtable)
  grDevices::dev.off()
}

save_factor_heatmap(
  file.path(output_dir, "factor_comparison.pdf"),
  type = "pdf"
)

save_factor_heatmap(
  file.path(output_dir, "factor_comparison.tiff"),
  type = "tiff"
)

# -------------------------------------------------------------------------
# Session information
# -------------------------------------------------------------------------

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "sessionInfo.txt")
)

message(
  "Done.\n",
  "Best model: ", best_model, "\n",
  "Output folder: ", output_dir
)