suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(rlang)
  library(readr)
})

# Input and output folders
rds_dir <- "RDS"
out_dir <- file.path("results", "sample_size")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Find all gene plot objects
rds_files <- list.files(
  rds_dir,
  pattern = ".rds$",
  full.names = TRUE
)

if (length(rds_files) == 0) {
  stop("No '*_plot.rds' files were found in: ", rds_dir)
}

# Find the plotted y-variable in each ggplot object
get_yvar <- function(p) {
  if (!is.null(p$mapping$y)) {
    return(as_name(p$mapping$y))
  }
  
  for (i in seq_along(p$layers)) {
    if (!is.null(p$layers[[i]]$mapping$y)) {
      return(as_name(p$layers[[i]]$mapping$y))
    }
  }
  
  stop("Could not detect the y variable in one RDS plot.")
}

# Find a column name safely
get_colname <- function(df, candidates) {
  found <- intersect(candidates, colnames(df))
  
  if (length(found) == 0) {
    stop(
      "None of these columns were found: ",
      paste(candidates, collapse = ", ")
    )
  }
  
  found[1]
}

# Extract sample size from every RDS plot
sample_size_table <- map_dfr(rds_files, function(rds_file) {
  
  p <- readRDS(rds_file)
  
  if (is.null(p$data)) {
    stop("No data were stored inside: ", basename(rds_file))
  }
  
  plot_data <- as.data.frame(p$data)
  
  gene <- basename(rds_file) |>
    tools::file_path_sans_ext() |>
    sub("_plot$", "", x = _)
  
  yvar <- get_yvar(p)
  
  diet_col <- get_colname(plot_data, c("diet", "Diet"))
  day_col  <- get_colname(plot_data, c("day", "Day", "day_label", "day_num"))
  
  plot_data %>%
    filter(
      !is.na(.data[[yvar]]),
      !is.na(.data[[diet_col]]),
      !is.na(.data[[day_col]])
    ) %>%
    count(
      gene = gene,
      day = as.character(.data[[day_col]]),
      diet = as.character(.data[[diet_col]]),
      name = "n"
    )
})

# Order days
sample_size_table <- sample_size_table %>%
  mutate(
    day_num = as.numeric(gsub("[^0-9]", "", day))
  ) %>%
  arrange(gene, day_num, diet) %>%
  select(gene, day, diet, n)

print(sample_size_table)

# Save one long-format CSV table
write_csv(
  sample_size_table,
  file.path(out_dir, "sample_size.csv")
)
