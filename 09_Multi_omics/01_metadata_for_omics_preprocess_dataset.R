library(data.table)
library(dplyr)

# -----------------------------------------------------------------------------

# Input files

# -----------------------------------------------------------------------------

input_dir <- "input_files"

metadata_file <- file.path(input_dir, "STPN2309_metadata_all.tsv")
glucose_file  <- file.path(input_dir, "plasma_glucose.tsv")
dlactate_file <- file.path(input_dir, "plasma_dlactate.tsv")
llactate_file <- file.path(input_dir, "plasma_llactate.tsv")
bw_file       <- file.path(input_dir, "body_weight.tsv")
hsi_file      <- file.path(input_dir, "hepatosomatic_index.tsv")
methylome_file <- file.path(input_dir, "methylome_for_dirichlet.tsv")

output_file <- "metadata_for_omics_data.rds"

required_files <- c(
  metadata_file,
  glucose_file,
  dlactate_file,
  llactate_file,
  bw_file,
  hsi_file,
  methylome_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "The following input files were not found:\n",
    paste(missing_files, collapse = "\n")
  )
}

# -----------------------------------------------------------------------------

# Metadata: retain design factors only

# -----------------------------------------------------------------------------

metadata <- fread(metadata_file) |> as.data.frame()

rownames(metadata) <- metadata[[1]]

metadata <- metadata |>
  select(
    day_categ,
    diet,
    group,
    phase,
    diet_phase
  ) |>
  mutate(
    across(
      everything(),
      as.factor
    )
  )

if (anyDuplicated(rownames(metadata))) {
  stop("Duplicated sample IDs detected in metadata.")
}

# -----------------------------------------------------------------------------

# Host physiological traits

# -----------------------------------------------------------------------------

plasma_glucose <- fread(glucose_file) |> as.data.frame()

replicate_cols <- grep("^replicate_", colnames(plasma_glucose), value = TRUE)

if (length(replicate_cols) == 0) {
  stop("No glucose replicate columns matching '^replicate_' were found.")
}

plasma_glucose_clean <- plasma_glucose |>
  mutate(
    plasma_glucose = rowMeans(
      across(all_of(replicate_cols)),
      na.rm = TRUE
    )
  ) |>
  select(SampleID_STPN2309, plasma_glucose)

HSI_clean <- fread(hsi_file) |>
  as.data.frame() |>
  select(SampleID_STPN2309, HSI)

BW_clean <- fread(bw_file) |>
  as.data.frame() |>
  select(SampleID_STPN2309, body_weight)

plasma_dlactate_clean <- fread(dlactate_file) |>
  as.data.frame() |>
  transmute(
    SampleID_STPN2309 = sampleID_STPN2309,
    plasma_DLactate = Plasma_D_Lactate
  )

plasma_llactate_clean <- fread(llactate_file) |>
  as.data.frame() |>
  transmute(
    SampleID_STPN2309 = sampleID_STPN2309,
    plasma_LLactate = Plasma_L_Lactate
  )

host_traits_merged <- list(
  plasma_glucose_clean,
  HSI_clean,
  BW_clean,
  plasma_dlactate_clean,
  plasma_llactate_clean
) |>
  Reduce(
    f = function(x, y) full_join(x, y, by = "SampleID_STPN2309")
  )

if (anyDuplicated(host_traits_merged$SampleID_STPN2309)) {
  stop("Duplicated sample IDs detected after merging host traits.")
}

host_traits_all <- as.data.frame(host_traits_merged)

rownames(host_traits_all) <- host_traits_all$SampleID_STPN2309
host_traits_all$SampleID_STPN2309 <- NULL

host_traits_all <- as.matrix(host_traits_all)

storage.mode(host_traits_all) <- "numeric"


# -----------------------------------------------------------------------------

# Midgut DNA cytosine-variant composition

# -----------------------------------------------------------------------------

methylome <- fread(methylome_file) |> as.data.frame()

rownames(methylome) <- methylome[[1]]

methylome <- methylome |>
  select(
    dC_rel,
    mdC_rel,
    hmdC_rel
  ) |>
  as.matrix()

storage.mode(methylome) <- "numeric"

if (anyDuplicated(rownames(methylome))) {
  stop("Duplicated sample IDs detected in methylome data.")
}

# -----------------------------------------------------------------------------

# Assemble and save

# -----------------------------------------------------------------------------

data <- list(
  metadata = metadata,
  host_traits = host_traits_all,
  Liver_DNA_methylation = methylome
)

message("Dimensions of exported objects:")
print(lapply(data, dim))

saveRDS(data, output_file)

message("Saved: ", output_file)
