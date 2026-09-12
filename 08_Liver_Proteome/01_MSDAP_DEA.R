library(msdap)

# 1) Load data files, output from upstream raw data processor and the exact same fasta file(s) used there
dataset = import_dataset_diann(filename = "260404_LiverTroutProteome.parquet")
dataset = import_fasta(dataset, files = "01_input_files/UP193380_110924_HaoUnivContam.fasta")

# 2) Create a template file that describes all samples. A new Excel table will be created at this path
# - note; you only have to do steps 2 and 3 once per dataset
# write_template_for_sample_metadata(dataset, "sample_metadata.xlsx")

# 3) Time to step away from R for a sec, and edit this template file in Excel or LibreOffice;
# - describe the sample group of each sample in the "group" column
# - add additional columns with any metadata that varies between samples (measurement order, gel, gel lane, batch, etc.) -->> QC figures will be auto generated
# - further documentation is available in the "instructions" tab within the Excel file

# 4) Load sample metadata from file you just edited (don't forget to save it first)
dataset = import_sample_metadata(dataset, filename = "01_input_files/metadata.xlsx")

# 5) Optionally, describe a statistical contrast; in this example we compare sample groups "WT" and "KO".
# - you should use exact same labels as "group" column in sample metadata table.
# - if you don't want to do stats, simply remove or comment this line (e.g. just look at QC report, or maybe your dataset has 1 only experimental group).
# - example for multiple contrasts; dataset = setup_contrasts(dataset, contrast_list = list( c("control", "condition_a"),  c("control", "condition_b")  ) )
# - alternatively, use the more flexible function add_contrast()  (requires MS-DAP version 1.2 or newer)
#     - check the extensive set of examples with R command: ?add_contrast
dataset = setup_contrasts(dataset, contrast_list = list(
  c("NC_day_01", "HC_day_01"),
  c("NC_day_02", "HC_day_02"),
  c("NC_day_03", "HC_day_03"),
  c("NC_day_04", "HC_day_04"),
  c("NC_day_10", "HC_day_10"),
  c("NC_day_15", "HC_day_15"),
  c("NC_day_22", "HC_day_22")
) , random_variables = c("trout_liver_batch")
)

# print an overview of all contrasts
print_contrasts(dataset)


# 6) Main function that runs the entire pipeline
# for DIA, recommended settings are defined below, selecting only peptides that were confidently detected/identified in most samples
# for DDA, 'confident detection' relies on MS/MS which may be more rare (relying on match-between-runs instead)
# so for DDA we recommend to set no or minimal requirements on 'detect' parameters; "filter_fraction_detect = 0" and "filter_min_detect = 0" (or 1 if you want at least 1 MS/MS detect per peptide per sample group)
dataset = analysis_quickstart(
  dataset,
  filter_min_detect = 3,            # each peptide must have a good confidence score in at least N samples per group
  filter_min_quant = 3,             # similarly, the number of reps where the peptide must have a quantitative value
  filter_fraction_detect = 0.5,    # each peptide must have a good confidence score in at least 50% of samples per group
  filter_fraction_quant = 0.5,     # analogous for quantitative values
  filter_min_peptide_per_prot = 1,  # minimum number of peptides per protein (after applying above filters) required for DEA. Set this to 2 to increase robustness, but note that'll discard approximately 25% of proteins in typical datasets (i.e. that many proteins are only quantified by 1 peptide)
  filter_by_contrast = TRUE,        # only relevant if dataset has 3+ groups. For DEA at each contrast, filters and normalization are applied on the subset of relevant samples within the contrast for efficiency, see further MS-DAP manuscript. Set to FALSE to disable and use traditional "global filtering" (filters are applied to all sample groups, same data table used in all statistics)
  norm_algorithm = c("vsn", "modebetween_protein"), # normalization; first vsn, then mode between on protein-level (applied sequentially so the MS-DAP mode between algorithm corrects scaling/balance between-sample-groups). "mwmb" is a good alternative for "vsn"
  dea_algorithm = c("deqms", "ebayes", "msqrobsum", "msqrob", "msempire"), # statistics; apply multiple methods in parallel/independently
  dea_qvalue_threshold = 0.05,                      # threshold for significance of adjusted p-values in figures and output tables
  dea_log2foldchange_threshold = NA,                # threshold for significance of log2 fold changes. 0 = disable, NA = automatically infer through bootstrapping
  diffdetect_min_peptides_observed = 2,             # 'differential detection' only for proteins with at least 2 peptides. The differential detection metric is a niche usecase and mostly serves to identify proteins identified near-exclusively in 1 sample group and not the other
  diffdetect_min_samples_observed = 3,              # 'differential detection' only for proteins observed in at least 3 samples per group
  diffdetect_min_fraction_observed = 0.5,           # 'differential detection' only for proteins observed in 50% of samples per group
  output_qc_report = TRUE,                          # optionally, set to FALSE to skip the QC report (not recommended for first-time use)
  output_abundance_tables = TRUE,                   # optionally, set to FALSE to skip the peptide- and protein-abundance table output files
  output_dir = "01_MSDAP_DEA",                     # output directory, here set to "msdap_results" within your working directory. Alternatively provide a full path, eg; output_dir="C:/path/to/myproject",
  output_within_timestamped_subdirectory = TRUE,
  multiprocessing_maxcores=2
)


print_dataset_summary(dataset)
