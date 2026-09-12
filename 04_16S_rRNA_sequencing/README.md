
## Purpose

This module processes 16S rRNA gene-sequencing data from midgut digesta- and mucus-associated microbiota of rainbow trout. It covers sequencing quality control, QIIME2 processing, diversity and community-structure analyses, differential abundance, microbiota–host trait associations, and preparation of microbiota matrices for multi-omics integration.

Digesta and mucus are analysed as separate microbial compartments after the shared import and initial quality-control stages.


## Amplicon library preparation and sequencing

The V3–V4 region of the bacterial 16S rRNA gene was amplified from midgut digesta and mucus gDNA using primers 341F (S-D-Bact-0341-b-S-17) and 785R (S-D-Bact-0785-a-A-21). Amplicons were generated using a two-step PCR workflow.

Library preparation, indexing, pooling, and sequencing were performed at the Genomic Platform of Bordeaux (PGTB; Univ. Bordeaux, INRAE, BIOGECO, F-33610 Cestas, France), with support from Erwan Guichoux and Prescillia Alves-Gomes. The pooled amplicon library was sequenced on an Illumina NextSeq 2000 P1 platform using 2 × 300 bp paired-end. Barcode-demultiplexed paired-end FASTQ files were used as input for the QIIME2 workflow.

PGTB is part of the INRAE Genomics Research Infrastructure and the France Génomique network (https://doi.org/10.15454/1.5572396583599417E12).


## Directory layout

- `QC/` — global sequencing quality-control reports for the original raw forward (`R1`) and reverse (`R2`) reads. This QC was performed after sequencing and before any preprocessing in QIIME2.
- `qiime2/` — QIIME2 shell workflow for preprocessing, denoising, contamination filtering, taxonomy assignment, feature filtering, and phylogenetic-tree construction. It also retains selected lightweight QIIME2 summaries and visualisations.
- `downstream_analysis/` — R and Python analyses based on the processed QIIME2 artifacts.


## Data availability and requirements

Raw sequencing reads are not stored in this Git repository. The European Nucleotide Archive (ENA) under study accession PRJEB126088 will be made publicly available upon publication.

The QIIME2 workflow requires demultiplexed paired-end reads in the format specified in `qiime2_pipeline.sh`, a QIIME2-compatible metadata file, and a QIIME2 2024.10 amplicon environment with the required plugins. The workflow uses QIIME2, DADA2, Cutadapt, RESCRIPt, FastQC, MultiQC, and SILVA 138.2 SSU NR99 resources.

Before rerunning the workflow, inspect the paths, input file names, available cores, and software environment specified in the scripts.


## QIIME2 workflow

### Main scripts

- `qiime2/qiime2_pipeline.sh` — main QIIME2 workflow.
- `qiime2/python_script_fastaQC_report.py` — helper script for per-sample FastQC reporting of exported reads during trimming checks.
- `qiime2/python_script_multiQC_report.py` — helper script that aggregates FastQC outputs into global forward- and reverse-read MultiQC reports during trimming checks.


### Processing steps

`qiime2_pipeline.sh` performs the following operations:

1. Imports demultiplexed paired-end reads and generates an initial QIIME2 demultiplexing summary.
2. Removes the V3–V4 primer sequences, Illumina adapter sequences, poly-G/poly-A tails, reads shorter than 180 bp, and low-quality 3′ bases using Cutadapt.
3. Exports selected trimmed reads for FastQC and MultiQC checks during trimming optimisation.
4. Denoises paired reads with DADA2 to generate an ASV feature table, representative ASV sequences, and denoising statistics.
5. Obtains and processes SILVA 138.2 SSU NR99 reference sequences with RESCRIPt, extracts the V3–V4 target region, and prepares the taxonomy classifier.
6. Identifies contaminants independently in digesta, mucus, and feed-control samples using the QIIME2 `decontam-identify` frequency method and DNA concentration metadata. Features with contaminant probability >0.1 are removed within each compartment.
7. Merges decontaminated compartment-specific tables and representative sequences, assigns taxonomy, retains features observed in at least two samples and assigned at phylum level, and constructs a rooted phylogenetic tree with MAFFT and FastTree.


## Downstream analyses

Run the R analyses from `downstream_analysis/`. The principal order is shown below.

### Import, normalisation, composition, diversity, and batch diagnostics

1. `01_16s_import_preprocess.R` — imports the final QIIME2 feature table, taxonomy, representative sequences, phylogenetic tree, and metadata into `microeco`; retains bacterial ASVs; creates raw and SRS-normalised datasets; calculates alpha and beta diversity; and writes RData objects and diversity tables. SRS normalisation is performed at 50,000 reads per sample for diversity analyses.
2. `01_ASVs_collapsed_tables_RA_summary.R` — produces relative-abundance summaries.
3. `02_microbiota_composition_microViz.R` — creates compositional figures with `microViz` and `phyloseq`.
4. `03_16S_diversity_analysis.R` — performs alpha-diversity analyses, beta-diversity analyses, PCoA/PERMANOVA, dbRDA, Euler-diagram analyses, and related figures. Alpha/beta diversity, ordination, and dbRDA are implemented with `microeco`; see the [microeco tutorial](https://chiliubio.github.io/microeco_tutorial/). Euler diagrams are generated with [`eulerr`](https://jolars.github.io/eulerr/articles/introduction.html) from ASV presence/absence sets, showing unique and shared ASVs between sample types or diets and reporting both feature-number (`Numratio`) and sequence-abundance (`Seqratio`) proportions.
5. `04_PLSDAbatch_exploration.R` — performs an exploratory, compartment-specific assessment of technical batch structure in digesta and mucus microbiota. The workflow is adapted from the [PLSDAbatch brief vignette](https://bioconductor.posit.co/packages/3.20/bioc/vignettes/PLSDAbatch/inst/doc/brief_vignette.html), specifically its preprocessing and batch-effect detection steps. Raw ASV counts are pre-filtered with `PLSDAbatch::PreFL()` and CLR-transformed with an offset of 1. Four technical variables are evaluated independently: DNA-extraction date; grinding method, representing the device used for mechanical cell lysis; PCR1 plate, corresponding to the first PCR amplification of the V3–V4 16S rRNA gene region; and PGTB plate, corresponding to the Illumina-platform sequencing plate. For each compartment and technical factor, the script generates PCA density panels, boxplots and density plots for the ASV with the largest absolute PC1 loading, a linear model for that ASV, and variation partitioning with `vegan::varpart()`. Linear models and variation partitioning account for the full `Diet × Day` biological structure through a combined diet-by-day grouping factor. This is an exploratory batch-QC workflow only: no `PLSDA_batch()` correction is performed, no samples are removed, and no corrected microbiota matrix is used in downstream analyses.
7. `05_export_files_for_BIRDMAn_DA_analysis.R` — exports the unnormalised ASV count table, taxonomy table, and phylogenetic tree required by the BIRDMAn branches.

### Primary ASV-level differential abundance: BIRDMAn

The primary differential-abundance analysis is performed separately for digesta and mucus:

- `DA_analysis_BIRDMAn/digesta/`
- `DA_analysis_BIRDMAn/mucus/`

Each branch contains:

1. `01_BIRDMAn_*.ipynb` — prepares the filtered compartment-specific BIOM table, fits the custom Bayesian Stan model through BIRDMAn, and writes posterior inference and diagnostics.
2. `02_export_BIRDMAn_results_for_R.ipynb` — derives HC–NC posterior contrasts for each sampling day and exports posterior summaries and contrast-specific diagnostics.
3. `03_DA_*_heatmap.R` — creates ASV-level differential-abundance heatmaps.
4. `04_spearman_correlation_*.ipynb` — calculates ASV–host trait correlations from CLR-transformed count profiles.
5. `05_spearman-corr_*.R` — creates correlation heatmaps and related R outputs.

#### BIRDMAn model

BIRDMAn was used for ASV-level Bayesian differential-abundance analysis. Separate digesta and mucus notebooks define compartment-specific custom Stan models and fit them through `birdman.TableModel`.

Unlike the `SingleFeatureModel approach`, which fits an independent model for each ASV, `TableModel` fits the complete ASV count table in a single joint Stan model. Thus, all ASVs are analysed within one model invocation, while retaining ASV-specific parameters (including abundance effects and dispersion).


##### Experimental design and compositional parameterisation

The regression design is `Diet * Day`, with NC and `day_01` as reference levels. The model is parameterised in additive log-ratio (ALR) space using the same selected reference ASV in both digesta and mucus.

Posterior coefficient arrays are subsequently centred to CLR coordinates before HC–NC contrasts are calculated separately for each sampling day.

##### Feature filtering

Filtering is compartment specific:

- Digesta: ASVs present in at least 20% of retained samples.
- Mucus: ASVs present in at least 20% of retained samples.

##### Model comparison, diagnostics, and reporting

The fitted `Diet * Day` model is compared with an intercept-only model using leave-one-out cross-validation.

The current notebooks use four MCMC chains, 500 warm-up iterations, and 500 retained draws per chain. Convergence and sampling quality are evaluated using R-hat and effective sample size. Exported contrast diagnostics require R-hat <1.05 and ESS-tail >400.

ASVs displayed in differential-abundance heatmaps additionally require probability of direction >0.95 and a 95% highest-density interval excluding zero.

##### Documentation

The implementation follows the BIRDMAn framework. Documentation covering custom Stan models, `TableModel`, posterior diagnostics, and negative-binomial models is available at [BIRDMAn documentation](https://birdman.readthedocs.io/en/latest/).

### Taxonomic-level sensitivity analysis: microViz/corncob

`06_DA_analysis_digesta_microViz.R` and `07_DA_analysis_mucus_microViz.R` provide a frequentist sensitivity analysis after taxonomic collapsing.

Each script imports the BIRDMAn-filtered BIOM table with sample metadata, repairs and prefixes taxonomic labels, and fits `corncob::bbdml` beta-binomial regression models through `microViz`. Models are run at phylum, class, order, family, and genus levels using `Diet` as the explanatory variable, with NC as the reference level. Thus, this analysis estimates an overall HC–NC association pooled across sampling days; it does not fit a `Diet * Day` interaction.

The beta-binomial model accounts for variable library sizes by modelling taxon counts relative to the total count per sample. P-values are adjusted by the Benjamini–Hochberg method separately within taxonomic rank.

These models should be interpreted as a taxonomic-level sensitivity analysis of the primary ASV-level BIRDMAn results, not as a replacement for the BIRDMAn analysis.

The taxonomic sensitivity-analysis workflow follows the microViz guidance for statistical modelling of individual taxa and taxonomic association trees: [microViz taxon-modelling vignette](https://david-barnett.github.io/microViz/articles/web-only/modelling-taxa.html).

### Microbiota–host trait associations

The Spearman-correlation notebooks reuse the BIRDMAn-filtered ASV BIOM table and the BIRDMAn model metadata for their respective compartment. Counts are zero-replaced with multiplicative replacement, CLR-transformed, and correlated with host traits using Spearman rank correlation. Multiple testing is controlled with the Benjamini–Hochberg procedure.

### Preparation for MOFA integration

1. `08_ASVs_tables_tsv_for_mofa.R` — creates ASV-level digesta and mucus matrices for multi-omics integration, writes ASV-to-taxonomy mappings, and applies CLR transformation.
2. `09_process_ASV_tables_for_MOFA.R` — harmonises sample identifiers and writes sample-by-feature microbiota RDS matrices for MOFA.

The MOFA export uses the same compartment-specific, prevalence-filtered ASV tables as the differential-abundance and Spearman-correlation workflows. These tables are subsequently CLR-transformed and harmonised for sample identifiers before multi-omics integration.

## Data lineage and filtering

| Analysis | Starting count table | Additional filtering/transformation |
|---|---|---|
| BIRDMAn DA | Compartment-specific ASV count table exported by `05_export_files_for_BIRDMAn_DA_analysis.R` | Digesta: prevalence ≥20%; mucus: prevalence ≥10%; raw counts retained for NB2 modelling |
| microViz/corncob DA | BIRDMAn-filtered BIOM table with sample metadata | Taxonomy repair and rank-prefixing; taxonomic aggregation; beta-binomial regression |
| Spearman correlations | BIRDMAn-filtered BIOM table with taxonomy plus BIRDMAn model metadata | CLR transformation |
| MOFA | Same compartment-specific prevalence-filtered ASV table used for DA and Spearman analyses | CLR transformation and sample-ID harmonisation |

Accordingly, the BIRDMAn and Spearman workflows use the same retained compartment-specific ASV set. The microViz workflow starts from that same BIRDMAn-filtered ASV table. The MOFA workflow is separate.

## Main outputs

- `QC/` — global raw-read R1/R2 quality-control reports.
- `qiime2/` summary folders — selected QIIME2 visualisations and lightweight tables documenting import, trimming, denoising, decontamination, and taxonomy assignment.
- `downstream_analysis/RData/` — `microeco` objects for raw, compartment-specific, and SRS-normalised data.
- `downstream_analysis/01_diversity_tables/` — alpha- and beta-diversity tables.
- `downstream_analysis/06_DA_analysis_digesta_microViz/` and `07_DA_analysis_mucus_microViz/` — taxonomic beta-binomial model objects, result tables, plots, and session information.
- `downstream_analysis/DA_analysis_BIRDMAn/*/output_files/` — filtered BIOM tables, model metadata, dispersion-prior tables, generated Stan code, posterior inference, LOO comparisons, and diagnostics.
- `downstream_analysis/DA_analysis_BIRDMAn/*/export_for_R/` — posterior contrast summaries, convergence tables, ASV-level DA heatmaps, and microbiota–host trait correlation outputs.
- `downstream_analysis/09_process_ASV_tables_for_MOFA/` — digesta and mucus sample-by-feature matrices and associated mappings for multi-omics integration.

## Main software

QIIME2, RESCRIPt, DADA2, FastQC, MultiQC, R (`microeco`, `microViz`, `phyloseq`, `vegan`, `PLSDAbatch`), and Python (`BIRDMAn`, `biom-format`, `arviz`, `xarray`, `scikit-bio`, `pandas`, `numpy`).
