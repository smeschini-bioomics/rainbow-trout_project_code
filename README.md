# Rainbow Trout Multi-Omics Analysis

Reproducible analysis code for comparison of two isolipidic and isoenergetic experimental diets differing in carbohydrate-to-protein ratio diet: NC (0% digestible carbohydrate / 60% crude protein) and HC (~30% gelatinized starch / 40% crude protein) in rainbow trout (*Oncorhynchus mykiss*). The repository contains workflows for host physiological and plasmatic parameters, targeted liver mRNA, global DNA methylation, 16S rRNA gene sequencing, histone post-translational modifications (hPTMs), bottom-up proteomics, and multi-omics integration.

## Study design

The study compares two diets, **NC** and **HC**, across seven post-prandial sampling days: **day 1, 2, 3, 4, 10, 15, and 22**. Samples are from independent fish; the experiment is not longitudinal.

## Repository structure

| Directory | Content |
|---|---|
| `01_Host_Physiological_Traits/` | Body weight, plasma glucose, hepatosomatic index, and L-/D-lactate analyses; trait correlations, coefficient-of-variation summaries, and sample-size exports. |
| `02_Liver_mRNA/` | Targeted liver mRNA relative-concentration analysis. |
| `03_Liver_Global_DNA_methylation/` | Preparation and Bayesian Dirichlet modelling of global cytosine relative composition in liver. |
| `04_16S_rRNA_sequencing/` | QIIME2 preprocessing, microbiota composition/diversity, differential-abundance analyses, and microbiota matrices for MOFA. |
| `05_Midgut_hPTMs/` | Midgut histone PTM preprocessing, Bayesian modelling, batch assessment/correction, multivariate analyses, and correlation structure. |
| `06_Liver_hPTMs/` | Liver histone PTM preprocessing, Bayesian modelling, batch assessment/correction, multivariate analyses, and correlation structure. |
| `07_Midgut_Proteome/` | Midgut DIA proteomics: MS-DAP differential expression, matrix processing, GSVA, batch assessment, and multivariate analyses. |
| `08_Liver_Proteome/` | Liver DIA proteomics: MS-DAP differential expression, matrix processing, GSVA, batch assessment, and multivariate analyses. |
| `09_Multi_omics/` | MOFA data assembly, quality control, model training/comparison, factor interpretation, enrichment, and MultiGroupPower analyses. |
| `histone_database/` | In-house workflow for assembling and documenting the custom rainbow trout histone FASTA background used for in-silico hPTM spectral library preparation. It also contains a downstream H3/H4 peptide-coverage and PTM-site alignment-visualisation workflow. |
| `proteome_STRING_annotation/` | Assessment of STRING functional annotation coverage against the UniProt rainbow trout reference proteome. Exports FASTA–STRING identifier-matching diagnostics, category-wise coverage, category-specific terms-per-protein statistics, term-level annotation tables, and summary figures. |

Each top-level module contains a dedicated `README.md` with its required files, script order, outputs, and methodological notes.

## Running the code

The scripts use relative paths and are intended to be run from the relevant analysis directory, not from an arbitrary working directory.

Run numbered scripts in ascending order unless a module README specifies otherwise. Scripts create their own result directories. Some scripts contain project-specific filenames, timestamped MS-DAP output paths, or local Python paths; update these settings only when reproducing the workflow from a fresh clone.

Raw inputs and some large intermediate objects are intentionally not stored in Git. Public data-accession links and any release-specific processed-data archive will be added before publication.

## Software

Primary analyses use R. Core packages include `glmmTMB`, `brms`, `QFeatures`, `msqrob2`, `limma`, `MS-DAP`, `GSVA`, `MOFA2`, `MultiPower`, `microeco`, `microViz`, and `vegan`. The 16S workflow additionally requires QIIME2, FastQC/MultiQC, and Python notebooks for BIRDMAn-based differential abundance.

Most major scripts write a `sessionInfo.txt` file into their output directory. These files should be retained with the corresponding released results.

## Data availability

Raw 16S rRNA sequencing reads and raw mass-spectrometry files are not distributed through GitHub.

- **16S rRNA gene sequencing**: The European Nucleotide Archive (ENA) under study accession PRJEB126088 and will be made publicly available upon publication.
- **Proteomics**: the ProteomeXchange/PRIDE accession and public download page will be added upon publication.
- **Processed matrices, metadata, and publication-ready tables**: •	Any additional information required to reanalyze the data reported in the associated publication is available from the lead contact upon request.

Until publication, this repository remains private and the data are available only to authorised collaborators.


## External reference proteome, custom histone database, and annotation resources

Large external FASTA resources, STRING annotation files, and the historical custom histone FASTA used for hPTM identification are not committed to this repository.

### UniProt reference proteome

The following FASTA database was used for the midgut and liver proteome workflows and for STRING annotation-coverage assessment:

- *Oncorhynchus mykiss* UniProt reference proteome FASTA database  
  `UP000193380-8022`  
  Genome assembly: `GCA_900005705.1`  
  46,436 protein sequences  
  Downloaded: 2 October 2024  

Expected file:

- `24_10_02_uniprot-ref_proteome_trout_UP000193380.fasta`

The `proteome_STRING_annotation/` workflow uses this complete UniProt reference FASTA as its annotation universe. It evaluates FASTA-to-STRING identifier overlap, overall and category-specific annotation coverage, and category-specific protein-to-term annotation statistics.

### Custom histone FASTA database

The `histone_database/` directory contains in-house scripts used to collect, merge, de-duplicate, and format candidate rainbow trout histone sequences for Mascot-compatible in silico spectral-library preparation in the liver and midgut hPTM workflows.

The scripts combine manually collected candidate H1, H2A, H2B, H3, and H4-related sequences from:

- UniProt reference proteome `UP000193380`;
- UniProt proteome `UP000694395`;
- NCBI RefSeq, USDA OmykA_1.1;
- Ensembl, USDA OmykA_1.1;
- MS HistoneDB 2.0.

The database was designed as a broad sequence background for peptide-spectrum matching, not as a fully validated catalogue of rainbow trout histone variants.

### Histone peptide and PTM alignment visualisation

The [`histone_database/histone_alignement/`](histone_database/histone_alignement/) subdirectory contains a separate downstream visualisation workflow for histone H3 and H4.

It uses the peptide-level tables prepared for the liver and midgut msqrobPTM analyses together with cross-species H3 and H4 multiple-sequence alignments. Rainbow trout sequences are displayed alongside homologues from selected model species to visualise conserved regions, experimental peptide sequence coverage, and detected PTM sites.

Detected peptide sequences are colour-coded across the alignment, whereas pre-annotated PTM sites are displayed as symbols at their corresponding residues. PTM coordinates are already expressed relative to canonical human H3 or H4 reference sequences following prior manual peptide-conservation and site-annotation checks. This visualisation workflow does not modify the Mascot database, infer histone variants, or independently assign PTM positions.

See the [dedicated alignment-workflow README](histone_database/histone_alignement/README.md) for the alignment generation procedure, prior annotation steps, notebook workflow, outputs, and reproducibility details.

### STRING functional annotations

Functional proteome and multi-omics enrichment analyses use STRING v12.0 annotations for rainbow trout:

- [STRING v12.0: *Oncorhynchus mykiss* (rainbow trout), organism `STRG0A55HWH`](https://version-12-0.string-db.org/organism/STRG0A55HWH)

Download the following files from STRING and place them in the relevant analysis module:

- `STRG0A55HWH.protein.info.v12.0.txt`
- `110079946.protein.enrichment.terms.v12.0.txt`

Repository DOI, and licence information will be added with the public release accompanying the manuscript.
