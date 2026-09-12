# Liver proteome

## Purpose

This module analyses liver data-independent acquisition (DIA) LC–MS/MS proteomics data quantified with DIA-NN. It performs protein-level differential-abundance analysis, prepares the liver sample-by-protein matrix, calculates GSVA pathway scores, assesses technical batch structure, and evaluates multivariate proteome profiles.

The computational workflow is the same as for the midgut proteome. See the [Midgut Proteome README](../07_Midgut_Proteome/README.md) for the full MS-DAP, GSVA, batch-correction, and PCA/PERMANOVA workflow description. This README documents only liver-specific aspects.

## Sample fractionation and LC–MS/MS analysis

Liver samples were processed and analysed by [ProGenTomics](https://www.progentomics.ugent.be/), Ghent, Belgium.

The liver proteome and liver hPTM datasets were generated from sequential fractions of the same extraction workflow. Following acid addition for histone extraction and centrifugation, the acid-soluble fraction was used for liver hPTM analysis. The remaining acid-insoluble pellet was retained for DIA proteomics.

Thus, this dataset represents proteins recovered from the non-histone, acid-insoluble fraction rather than from the acid extract used for liver hPTM analysis. The proteome and hPTM measurements are complementary molecular readouts derived from sequential fractions of the same biological sample.

## Liver-specific inputs

Run scripts from `08_Liver_Proteome/`.

The module uses the same DIA-NN/MS-DAP input structure as the midgut workflow:

- DIA-NN quantitative report (`260404_LiverTroutProteome.parquet`);
- FASTA database used for DIA-NN searching;
- MS-DAP sample metadata;
- `STPN2309_metadata_all.tsv`;
- STRING v12 rainbow trout enrichment-term annotation file.

The large DIA-NN `260404_LiverTroutProteome.parquet` file is not committed. It will be deposited in an external public repository before publication; the accession number and download link will be added in the public repository release.

The STRING annotation file is also external and is not committed. Download it from the [STRING v12 rainbow trout organism page](https://version-12-0.string-db.org/organism/STRG0A55HWH) and place it in this module directory.

## Liver-specific metadata and technical batch

The liver workflow maps proteomics sample labels to project sample identifiers using the liver proteomics sample-ID field in `STPN2309_metadata_all.tsv`.

The technical batch variable is:

```text
trout_liver_batch
```

It is used for liver-specific batch assessment and, where applicable, for batch adjustment while protecting the biological design:

```r
~ diet * day_categ
```

As in the midgut workflow, protein-level MS-DAP differential-abundance analysis remains based on the MS-DAP contrast-specific processing. Batch-corrected matrices are used only for profile-level visualisation and multivariate analyses, not as replacement input for MS-DAP DEA.

## Liver-specific analysis outputs

The same analysis sequence is used as for the midgut proteome:

1. MS-DAP QC, filtering, normalisation, and day-specific HC-versus-NC DEA.
2. Cross-day DEA summaries and count tables.
3. Processed liver sample-by-protein matrix.
4. STRING-based GSVA and Biological Process heatmaps.
5. Liver batch exploration and, where configured, limma batch correction.
6. PCA, PERMANOVA, dispersion testing, and PC-score comparisons.

The central liver matrices are saved under the module `rds/` and batch-correction output directories, using liver-specific file names such as:

```text
liver_proteome_filtered.rds
liver_proteome_filtered_batch_corrected_limma.rds
```

## Interpretation

Day-specific MS-DAP contrasts estimate HC-versus-NC protein abundance differences independently at days 1, 2, 3, 4, 10, 15, and 22.

GSVA scores represent relative pathway activity inferred from the processed liver protein matrix. PCA/PERMANOVA analyses describe global liver proteome structure after the configured handling of technical batch variation.

## Main R packages

The liver workflow uses the same core packages as the midgut proteome module, including `msdap`, `data.table`, `dplyr`, `tidyr`, `GSVA`, `limma`, `sva`, `vegan`, `ComplexHeatmap`, `ggplot2`, and associated data-handling and plotting packages.
