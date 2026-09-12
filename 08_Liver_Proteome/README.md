# Liver proteome

## Purpose

This module analyses liver data-independent acquisition (DIA) LC–MS/MS proteomics data quantified with DIA-NN. It performs protein-level differential-abundance analysis, prepares the liver sample-by-protein matrix, calculates GSVA pathway scores, assesses technical batch structure, and evaluates multivariate proteome profiles.

The computational workflow is the same as for the midgut proteome. See the [Midgut Proteome README](../07_Midgut_Proteome/README.md) for the full MS-DAP, GSVA, batch-correction, and PCA/PERMANOVA workflow description.

## Sample fractionation and LC–MS/MS analysis

Liver samples were processed and analysed by [ProGenTomics](https://www.progentomics.ugent.be/), Ghent, Belgium.

The liver proteome and liver hPTM datasets were generated from sequential fractions of the same extraction workflow. Following acid addition for histone extraction and centrifugation, the acid-soluble fraction was used for liver hPTM analysis. The remaining acid-insoluble pellet was retained for DIA proteomics.

Thus, this dataset represents proteins recovered from the non-histone, acid-insoluble fraction rather than from the acid extract used for liver hPTM analysis. The proteome and hPTM measurements are complementary molecular readouts derived from sequential fractions of the same biological sample.


## Main R packages

The liver workflow uses the same core packages as the midgut proteome module, including `msdap`, `data.table`, `dplyr`, `tidyr`, `GSVA`, `limma`, `sva`, `vegan`, `ComplexHeatmap`, `ggplot2`, and associated data-handling and plotting packages.
