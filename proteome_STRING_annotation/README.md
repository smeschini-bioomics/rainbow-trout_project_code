# STRING Annotation Coverage and Statistics

## Purpose

This workflow evaluates the coverage and structure of STRING functional annotations for the complete *Oncorhynchus mykiss* UniProt reference proteome.

It is a descriptive quality-control and documentation workflow.

The analysis:

- matches protein identifiers between the UniProt reference FASTA and the STRING enrichment-term file;
- quantifies overall annotation coverage of the reference proteome;
- quantifies annotation coverage separately for each STRING category;
- summarises the number of distinct terms assigned to each annotated protein within each category;
- summarises the number of proteins assigned to each term within each category

## Generation of the custom STRING proteome

The complete *Oncorhynchus mykiss* UniProt reference proteome was uploaded to STRING using the **“Annotate your proteome”** / **“Adding a new species to STRING”** workflow.

The uploaded FASTA provided the protein sequences used by STRING to create a custom rainbow trout proteome, infer functional and physical interaction networks, and assign functional annotations, including Gene Ontology and KEGG terms. STRING assigned this uploaded proteome the custom organism identifier:

```text
STRG0A55HWH
```

## Input files

Place the following files in this directory before running the script.

| File | Description |
|---|---|
| `24_10_02_uniprot-ref_proteome_trout_UP000193380.fasta` | Complete *O. mykiss* UniProt reference proteome used as the annotation universe. |
| `110079946.protein.enrichment.terms.v12.0.txt` | STRING v12.0 protein-to-term annotation file for rainbow trout. |

The script expects standard UniProt FASTA headers, for example:

```text
>sp|A0A123ABC4|PROTEIN_NAME ...
>tr|A0A123ABC4|PROTEIN_NAME ...
```

Protein identifiers are extracted from the second pipe-delimited field of the UniProt header. STRING identifiers are converted to protein accessions by removing the organism-specific prefix before the first period.

## Script

Run:

```text
annotation_statistics_reference_proteome.R
```

The script uses the following R packages:

- `readr`
- `dplyr`
- `stringr`
- `tidyr`
- `ggplot2`
- `patchwork`

## Workflow

### 1. Import and clean STRING annotations

The STRING enrichment-term file is imported with the following fields:

```text
string_protein_id
category
term
description
```

The script removes:

- empty protein identifiers, categories, and terms;
- contaminants, reverse sequences, and decoy identifiers beginning with `Cont_`, `REV_`, or `DECOY_`;
- duplicated protein-category-term combinations.

### 2. Import the UniProt reference proteome

All FASTA headers are read from the complete UniProt reference proteome.

The script records:

- total FASTA header entries;
- number of unique parsed protein identifiers;
- number of duplicated FASTA identifiers, if present.

The unique protein identifiers in this FASTA define the reference annotation universe.

### 3. Restrict STRING annotations to the reference proteome

The workflow then calculates:

- reference proteins with at least one STRING annotation;
- reference proteins without STRING annotation;
- the percentage of the reference proteome with at least one annotation;

### 4. Calculate category-specific annotation statistics

All annotation summaries are calculated separately within each STRING annotation category.

For each protein-category combination, the script calculates:

```text
n_terms
```

which represents the number of distinct terms assigned to that protein within that annotation category.

For each term-category combination, the script calculates:

```text
n_proteins
```

which represents the number of distinct reference proteins assigned to that term within that annotation category.


### 5. Generate tables and plots

The workflow exports global, identifier-matching, category-specific, protein-level, and term-level statistics.

It also generates plots describing:

1. overall reference-proteome coverage;
2. reference-proteome coverage by annotation category;
3. number of unique terms by annotation category;
4. distribution of distinct terms per annotated protein within each category;
5. a combined summary figure containing all four panels.

## Output directory

All outputs are written to:

```text
annotation_statistics_reference_proteome/
```

### Tables

| File | Description |
|---|---|
| `global_annotation_statistics.tsv` | Overall FASTA, STRING, and annotation-coverage statistics. |
| `identifier_overlap_diagnostic.tsv` | Identifier-matching diagnostic comparing the uploaded reference FASTA with the STRING term file; reports shared protein IDs and any unmatched IDs in either file. |
| `annotation_statistics_by_category.tsv` | Coverage, term counts, and protein/term summary statistics for each STRING category. |
| `protein_annotation_statistics_by_category.tsv` | Number and identity of terms assigned to each annotated protein within each category. |
| `term_annotation_statistics_by_category.tsv` | Number and identity of proteins assigned to each term within each category. |
| `unannotated_reference_proteins.tsv` | Reference FASTA proteins without a retained STRING annotation. |
| `string_proteins_not_in_reference_fasta.tsv` | STRING proteins not found in the supplied reference FASTA. |

### Figures

Each figure is exported as both PDF and TIFF at 1200 dpi with LZW compression.

| File stub | Content |
|---|---|
| `01_reference_proteome_coverage` | Overall reference-proteome annotation coverage. |
| `02_coverage_by_category` | Percentage of reference proteins with at least one term per category. |
| `03_unique_terms_by_category` | Number of represented terms per category. |
| `04_terms_per_protein_by_category` | Distribution of distinct terms per annotated protein within each category. |
| `05_combined_annotation_statistics` | Combined four-panel annotation-statistics figure. |
