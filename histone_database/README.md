# Custom Histone FASTA Database

## Purpose

This directory contains in-house utility scripts used to assemble a practical custom rainbow trout histone FASTA background for Mascot-compatible in silico spectral-library preparation in the liver and midgut hPTM analyses.

The database combines manually collected candidate H1, H2A, H2B, H3, and H4-related sequences from four rainbow trout proteome resources and MS HistoneDB 2.0.

These scripts are intended to facilitate sequence collection, source tracking, exact-sequence de-duplication, and FASTA-header formatting. They do **not** provide a fully manually curated, experimentally validated, or exhaustive reference catalogue of rainbow trout histone variants. As rainbow trout is a non-model species with incomplete and sometimes inconsistent histone database annotation, the input sequence set and selection rules can be adapted according to the specific analytical objective.

## Directory structure

```text
histone_database/
├── 01_FASTA_MS_HistoneDB/
│   ├── histone_variants_H1_MS_HistoneDB.fasta
│   ├── histone_variants_H2A_MS_HistoneDB.fasta
│   ├── histone_variants_H2B_MS_HistoneDB.fasta
│   ├── histone_variants_H3_MS_HistoneDB.fasta
│   ├── histone_variants_H4_MS_HistoneDB.fasta
│   └── 24_10_08_histone_seq_MS_HistoneDB.fasta
│
├── 02_input_files/
│   ├── 24_10_02_histone_prot_seq_annotation_from_UP000193380.txt
│   ├── 24_10_02_histone_prot_seq_annotation_from_UP000694395.txt
│   ├── 24_10_02_histone_prot_seq_annotation_from_RefSEQ_NCBI_USDA_OmykA_1_1.txt
│   └── 24_10_04_histone_prot_seq_annotation_from_ENSEMBL.txt
|
├── histone_alignement/
│   ├── README.md
│   ├── pymsaviz.ipynb
│   ├── H3_alignement.fas
│   ├── H4_alignement.fas
│   ├── 20260420_TroutGut_peptides_ions_raw_ab_rank3_msqrob-input.csv
│   ├── 260418_TroutLiver_H3H4_raw_ab_msqrob_input.csv
│   └── <dataset-specific output directories>
│
├── output_files/
|
├── 01_merging_fasta_files_from_MS_HistoneDB_2_0.py
├── 01_merging_fasta_files_from_MS_HistoneDB_2_0.py
├── 02_non_redundant_histone_database_from_fasta_files.py
└── 03_simplify_histones_fasta_header_for_input_Mascot.py
```
## Histone peptide and PTM alignment visualisation

The [`histone_alignement/`](histone_alignement/) subdirectory contains a separate downstream visualisation workflow for histone H3 and H4 peptide coverage and PTM-site annotation.

It uses the peptide-level input tables prepared for the midgut and liver msqrobPTM workflows, together with cross-species H3 and H4 multiple-sequence alignments. The workflow visualises experimentally detected peptide sequences by colour-coding their aligned residues and overlays previously annotated PTM sites as symbols.

This alignment workflow does not modify the Mascot search database, infer histone variants, or independently assign PTM coordinates. Peptide conservation and canonical human-reference PTM-site annotation were manually verified before generating these figures.

See the [dedicated alignment-workflow README](histone_alignement/README.md) for its inputs, prior annotation procedure, alignment generation, notebook workflow, outputs, and reproducibility information.


## Input sequence resources

Candidate histone sequences were manually selected before this workflow from the following rainbow trout protein resources:

| Source | Input file |
|---|---|
| UniProt reference proteome `UP000193380` | `02_input_files/24_10_02_histone_prot_seq_annotation_from_UP000193380.txt` |
| UniProt proteome `UP000694395` | `02_input_files/24_10_02_histone_prot_seq_annotation_from_UP000694395.txt` |
| NCBI RefSeq, USDA OmykA_1.1 | `02_input_files/24_10_02_histone_prot_seq_annotation_from_RefSEQ_NCBI_USDA_OmykA_1_1.txt` |
| Ensembl, USDA OmykA_1.1 | `02_input_files/24_10_04_histone_prot_seq_annotation_from_ENSEMBL.txt` |
| MS HistoneDB 2.0 | `01_FASTA_MS_HistoneDB/24_10_08_histone_seq_MS_HistoneDB.fasta` |

The files in `02_input_files/` are FASTA-formatted despite using a `.txt` extension.

Manual selection of histone candidates from the complete rainbow trout proteome resources occurred before this workflow. This extraction step is not automated by the scripts in this directory.

## Workflow

### 1. Merge MS HistoneDB classes

`01_merging_fasta_files_from_MS_HistoneDB_2_0.py` merges the five class-specific MS HistoneDB FASTA files:

- `histone_variants_H1_MS_HistoneDB.fasta`
- `histone_variants_H2B_MS_HistoneDB.fasta`
- `histone_variants_H2A_MS_HistoneDB.fasta`
- `histone_variants_H3_MS_HistoneDB.fasta`
- `histone_variants_H4_MS_HistoneDB.fasta`

The historical merge order is retained: H1, H2B, H2A, H3, then H4.

The script does not modify headers or remove duplicate records.

**Output**

```text
01_FASTA_MS_HistoneDB/24_10_08_histone_seq_MS_HistoneDB.fasta
```

### 2. Collapse exact duplicate sequences across sources

`02_non_redundant_histone_database_from_fasta_files.py` combines the four manually selected rainbow trout histone sequence subsets with the merged MS HistoneDB FASTA.

For each source, the script extracts a concise source-specific identifier and protein description. Complete amino-acid sequences are then used as the grouping key:

- entries with 100% identical amino-acid sequences are represented by one FASTA record;
- source identifiers and descriptions from identical sequences are joined using ` | `;
- sequences differing by one or more amino acids are retained as separate entries.

The script also generates class-specific FASTA files and a sequence statistics file.

**Outputs**

```text
output_files/24_10_09_new_non_redundant_histone_DB_sequences.fasta
output_files/H1_sequences.fasta
output_files/H2A_sequences.fasta
output_files/H2B_sequences.fasta
output_files/H3_sequences.fasta
output_files/H4_sequences.fasta
output_files/24_10_09_sequence_statistics.txt
```

### 3. Simplify FASTA headers for Mascot input

`03_simplify_histones_fasta_header_for_input_Mascot.py` evaluates the merged source annotations retained in each non-redundant FASTA header and selects one representative annotation.

The scoring heuristic prioritises descriptions containing:

- histone-variant notation;
- `variant` or `isoform`;
- `centromeric`, `family`, or `like`.

Description length is used as a secondary ranking criterion.

This step changes FASTA headers only. It does not alter amino-acid sequences and does not perform further sequence de-duplication.

**Output**

```text
output_files/24_10_09_histones_sequences_header_simplified_mascot_input.fasta
```


## Database scope and peptide-level interpretation

The custom FASTA was assembled to provide Mascot with a broad set of plausible histone-related amino-acid sequences for peptide-spectrum matching and in silico spectral-library preparation. Its purpose is to maximise the sequence space considered during matching, including closely related sequences that differ by one or a few residues, because these differences can substantially change the theoretical peptidoforms and fragment-ion patterns generated when combinatorial histone PTMs are considered.

The database was therefore not constructed to establish definitive protein, canonical histone, or histone-variant identity. The source annotation attached to a matching sequence is not used as the analytical basis of the study.

Because the DDA bottom-up hPTM workflow cannot resolve histone-variant-specific abundance, no variant-level protein inference was attempted. Many histone tail peptides are shared across closely related variants, preventing reliable assignment of their signal to an individual variant.

Thus, the Mascot protein label is only the representative database annotation associated with a peptide-spectrum match. At post-acquisition review, labels were corrected only when needed to place the appropriate histone family set.

Importantly, as with any database-search approach, a peptide-spectrum match may be explained equally well or even better by a biologically relevant sequence that is not represented in the search database. In addition, the set of considered PTM combinations is necessarily restricted because the theoretical search space increases rapidly as the number of potentially modified sites and co-occurring modifications increases.

## File flow

```text
MS HistoneDB H1/H2B/H2A/H3/H4 FASTA files
        |
        |  01_merging_fasta_files_from_MS_HistoneDB_2_0.py
        v
24_10_08_histone_seq_MS_HistoneDB.fasta
        |
        +---------------------------------------------+
        |                                             |
        v                                             |
Candidate histone subsets from UP000193380,           |
UP000694395, RefSeq, and Ensembl                      |
        |                                             |
        |  02_non_redundant_histone_database_from_fasta_files.py
        v
24_10_09_new_non_redundant_histone_DB_sequences.fasta
        |
        |  03_simplify_histones_fasta_header_for_input_Mascot.py
        v
24_10_09_histones_sequences_header_simplified_mascot_input.fasta

Archived historical FASTA used for hPTM analysis:
24_10_09_histones_sequences_header_simplified_mascot_input_older.fasta
```

## Running the workflow

Run the scripts from the `histone_database/` directory:

```powershell
py .\01_merging_fasta_files_from_MS_HistoneDB_2_0.py
py .\02_non_redundant_histone_database_from_fasta_files.py
py .\03_simplify_histones_fasta_header_for_input_Mascot.py
```

## Requirements

- Python 3
- Biopython

```powershell
py -m pip install biopython
```

## Availability

The complete source proteome FASTA databases and large derived FASTA resources are not stored in the Git repository because of their size and external-resource provenance.

The exact archived Mascot-input FASTA used for hPTM identification, together with the associated large FASTA resources, is available in the associated Zenodo archive: DOI `XXXXX`.
