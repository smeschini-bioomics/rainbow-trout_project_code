# Histone H3/H4 peptide and PTM alignment visualisation

## Purpose

This archive contains the input tables, protein multiple-sequence alignments, Python notebook, and generated figures used to visualise experimentally detected rainbow trout (*Oncorhynchus mykiss*) histone peptides and their pre-annotated post-translational modification (PTM) sites in cross-species H3 and H4 alignments.

The workflow has three complementary objectives:

1. **Assess histone sequence conservation.**  
   Rainbow trout H3 and H4 sequences are aligned with homologous sequences from commonly used model species, including human (*Homo sapiens*), mouse (*Mus musculus*), *Arabidopsis thaliana*, *Caenorhabditis elegans*, *Xenopus laevis*, *Drosophila melanogaster*, and *Saccharomyces cerevisiae*. This enables visual inspection of conserved and divergent histone regions.

2. **Visualise experimentally detected peptides and PTM sites.**  
   The workflow searches the amino-acid sequences of H3 and H4 peptides identified in the liver and midgut experiments against each ungapped alignment sequence. Matched peptide residues are displayed using a colour code, allowing direct visual assessment of experimental sequence coverage. Previously annotated PTM sites are overlaid as symbols at their corresponding aligned residues.

3. **Display PTM sites using a standardised nomenclature.**  
   PTM coordinates in the input tables are already assigned relative to canonical human H3 or H4 reference sequences. Human histone numbering is widely used in the literature and in commercial histone PTM antibodies. Displaying sites using this convention supports comparison of identified rainbow trout PTM sites with previously reported histone marks across species.

The workflow is designed for four datasets:

| Dataset tag | Tissue | Histone | Input table | Alignment FASTA | Human reference sequence | Reference description |
|---|---|---|---|---|---|---|
| `midgut_H3` | Midgut | H3 | `20260420_TroutGut_peptides_ions_raw_ab_rank3_msqrob-input.csv` | `H3_alignement.fas` | <code>sp&#124;P68431&#124;H31</code> | Human histone H3.1 |
| `midgut_H4` | Midgut | H4 | `20260420_TroutGut_peptides_ions_raw_ab_rank3_msqrob-input.csv` | `H4_alignement.fas` | <code>sp&#124;P62805&#124;H4</code> | Human histone H4 |
| `liver_H3` | Liver | H3 | `260418_TroutLiver_H3H4_raw_ab_msqrob_input.csv` | `H3_alignement.fas` | <code>sp&#124;P68431&#124;H31</code> | Human histone H3.1 |
| `liver_H4` | Liver | H4 | `260418_TroutLiver_H3H4_raw_ab_msqrob_input.csv` | `H4_alignement.fas` | <code>sp&#124;P62805&#124;H4</code> | Human histone H4 |

## Archive structure

```text
histone_alignement/
├── README.md
├── pymsaviz.ipynb
├── H3_alignement.fas
├── H4_alignement.fas
├── 20260420_TroutGut_peptides_ions_raw_ab_rank3_msqrob-input.csv
├── 260418_TroutLiver_H3H4_raw_ab_msqrob_input.csv
├── midgut_H3_MSA_outputs/
├── midgut_H4_MSA_outputs/
├── liver_H3_MSA_outputs/
└── liver_H4_MSA_outputs/
```

## Inputs and prior annotation

- `H3_alignement.fas` and `H4_alignement.fas` are aligned histone reference sequences exported in FASTA format.
- `*_raw_ab_msqrob_input.csv` files are the input tables used for the [midgut hPTM](../../05_Midgut_hPTMs/) and [liver hPTM](../../06_Liver_hPTMs/) msqrobPTM workflows. The same tables are used here to retrieve experimentally detected peptide sequences and their PTM annotations.

Before analysis, bare Mascot `Me` modification labels were renamed to `Me1`. This makes monomethylation explicit and distinguishes it unambiguously from `Me2` and `Me3`. This nomenclature correction did not alter peptide sequences, PTM coordinates, or quantitative values.

For this alignment visualisation, all detected and recognised PTM classes from the Mascot search are displayed before biological feature filtering used in downstream analyses.

### Prior peptide-conservation and PTM-site annotation

This notebook does not independently infer peptide conservation or assign canonical PTM positions. These steps were completed manually before visualisation.

All peptidoforms detected in this study were manually inspected against the corresponding rainbow trout and canonical human H3/H4 sequences. Peptide-level conservation and peptidoform assignments were cross-referenced with the comprehensive human histone peptide reference dataset reported by Fernandez-Rojas et al., including the standard and HeLa S3 peptide LC-TIMS-ToF MS/MS characteristics in Table 1.

After peptide conservation was verified, peptidoform-level observations were manually annotated to canonical histone PTM sites commonly reported in the literature using the histone PTM compendium of Huang et al. The input CSV files therefore already contain the validated human-reference position and residue associated with each PTM annotation.

The present notebook uses these pre-assigned coordinates to place PTM symbols on the alignment and verifies that each reported residue matches the selected human reference sequence. The cross-species alignment provides a visual confirmation of conservation and coverage; it is not used to infer or reassign PTM-site positions.

## Multiple-sequence alignment generation

Histone H3 and H4 protein sequences were aligned separately using the MUSCLE algorithm implemented in MEGA version 11.0.13.

Each multiple-sequence alignment was exported from MEGA in aligned FASTA format (`*.fas`). These aligned FASTA files were used directly as inputs to `pymsaviz.ipynb` for peptide-coverage and PTM-site visualisation.

## Workflow

The notebook performs the following steps:

1. Selects rows matching the requested histone (`Histone H3` or `Histone H4`) from the selected input CSV.
2. Retains valid unmodified peptide strings and optionally removes peptides listed in `excluded_peptides`.
3. Parses Mascot-style PTM annotations from `Variable modifications ([position] description)`, expected in the form `[position] (residue) modification`.
4. Normalises PTM labels to explicit plotting classes: `Ac`, `Me1`, `Me2`, `Me3`, `Bu`, `Prop`, `Cr`, `La`, `Ox`, `Deam`, and `Fo`.
5. Requires explicit methylation classes (`Me1`, `Me2`, or `Me3`) and stops if a bare `Me` annotation is detected.
6. Uses the pre-assigned canonical human H3/H4 PTM coordinates present in the input table to place PTM markers on the alignment, and verifies that the reported residue matches the selected human reference sequence.
7. Searches every retained peptide against every ungapped alignment sequence and colours the corresponding alignment residues.
8. Draws one PTM symbol for each distinct site × PTM-class combination. Multiple PTM classes at the same residue are vertically stacked.
9. Colours sequence labels by species while preserving the original protein identifier in the displayed label.

Only these input columns are used to construct the alignment figure:

- `Protein`
- `Sequence`
- `Variable modifications ([position] description)`

## Running the notebook

Open `pymsaviz.ipynb` in a Python/Jupyter environment and run all cells from top to bottom.

Required Python packages:

```text
pandas
matplotlib
biopython
pymsaviz
session-info
ipykernel
```

The notebook metadata records Python 3.12.6. Package versions are environment-specific and should be recorded when rerunning the analysis.

### Dataset-specific settings

Edit only the user-settings cell before each run:

```python
dataset_tag = "liver_H4"
csv_file = Path("260418_TroutLiver_H3H4_raw_ab_msqrob_input.csv")
msa_file = Path("H4_alignement.fas")
protein_name = "Histone H4"
site_prefix = "H4"
reference_id = "sp|P62805|H4"
```

Use `excluded_peptides = set()` when no additional peptide exclusion is required.

By default, `plot_ptm_types` contains all recognised PTM classes. It can be changed for a targeted display, while `*_PTM_annotation_audit.tsv` always records the raw-to-normalised PTM mapping.

## Generated outputs

Each run writes an output directory named `<dataset_tag>_MSA_outputs/`.

| File | Contents |
|---|---|
| `*_identified_peptides_and_PTM_sites.pdf` | Main alignment figure: experimentally detected peptides are colour-coded and PTM sites are marked by symbol. |
| `*_PTM_symbol_legend.pdf` | Standalone PTM symbol legend. |
| `*_identified_peptides.tsv` | Unique retained peptide sequences. |
| `*_all_PTM_annotations.tsv` | All parsed PTM annotations before PTM-class display filtering. |
| `*_PTM_annotation_audit.tsv` | Mapping from raw modification labels to normalised PTM classes. |
| `*_plotted_PTM_annotations.tsv` | PTM annotations retained for plotting, including reference residues and alignment columns. |
| `*_plotted_PTM_sites.tsv` | Distinct site × PTM combinations, marker order, and vertical stack offsets. |
| `*_alignment_sequence_info.tsv` | FASTA record IDs, complete headers, species assignments, display labels, and ungapped lengths. |
| `*_peptide_alignment_matches.tsv` | All exact peptide matches across alignment sequences and ungapped coordinates. |
| `*_peptide_colour_key.tsv` | Colour assigned to each peptide. |

The notebook also contains TIFF export commands at 1200 dpi. TIFF files are not included in this archive.

## Archived output summary

| Dataset | Unique peptides | Alignment sequences | Distinct PTM residues | Site × PTM combinations |
|---|---:|---:|---:|---:|
| `liver_H3` | 10 | 20 | 13 | 24 |
| `liver_H4` | 10 | 10 | 9 | 14 |
| `midgut_H3` | 10 | 20 | 15 | 29 |
| `midgut_H4` | 10 | 10 | 11 | 21 |

All retained peptides in each archived dataset matched at least one sequence in the corresponding alignment.

## Important checks and limitations

- PTM positions are pre-assigned relative to the selected canonical human H3 or H4 reference sequence during prior manual annotation; this notebook does not infer or reassign PTM coordinates.
- A PTM annotation whose reported residue differs from the selected human reference residue triggers an error rather than being plotted.
- Peptides are colour-mapped by exact ungapped amino-acid sequence match. This is an annotation and coverage visualisation; it does not quantify peptide abundance or PTM occupancy.
- At overlapping peptide positions, the displayed residue colour is assigned to the longest peptide. The complete set of matches remains available in `*_peptide_alignment_matches.tsv`.
- The sequence alignments support visual inspection of conservation. They do not establish histone-variant-specific abundance or protein inference.

## Reproducibility

The Python and package environment is printed at the end of the notebook using the `session-info` package:

```python
# Install once, if needed:
# pip install session-info

import session_info
session_info.show()
```

This records the Python version and versions of packages loaded in the current notebook session. Retain the resulting output with the corresponding archived analysis results.

## References

Fernandez-Rojas M, Fuller CN, Tose LV, et al. 2024. Histone Modification Screening Using Liquid Chromatography, Trapped Ion Mobility Spectrometry, and Time-of-Flight Mass Spectrometry. *Journal of Visualized Experiments*, 203: e65589. doi:10.3791/65589.

Huang H, Sabari BR, Garcia BA, Allis CD, Zhao Y. 2014. SnapShot: Histone Modifications. *Cell* 159(2): 458–458.e1. doi:10.1016/j.cell.2014.09.037.
