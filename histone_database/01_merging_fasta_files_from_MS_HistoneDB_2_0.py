from pathlib import Path
from Bio import SeqIO


# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------

ROOT_DIR = Path(__file__).resolve().parent
FASTA_DIR = ROOT_DIR / "01_FASTA_MS_HistoneDB"

# Keep the original historical order: H1, H2B, H2A, H3, H4
fasta_files = [
    FASTA_DIR / "histone_variants_H1_MS_HistoneDB.fasta",
    FASTA_DIR / "histone_variants_H2B_MS_HistoneDB.fasta",
    FASTA_DIR / "histone_variants_H2A_MS_HistoneDB.fasta",
    FASTA_DIR / "histone_variants_H3_MS_HistoneDB.fasta",
    FASTA_DIR / "histone_variants_H4_MS_HistoneDB.fasta",
]

output_file = (
    FASTA_DIR /
    "24_10_08_histone_seq_MS_HistoneDB.fasta"
)


# -----------------------------------------------------------------------------
# Check inputs
# -----------------------------------------------------------------------------

for fasta_file in fasta_files:
    if not fasta_file.exists():
        raise FileNotFoundError(f"Input FASTA not found:\n{fasta_file}")


# -----------------------------------------------------------------------------
# Merge FASTA files
# -----------------------------------------------------------------------------

n_records = 0

with output_file.open("w", encoding="utf-8") as outfile:

    for fasta_file in fasta_files:

        for record in SeqIO.parse(fasta_file, "fasta"):
            SeqIO.write(record, outfile, "fasta")
            n_records += 1


print(f"FASTA files merged into: {output_file}")
print(f"Sequences written: {n_records}")