from pathlib import Path
from Bio import SeqIO
from collections import defaultdict

# Paths
ROOT_DIR = Path(__file__).resolve().parent
INPUT_DIR = ROOT_DIR / "02_input_files"
MS_HISTONEDB_DIR = ROOT_DIR / "01_FASTA_MS_HistoneDB"
OUTPUT_DIR = ROOT_DIR / "output_files"

# File paths to your input files
files = [
    INPUT_DIR / "24_10_04_histone_prot_seq_annotation_from_ENSEMBL.txt",
    INPUT_DIR / "24_10_02_histone_prot_seq_annotation_from_UP000694395.txt",
    INPUT_DIR / "24_10_02_histone_prot_seq_annotation_from_RefSEQ_NCBI_USDA_OmykA_1_1.txt",
    MS_HISTONEDB_DIR / "24_10_08_histone_seq_MS_HistoneDB.fasta",
    INPUT_DIR / "24_10_02_histone_prot_seq_annotation_from_UP000193380.txt",
]

# Dictionary to store sequences with headers (key=sequence, value=list of headers)
sequence_dict = defaultdict(list)

# Function to process files and store sequences in the dictionary
def process_file(filename, source):
    current_header = ""
    current_sequence = []

    with open(filename, "r") as file:
        for line in file:
            line = line.strip()

            if line.startswith("#"):  # Skip lines starting with #
                continue

            if line.startswith(">"):  # New header found
                if current_sequence:  # Save the previous sequence
                    seq_str = ''.join(current_sequence)
                    sequence_dict[seq_str].append(current_header)
                    current_sequence = []  # Reset sequence

                # Process the header based on the source
                if source == "uniprot_ref":  # For UP000193380
                    header_parts = line.split(" ", 1)
                    identifier = header_parts[0].split("|")[1]
                    description = (
                        header_parts[1].split(" OS=")[0]
                        if len(header_parts) > 1
                        else ""
                    )
                    current_header = (
                        f"{identifier} {description} [UniProt_ref]"
                    )

                elif source == "uniprot_ens":  # For UP000694395
                    header_parts = line.split(" ", 1)
                    identifier = header_parts[0].replace(">", "")
                    description = (
                        header_parts[1].split(" OS=")[0]
                        if len(header_parts) > 1
                        else ""
                    )
                    current_header = (
                        f"{identifier} {description} [UniProt_ens]"
                    )

                elif source == "MS_HistoneDB":  # For MS_HistoneDB
                    header_parts = line.split("|")
                    identifier = header_parts[0].replace(">", "")
                    description = (
                        header_parts[2]
                        if len(header_parts) > 2
                        else ""
                    )
                    current_header = (
                        f"{identifier} {description} [MS_HistoneDB]"
                    )

                elif source == "refseq":  # For RefSeq NCBI
                    header_parts = line.split(" ", 1)
                    identifier = header_parts[0].replace(">", "")
                    description = (
                        header_parts[1].split(" [")[0]
                        if len(header_parts) > 1
                        else ""
                    )
                    current_header = (
                        f"{identifier} {description} [RefSeq]"
                    )

                elif source == "ensembl":  # For Ensembl
                    header_parts = line.split(" ", 1)
                    identifier = (
                        header_parts[0]
                        .split(".")[0]
                        .replace(">", "")
                    )

                    description = ""

                    if "description:" in line:
                        desc_start = (
                            line.find("description:")
                            + len("description:")
                        )
                        description = (
                            line[desc_start:]
                            .split("[")[0]
                            .strip()
                        )

                    current_header = (
                        f"{identifier} {description} [Ensembl]"
                    )

            else:
                current_sequence.append(line)

        # Handle the last sequence in the file
        if current_sequence:
            seq_str = ''.join(current_sequence)
            sequence_dict[seq_str].append(current_header)


# Process all input files with their respective sources
process_file(
    INPUT_DIR / "24_10_02_histone_prot_seq_annotation_from_UP000193380.txt",
    "uniprot_ref"
)

process_file(
    INPUT_DIR / "24_10_02_histone_prot_seq_annotation_from_UP000694395.txt",
    "uniprot_ens"
)

process_file(
    INPUT_DIR / "24_10_02_histone_prot_seq_annotation_from_RefSEQ_NCBI_USDA_OmykA_1_1.txt",
    "refseq"
)

process_file(
    MS_HISTONEDB_DIR / "24_10_08_histone_seq_MS_HistoneDB.fasta",
    "MS_HistoneDB"
)

process_file(
    INPUT_DIR / "24_10_04_histone_prot_seq_annotation_from_ENSEMBL.txt",
    "ensembl"
)


# Write the non-redundant sequences to a new FASTA file
with open(
    OUTPUT_DIR / "24_10_09_new_non_redundant_histone_DB_sequences.fasta",
    "w"
) as output_file:

    for seq, headers in sequence_dict.items():
        merged_header = " | ".join(headers)
        output_file.write(f">{merged_header}\n{seq}\n\n")


# Organize sequences into classes (H1, H2A, H2B, H3, H4)
classes = defaultdict(list)

for seq, headers in sequence_dict.items():
    merged_header = " | ".join(headers)
    description = merged_header.lower()

    if "h1" in description:
        classes["H1"].append((merged_header, seq))

    elif "h2a" in description:
        classes["H2A"].append((merged_header, seq))

    elif "h2b" in description:
        classes["H2B"].append((merged_header, seq))

    elif "h3" in description:
        classes["H3"].append((merged_header, seq))

    elif "h4" in description:
        classes["H4"].append((merged_header, seq))


# Output the organized sequences by class
for class_name, sequences in classes.items():
    with open(
        OUTPUT_DIR / f"{class_name}_sequences.fasta",
        "w"
    ) as class_file:

        for header, seq in sequences:
            class_file.write(f">{header}\n{seq}\n\n")


# Generate statistics
total_sequences = len(sequence_dict)

unique_sequences = sum(
    1
    for headers in sequence_dict.values()
    if len(headers) == 1
)

duplicate_sequences = total_sequences - unique_sequences

print(f"Total sequences: {total_sequences}")
print(f"Unique sequences: {unique_sequences}")
print(
    "Sequences with 100% identity match in multiple files: "
    f"{duplicate_sequences}"
)


# Output summary of duplicates
with open(
    OUTPUT_DIR / "24_10_09_sequence_statistics.txt",
    "w"
) as stat_file:

    stat_file.write(f"Total sequences: {total_sequences}\n")

    stat_file.write(
        f"Unique sequences: {unique_sequences}\n"
    )

    stat_file.write(
        "Sequences with 100% identity match in multiple files: "
        f"{duplicate_sequences}\n"
    )

    for seq, headers in sequence_dict.items():
        if len(headers) > 1:
            stat_file.write(f"Sequence: {seq}\n")
            stat_file.write(
                f"Headers: {' | '.join(headers)}\n\n"
            )