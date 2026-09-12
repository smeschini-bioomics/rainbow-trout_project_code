import re
from Bio import SeqIO
from pathlib import Path

# Define pattern to identify histone variants
HISTONE_VARIANT_PATTERN = re.compile(r"H\d[A-Z]\.\d|H\d[A-Z]|H\d\.\d|H\d")

# Define keywords that indicate specificity
SPECIFICITY_KEYWORDS = {"variant": 10, "isoform": 10, "centromeric": 5, "family": 3, "like": 2}

# Scoring function for specificity
def calculate_score(description):
    score = 0
    for keyword, points in SPECIFICITY_KEYWORDS.items():
        if keyword.lower() in description.lower():
            score += points
    if HISTONE_VARIANT_PATTERN.search(description):
        score += 20
    score += len(description) / 10
    return score

# Function to parse descriptions and extract ID, description, and database
def parse_description(description):
    match = re.match(r"(\S+)\s(.+?)\s\[(\w+)\]", description)
    if match:
        id_part = match.group(1)  # Keep the original ID exactly as it is
        desc_part = match.group(2)
        db_part = match.group(3)
    else:
        # When no match is found, try to separate by space, taking the first as the ID
        id_part = description.split()[0]
        desc_part = " ".join(description.split()[1:])
        db_part = ""
    return id_part, desc_part, db_part

# Function to choose the best description
def choose_best_description(descriptions):
    scored_descriptions = []
    for desc in descriptions:
        id_part, desc_part, db_part = parse_description(desc)
        score = calculate_score(desc_part)
        scored_descriptions.append((id_part, desc_part, db_part, score))
    scored_descriptions.sort(key=lambda x: (x[3], len(x[1])), reverse=True)
    best_id, best_desc, best_db, _ = scored_descriptions[0]
    return best_id, best_desc, best_db

# Function to simplify FASTA headers
def simplify_fasta_headers(input_path, output_path):
    with open(input_path, "r") as input_file, open(output_path, "w") as output_file:
        for record in SeqIO.parse(input_file, "fasta"):
            # Keep original description parts separated by " | " and find the best
            original_description = record.description.split(" | ")
            best_id, best_description, best_db = choose_best_description(original_description)
            # Only include the database part if it exists
            db_info = f" [{best_db}]" if best_db else ""
            # Construct the simplified header, ensuring it's properly formatted
            simplified_header = f">{best_id} {best_description}{db_info}"
            output_file.write(simplified_header + "\n")
            output_file.write(str(record.seq) + "\n\n")

# Input and output paths
# Input and output paths
ROOT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = ROOT_DIR / "output_files"

input_file_path = (
    OUTPUT_DIR /
    "24_10_09_new_non_redundant_histone_DB_sequences.fasta"
)

output_file_path = (
    OUTPUT_DIR /
    "24_10_09_histones_sequences_header_simplified_mascot_input.fasta"
)

# Run the simplification
simplify_fasta_headers(input_file_path, output_file_path)
