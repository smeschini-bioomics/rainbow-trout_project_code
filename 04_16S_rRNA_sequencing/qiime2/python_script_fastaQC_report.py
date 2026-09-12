import os
import subprocess

# Define directories!
trimmed_dir = "/mnt/c/Users/smeschini/Documents/DATA_INRAE/THESE/RESULTS/STPN2309_16S/STPN2309_16S_QIIME2/exported_fastqs_primers_adapters_poly_tail_180_ml_multiqc-trimmed"  # Directory with trimmed FASTQ files
fastqc_output_dir = "/mnt/c/Users/smeschini/Documents/DATA_INRAE/THESE/RESULTS/STPN2309_16S/STPN2309_16S_QIIME2/fastqc_reports_trimmed_seqs_2"  # Directory to save FastQC reports

# Create FastQC output directory if it doesn't exist
os.makedirs(fastqc_output_dir, exist_ok=True)

# Find and separate R1 and R2 files
r1_files = sorted([os.path.join(trimmed_dir, f) for f in os.listdir(trimmed_dir) if "_R1_" in f])
r2_files = sorted([os.path.join(trimmed_dir, f) for f in os.listdir(trimmed_dir) if "_R2_" in f])

# List of all files to process with FastQC (R1 and R2 files only)
all_files = r1_files + r2_files

# Run FastQC on each individual file with specified threads
for fastq_file in all_files:
    print(f"Running FastQC on {fastq_file} with 12 threads...")
    fastqc_cmd = [
        "fastqc",
        fastq_file,
        "--outdir", fastqc_output_dir,
        "--threads", "12"  # Specify 12 threads
    ]
    subprocess.run(fastqc_cmd)

print("FastQC analysis complete. Reports are saved in:", fastqc_output_dir)