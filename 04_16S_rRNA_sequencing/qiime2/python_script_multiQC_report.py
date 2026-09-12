import os
import subprocess
import shutil

# Original FastQC reports directory
fastqc_output_dir = "/mnt/c/Users/smeschini/Documents/DATA_INRAE/THESE/RESULTS/STPN2309_16S/STPN2309_16S_QIIME2/fastqc_reports_trimmed_seqs_2"
# Output directories for MultiQC reports
multiqc_output_dir_fwd = "/mnt/c/Users/smeschini/Documents/DATA_INRAE/THESE/RESULTS/STPN2309_16S/STPN2309_16S_QIIME2/multiqc_report_trimmed-seqs_2_forward"
multiqc_output_dir_rev = "/mnt/c/Users/smeschini/Documents/DATA_INRAE/THESE/RESULTS/STPN2309_16S/STPN2309_16S_QIIME2/multiqc_report_trimmed-seqs_2_reverse"

# Temporary directories to store separated FastQC reports
temp_fwd_dir = "/mnt/c/Users/smeschini/Documents/DATA_INRAE/THESE/RESULTS/STPN2309_16S/STPN2309_16S_QIIME2/temp_fastqc_fwd"
temp_rev_dir = "/mnt/c/Users/smeschini/Documents/DATA_INRAE/THESE/RESULTS/STPN2309_16S/STPN2309_16S_QIIME2/temp_fastqc_rev"

# Create directories if they don’t exist, cleaning up any old content
for dir_path in [temp_fwd_dir, temp_rev_dir, multiqc_output_dir_fwd, multiqc_output_dir_rev]:
    if os.path.exists(dir_path):
        shutil.rmtree(dir_path)  # Clear out the directory if it exists
    os.makedirs(dir_path, exist_ok=True)

# Separate FastQC reports into forward and reverse directories
print("Organizing files into forward and reverse directories...")
for filename in os.listdir(fastqc_output_dir):
    file_path = os.path.join(fastqc_output_dir, filename)
    if filename.endswith("_R1_001_fastqc.html") or filename.endswith("_R1_001_fastqc.zip"):
        print(f"Copying {filename} to forward directory")
        shutil.copy(file_path, temp_fwd_dir)
    elif filename.endswith("_R2_001_fastqc.html") or filename.endswith("_R2_001_fastqc.zip"):
        print(f"Copying {filename} to reverse directory")
        shutil.copy(file_path, temp_rev_dir)

# Run MultiQC for forward reads
print("Running MultiQC for forward reads...")
multiqc_cmd_fwd = ["multiqc", temp_fwd_dir, "-o", multiqc_output_dir_fwd]
subprocess.run(multiqc_cmd_fwd)

# Run MultiQC for reverse reads
print("Running MultiQC for reverse reads...")
multiqc_cmd_rev = ["multiqc", temp_rev_dir, "-o", multiqc_output_dir_rev]
subprocess.run(multiqc_cmd_rev)

# Cleanup temporary directories
shutil.rmtree(temp_fwd_dir)
shutil.rmtree(temp_rev_dir)

print(f"MultiQC report for forward reads generated at: {multiqc_output_dir_fwd}")
print(f"MultiQC report for reverse reads generated at: {multiqc_output_dir_rev}")
