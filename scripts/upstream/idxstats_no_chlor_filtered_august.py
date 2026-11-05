# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Counts primary, properly paired reads with mapping quality >= 60 per giant sequoia
# chromosome (samtools view -F 256 -f 2 -q 60; samtools idxstats). These counts are n_ij in alpha_ij (Methods:
# Chromosome copy number; Figure 3; Supplementary Table S1). The input folder in the script was edited after the
# survey run.
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
import os
import subprocess
import pandas as pd

#Same script as for MQ30, but I changed the quality on August 15th 

# Directories
bam_dir = "/space/s1/sashanikolaeva/Results/Alignments/SESE/Amanda_to_SEGI"
output_dir = "/space/s1/sashanikolaeva/Results/Analysis/Ploidy_results/idxstats_results_AmandaMQ60"
output_file = os.path.join(output_dir, "read_counts_per_chromosome_AmandaMQ60.csv")

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# List BAM files
bam_files = [f for f in os.listdir(bam_dir) if f.endswith(".bam")]

# Initialize an empty DataFrame to store the results
results_df = pd.DataFrame()

# Iterate over each BAM file and run idxstats
for bam_file in bam_files:
    bam_path = os.path.join(bam_dir, bam_file)
    sample_name = bam_file.split(".bam")[0]
    
    # Filter reads by MQ and proper pairing, exclude secondary alignments
    cmd_filter = f"samtools view -b -F 256 -f 2 -q 60 {bam_path} | samtools idxstats -"
    result = subprocess.run(cmd_filter, shell=True, capture_output=True, text=True)
    
    # Parse the idxstats output
    sample_counts = {}
    for line in result.stdout.splitlines():
        fields = line.split("\t")
        scaffold = fields[0]
        mapped_reads = int(fields[2])
        sample_counts[scaffold] = mapped_reads
    
    # Convert the sample_counts dictionary to a DataFrame
    sample_df = pd.DataFrame.from_dict(sample_counts, orient='index', columns=[sample_name])
    
    # Transpose the DataFrame to have scaffolds as columns
    sample_df = sample_df.T
    
    # Append the sample_df to the results_df
    results_df = pd.concat([results_df, sample_df], axis=0)

# Fill any NaN values with 0
results_df = results_df.fillna(0)

# Save the final DataFrame to CSV
results_df.to_csv(output_file)

print(f"Read counts table saved to {output_file}")
