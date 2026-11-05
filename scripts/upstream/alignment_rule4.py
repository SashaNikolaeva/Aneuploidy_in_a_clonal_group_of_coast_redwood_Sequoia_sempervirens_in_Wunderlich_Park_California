# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Aligns the survey reads to the giant sequoia genome v2.0 (BWA-MEM; samtools markdup -r
# removes duplicates). Its BAMs are the input of the per-chromosome read counts (Methods: Chromosome copy
# number). Identified from the output file names; the BAMs have since been deleted.
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
import os
import subprocess
import re
from concurrent.futures import ThreadPoolExecutor

# Define paths
raw_reads_dir = "/space/s1/sashanikolaeva/Data/Full_sequencing_data"
reference_file = "/space/s1/sashanikolaeva/Data/SEGI_genome/SEGI.2.0.fa"
output_dir = "/space/s1/sashanikolaeva/Results/Alignments/SESE/Full_data_to_SEGI"

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Clear old output files
def clear_output_dir(directory):
    for f in os.listdir(directory):
        try:
            os.remove(os.path.join(directory, f))
        except Exception as e:
            print(f"Error removing file {f}: {e}")

clear_output_dir(output_dir)

# Index the reference genome (do this once before processing reads)

#subprocess.run(["bwa", "index", reference_file])

# Function to extract the common prefix from the file names
def extract_prefix(filename):
    match = re.match(r"(.*)(_R[12]_concatenated\.fastq\.gz)", filename)
    return match.group(1) if match else None

# Find all read files and group them into pairs
read_files = [f for f in os.listdir(raw_reads_dir) if f.endswith('.fastq.gz')]
paired_end_reads = {}
for f in read_files:
    prefix = extract_prefix(f)
    if prefix:
        if prefix not in paired_end_reads:
            paired_end_reads[prefix] = []
        paired_end_reads[prefix].append(f)

# Function to generate the command for processing a pair of reads
def generate_command(read_prefix, read1, read2):
    final_bam = os.path.join(output_dir, f"{read_prefix}_SEGI_full.bam")
    command = (
        f"bwa mem -t 30 {reference_file} {read1} {read2} | "
        f"samtools view -@ 30 -b - | "
        f"samtools collate -@ 30 -O - | "
        f"samtools fixmate -@ 30 -m - - | "
        f"samtools sort -@ 30 -o - | "
        f"samtools markdup -@ 30 -r - {final_bam} && "
        f"samtools index {final_bam}"
    )
    return command

# Function to process a pair of reads
def process_reads(read_prefix, reads):
    if len(reads) != 2:
        print(f"Warning: Expected 2 files for {read_prefix}, found {len(reads)}")
        return

    read1, read2 = sorted(reads)
    read1 = os.path.join(raw_reads_dir, read1)
    read2 = os.path.join(raw_reads_dir, read2)
    command = generate_command(read_prefix, read1, read2)
    
    print(f"Processing {read_prefix}...")
    subprocess.run(command, shell=True)
    print(f"Completed {read_prefix}.")

# Process reads in parallel
with ThreadPoolExecutor(max_workers=2) as executor:
    futures = [executor.submit(process_reads, read_prefix, reads) for read_prefix, reads in paired_end_reads.items()]

# Ensure all threads are completed
for future in futures:
    future.result()

print("Alignment and processing complete.")
