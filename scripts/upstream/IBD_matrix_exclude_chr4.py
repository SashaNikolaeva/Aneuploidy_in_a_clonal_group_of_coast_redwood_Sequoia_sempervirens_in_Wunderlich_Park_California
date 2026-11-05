#!/usr/bin/env python3
# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: The same distance matrix without SNPs placed on chromosome 4 (Methods: check that the
# aneuploid chromosome does not drive the distances).
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
"""
Calculate IBD distance matrix excluding SNPs that map to SEGI chr4.
Based on IBD_matrix.py but filters out chr4 SNPs before calculation.
"""
import numpy as np
import pandas as pd

# --- Paths ---
dosage_path = "/space/s1/sashanikolaeva/Results/Population_structure/April2025files/dosage_matrix.tsv"
chr4_snps_path = "/space/s1/sashanikolaeva/Results/clone/snp_segi_alignment/chr4_snps.txt"
out_path = "/space/s1/sashanikolaeva/Results/clone/genetic_distance_matrix_IBD_no_chr4.csv"

print("Loading chr4 SNPs to exclude...")
chr4_df = pd.read_csv(chr4_snps_path, sep='\t')
# Create SNP IDs in the format used by dosage matrix: scaffold_position
chr4_snps = set(chr4_df['scaffold'] + '_' + chr4_df['position'].astype(str))
print(f"Loaded {len(chr4_snps)} chr4 SNPs to exclude")

print("Loading dosage matrix...")
dosage_df = pd.read_csv(dosage_path, sep="\t", index_col=0)
print(f"Original dosage matrix: {dosage_df.shape[0]} samples x {dosage_df.shape[1]} SNPs")

# Filter out chr4 SNPs
snps_to_keep = [col for col in dosage_df.columns if col not in chr4_snps]
dosage_df = dosage_df[snps_to_keep]
print(f"After excluding chr4: {dosage_df.shape[0]} samples x {dosage_df.shape[1]} SNPs")
print(f"Excluded {len(chr4_snps) - (len(chr4_snps) - (334029 - len(snps_to_keep)))} SNPs")

sample_names = dosage_df.index.to_numpy()
X = dosage_df.to_numpy(dtype=np.float64)        # (n_samples, n_snps)

print("Calculating IBD distance matrix...")

# --- Masks and NaN handling ---
V = (~np.isnan(X)).astype(np.float64)
X0 = np.nan_to_num(X, nan=0.0)

# --- Pairwise counts: n_ij (# loci where both valid) ---
n_ij = V @ V.T                                  # (n_samples, n_samples)

# --- Sums over intersection for each pair ---
A = X0 @ V.T                                    # (n_samples, n_samples)
P = X0 @ X0.T                                   # (n_samples, n_samples)
numerator = A + A.T - 2.0 * P                   # (n_samples, n_samples)

# --- Distance matrix: D_ij = numerator / n_ij ---
with np.errstate(invalid="ignore", divide="ignore"):
    D = numerator / n_ij

# Set pairs with no shared valid sites to NaN
D[~np.isfinite(D)] = np.nan

# Enforce symmetry and zero diagonal
np.fill_diagonal(D, 0.0)

# --- Save as CSV ---
dist_df = pd.DataFrame(D, index=sample_names, columns=sample_names)
dist_df.to_csv(out_path)

print(f"\nDistance matrix written to: {out_path}")
print(f"Matrix size: {dist_df.shape[0]} x {dist_df.shape[1]}")
