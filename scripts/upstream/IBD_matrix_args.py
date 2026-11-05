#!/usr/bin/env python3
# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Identity-by-descent distance matrix from the allele-dosage matrix (PCR-corrected and
# uncorrected). Distances of Figure 2 and of the clonality tests (Methods: Genetic distance matrix).
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
# IBD distance matrix from an allele-dosage matrix. Same computation as
# Coast_Redwood_Population_Structure_dev/scripts/publication_224samples/IBD_matrix.py; the only change is that the
# input and output paths are arguments (the original has them hard-coded, and the PCR-corrected matrix was made
# by editing them by hand without saving the edit).
#
#   uncorrected matrix (genetic_distance_matrix_IBD.csv, 2025-08-19):
#     IBD_matrix_args.py <April2025files>/dosage_matrix.tsv <April2025files>/genetic_distance_matrix_IBD.csv
#   PCR-corrected matrix used for the clone paper (genetic_distance_matrix_IBD_corrected.csv, 2026-05-14):
#     IBD_matrix_args.py <April2025files>/dosage_matrix_pool_corrected.tsv <April2025files>/genetic_distance_matrix_IBD_corrected.csv
#   where <April2025files> = /space/s1/sashanikolaeva/Results/Population_structure/April2025files
# The dosage_matrix_pool_corrected.tsv (2026-05-09) is written by pool_correction.py in
# Coast_Redwood_Population_Structure_dev/scripts/02_quality_control/.
import sys
import numpy as np
import pandas as pd

if len(sys.argv) != 3:
    sys.exit("usage: IBD_matrix_args.py DOSAGE_MATRIX.tsv OUT_DISTANCE_MATRIX.csv")
dosage_path, out_path = sys.argv[1], sys.argv[2]

# --- Load dosage matrix ---
# rows: samples, columns: loci/SNPs; entries are dosage (can contain NaN)
dosage_df = pd.read_csv(dosage_path, sep="\t", index_col=0)

sample_names = dosage_df.index.to_numpy()
X = dosage_df.to_numpy(dtype=np.float64)        # (n_samples, n_snps)

# --- Masks and NaN handling ---
V = (~np.isnan(X)).astype(np.float64)           
X0 = np.nan_to_num(X, nan=0.0)                  

# --- Pairwise counts: n_ij (# loci where both valid) ---
# n = V @ V^T
n_ij = V @ V.T                                  # (n_samples, n_samples)

# --- Sums over intersection for each pair ---
# sum_i_over_intersection = X0[i] dot V[j]
# A[i,j] = sum_s X0[i,s] * V[j,s]
A = X0 @ V.T                                    # (n_samples, n_samples)
# A.T[j,i] = sum_j_over_intersection
P = X0 @ X0.T                                   # (n_samples, n_samples)
numerator = A + A.T - 2.0 * P                   # (n_samples, n_samples)

# --- Distance matrix: D_ij = numerator / n_ij ---
with np.errstate(invalid="ignore", divide="ignore"):
    D = numerator / n_ij

# Set pairs with no shared valid sites to NaN 
D[~np.isfinite(D)] = np.nan

# Enforce symmetry and zero diagonal (distance of a sample to itself)
# Note: the raw formula yields 2*f*(1-f) for i==j (not zero), so we zero it explicitly.
np.fill_diagonal(D, 0.0)

# --- Save as CSV ---
dist_df = pd.DataFrame(D, index=sample_names, columns=sample_names)
dist_df.to_csv(out_path)

print(f"Distance matrix written to: {out_path}")
