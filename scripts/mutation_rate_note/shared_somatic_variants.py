#!/usr/bin/env python3
# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Alleles carried by all nine stems at low fraction and rare in other trees, candidate
# mutations of the founding tree (Supplementary Note S1: 14 candidates, 13 inherited).
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
"""
Step 2: Identify shared clone-specific somatic variants.

Scans the population VCF for sites where all 9 clone ramets carry an ALT allele
at ~1/6 dosage (consistent with one copy out of six in a hexaploid), while the
allele is absent or very rare in non-clone samples.

These candidates represent somatic mutations that arose in the founding tree before
clonal reproduction began. Their count (after FP correction) divided by the
self-calibrated mutation rate gives an estimate of founding tree age.

Usage:
    conda run -n ntsynt-310 python3 shared_somatic_variants.py

Runtime: ~2–4 hours (VCF traversal, 284 samples, 334,029 SNPs)
         Step 3 permutation null: ~1–2 days (parallelized)

Outputs in /space/s1/sashanikolaeva/Results/clone/somatic_mutations_neale/:
    shared_clone_variants.tsv       - candidate shared somatic variants
    shared_clone_summary.txt        - N_shared + FP correction + age estimate
    shared_clone_af_dist.png        - ALT freq distribution in non-clone samples
"""

import os
import subprocess
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ─── Paths ────────────────────────────────────────────────────────────────────
VCF     = "/space/s1/santang3/prg/results/variants/vcf/allScaffs_allRegions_filtered.vcf.gz"
OUT_DIR = "/space/s1/sashanikolaeva/Results/clone/somatic_mutations_neale"

# ─── Sample lists ─────────────────────────────────────────────────────────────
CLONE_SAMPLES = [
    "RedwoodCap_13_WP1",
    "RedwoodCap_6_WP5",
    "RedwoodCap_31_c3",
    "RedwoodCap_31_WP2_1",
    "RedwoodCap_28_WP7",
    "RedwoodCap_28_WP6",
    "RedwoodCap_14_WP3",
    "RedwoodCap_22_WPd2",
    "RedwoodCap_33_c6",
]

# Exclude clone ramets + same-site samples from non-clone set
EXCLUDE_FROM_NONCLONE = set(CLONE_SAMPLES) | {
    "RedwoodCap_22_WPR",  # same site, not a clone ramet
    "RedwoodCap_16_C1",   # control from Wunderlich area
}

# ─── Filters ──────────────────────────────────────────────────────────────────
MIN_DP_CLONE    = 10    # minimum depth in each clone ramet
MIN_AO_CLONE    = 2     # minimum ALT reads in each clone ramet
MIN_AF_CLONE    = 0.08  # minimum ALT freq in clone ramet (1/6 = 0.167; allow 0.08 floor)
MAX_AF_CLONE    = 0.30  # maximum ALT freq in clone ramet (>2 copies → population variant)
MAX_AO_NONCLONE = 1     # ALT reads allowed in non-clone samples (0 = strict, 1 = lenient)
MAX_AC_NONCLONE = 5     # max number of non-clone samples that can carry the ALT


def get_sample_indices(header_line):
    """Return dicts: sample_name -> column index (0-based from sample columns)."""
    cols = header_line.rstrip("\n").split("\t")
    # cols[0:9] are CHROM POS ID REF ALT QUAL FILTER INFO FORMAT
    sample_cols = cols[9:]
    return {s: i for i, s in enumerate(sample_cols)}, sample_cols


def parse_ad(fmt_str, gt_str):
    """
    Extract (DP, AO) from a FreeBayes GT string.
    FORMAT: GT:DP:AD:RO:QR:AO:QA
    Returns (dp, ao) as ints. Returns (0, 0) on parse failure.
    """
    fields = gt_str.split(":")
    if len(fields) < 7:
        return 0, 0
    try:
        dp = int(fields[1])
        ao = int(fields[5])   # AO = alt observations
        return dp, ao
    except ValueError:
        return 0, 0


def run_scan():
    print("Opening VCF with bcftools view...")
    cmd = ["bcftools", "view", VCF]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)

    sample_idx = None
    nonclone_idx = None
    clone_idx = None
    all_samples = None

    candidates = []
    n_sites = 0
    n_pass_depth = 0

    for line in proc.stdout:
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            sample_idx, all_samples = get_sample_indices(line)
            clone_idx = [sample_idx[s] for s in CLONE_SAMPLES if s in sample_idx]
            nonclone_idx = [
                i for s, i in sample_idx.items()
                if s not in EXCLUDE_FROM_NONCLONE
            ]
            n_nonclone = len(nonclone_idx)
            print(f"  {len(clone_idx)} clone samples | {n_nonclone} non-clone samples")
            assert len(clone_idx) == 9, f"Expected 9 clone samples, got {len(clone_idx)}"
            continue

        parts = line.rstrip("\n").split("\t")
        chrom, pos, ref, alt = parts[0], int(parts[1]), parts[3], parts[4]
        fmt = parts[8]
        gt_cols = parts[9:]

        n_sites += 1
        if n_sites % 50_000 == 0:
            print(f"  {n_sites:,} sites | {len(candidates)} candidates so far")

        # ── Clone ramet filters ──────────────────────────────────────────────
        clone_ok = True
        clone_ao_sum = 0
        clone_dp_sum = 0
        for ci in clone_idx:
            dp, ao = parse_ad(fmt, gt_cols[ci])
            if dp < MIN_DP_CLONE:
                clone_ok = False
                break
            af = ao / dp if dp > 0 else 0
            if ao < MIN_AO_CLONE or af < MIN_AF_CLONE or af > MAX_AF_CLONE:
                clone_ok = False
                break
            clone_ao_sum += ao
            clone_dp_sum += dp

        if not clone_ok:
            continue
        n_pass_depth += 1

        # ── Non-clone filter ─────────────────────────────────────────────────
        n_nonclone_carry = 0
        max_ao_nonclone = 0
        for ni in nonclone_idx:
            dp, ao = parse_ad(fmt, gt_cols[ni])
            if dp < MIN_DP_CLONE:
                continue  # not callable — don't penalize
            if ao > MAX_AO_NONCLONE:
                n_nonclone_carry += 1
            if ao > max_ao_nonclone:
                max_ao_nonclone = ao

        if n_nonclone_carry > MAX_AC_NONCLONE:
            continue

        mean_clone_af = (clone_ao_sum / clone_dp_sum) if clone_dp_sum > 0 else 0
        candidates.append({
            "chrom":             chrom,
            "pos":               pos,
            "ref":               ref,
            "alt":               alt,
            "mean_clone_af":     round(mean_clone_af, 3),
            "n_nonclone_carry":  n_nonclone_carry,
            "max_ao_nonclone":   max_ao_nonclone,
        })

    proc.wait()
    print(f"\n  Total sites: {n_sites:,} | "
          f"Pass clone depth+AF: {n_pass_depth:,} | "
          f"Final candidates: {len(candidates)}")
    return candidates, len(nonclone_idx)


def permutation_null(n_perm=10_000):
    """
    Estimate false positives from rare inherited alleles using a permutation null.
    Draw 9 random non-clone samples and apply the same filters.

    Vectorized: all sites processed as matrices per permutation.
    Memory: ~350 MB for (n_sites × n_nonclone) int16 matrices.
    """
    print(f"\n=== Permutation null ({n_perm:,} iterations) ===")
    print("  Loading VCF into memory (334K SNPs × non-clone samples)...")

    cmd = ["bcftools", "view", VCF]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)

    sample_idx = None
    nonclone_idx_list = None
    dp_rows = []
    ao_rows = []

    for line in proc.stdout:
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            sample_idx, all_samples = get_sample_indices(line)
            nonclone_idx_list = [
                i for s, i in sample_idx.items()
                if s not in EXCLUDE_FROM_NONCLONE
            ]
            continue

        parts = line.rstrip("\n").split("\t")
        fmt = parts[8]
        gt_cols = parts[9:]

        dp_row = []
        ao_row = []
        for ni in nonclone_idx_list:
            dp, ao = parse_ad(fmt, gt_cols[ni])
            dp_row.append(dp)
            ao_row.append(ao)

        dp_rows.append(dp_row)
        ao_rows.append(ao_row)

    proc.wait()

    # Shape: (n_sites, n_nonclone) — use int16 to cap memory
    dp_mat = np.array(dp_rows, dtype=np.int16)
    ao_mat = np.array(ao_rows, dtype=np.int16)
    n_sites, n_nonclone = dp_mat.shape
    print(f"  Matrix: {n_sites:,} sites × {n_nonclone} non-clone samples "
          f"({dp_mat.nbytes + ao_mat.nbytes} bytes)")

    rng = np.random.default_rng(42)
    perm_counts = []

    for perm in range(n_perm):
        if perm % 1000 == 0:
            print(f"  Permutation {perm:,}/{n_perm:,}")

        chosen = rng.choice(n_nonclone, 9, replace=False)

        # (n_sites, 9) for the 9 "pseudo-clone" samples
        dp9 = dp_mat[:, chosen].astype(np.float32)
        ao9 = ao_mat[:, chosen].astype(np.float32)

        callable_9 = np.all(dp9 >= MIN_DP_CLONE, axis=1)        # (n_sites,)
        ao_ok      = np.all(ao9 >= MIN_AO_CLONE, axis=1)
        af9 = ao9 / np.maximum(dp9, 1)
        af_ok = np.all((af9 >= MIN_AF_CLONE) & (af9 <= MAX_AF_CLONE), axis=1)
        clone_pass = callable_9 & ao_ok & af_ok                  # (n_sites,)

        # Non-clone: everyone except chosen
        mask = np.ones(n_nonclone, dtype=bool)
        mask[chosen] = False
        dp_rest = dp_mat[:, mask].astype(np.float32)             # (n_sites, n_nonclone-9)
        ao_rest = ao_mat[:, mask].astype(np.float32)
        callable_rest = dp_rest >= MIN_DP_CLONE
        n_carry = np.sum((ao_rest > MAX_AO_NONCLONE) & callable_rest, axis=1)  # (n_sites,)

        n_pass = int(np.sum(clone_pass & (n_carry <= MAX_AC_NONCLONE)))
        perm_counts.append(n_pass)

    mean_fp = float(np.mean(perm_counts))
    ci_lo   = float(np.percentile(perm_counts, 2.5))
    ci_hi   = float(np.percentile(perm_counts, 97.5))
    print(f"  Mean FP (inherited alleles): {mean_fp:.1f} (95% CI: {ci_lo:.1f}–{ci_hi:.1f})")
    return mean_fp, ci_lo, ci_hi


def estimate_founding_age(n_shared, n_fp_perm, rate_summary_tsv):
    """
    Compute founding tree age from corrected shared variant count.
    Uses mutation rates from somatic_mutations_neale.py output.
    """
    if not os.path.exists(rate_summary_tsv):
        print(f"\n  WARNING: {rate_summary_tsv} not found — run somatic_mutations_neale.py first")
        return None

    df_rate = pd.read_csv(rate_summary_tsv, sep="\t")
    results = []
    n_true_shared = max(0, n_shared - n_fp_perm)

    EXOME_BP = 152_537_171  # callable denominator: full exome target (Step 1 will refine)
    PLOIDY = 6

    print(f"\n=== Age Estimation ===")
    print(f"  N_shared candidates: {n_shared}")
    print(f"  FP (permutation):    {n_fp_perm:.1f}")
    print(f"  N_true_shared:       {n_true_shared:.1f}")

    for _, row in df_rate.iterrows():
        try:
            mu = float(row["mu_per_site_per_yr"])
        except (ValueError, KeyError):
            continue
        # age = n_true / (mu * callable_bases * ploidy)
        # But mu here already accounts for ploidy (see somatic_mutations_neale.py)
        # age = n_true / (mu * ploidy * callable_bases)
        callable_bases = float(row.get("mean_callable", EXOME_BP))
        age = n_true_shared / (mu * PLOIDY * callable_bases) if mu > 0 and callable_bases > 0 else None
        results.append({
            "ramet_age_assumption": row["ramet_age_yr"],
            "exclude_pool31": row["exclude_pool31"],
            "mu_per_site_per_yr": mu,
            "n_true_shared": n_true_shared,
            "founding_tree_age_yr": round(age, 0) if age else "N/A",
        })
        if age:
            print(f"  Age (ramet={row['ramet_age_yr']}yr, excl_pool31={row['exclude_pool31']}): "
                  f"{age:.0f} years")

    return pd.DataFrame(results)


def plot_af_dist(candidates, out_path):
    if not candidates:
        return
    df = pd.DataFrame(candidates)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(df["mean_clone_af"], bins=30, color="#2980b9", edgecolor="white")
    ax.axvline(1/6, color="red", linestyle="--", label="Expected 1/6 dosage")
    ax.set_xlabel("Mean ALT frequency in clone ramets")
    ax.set_ylabel("Sites")
    ax.set_title(f"Shared clone-specific variant candidates (N={len(df)})\n"
                 f"mean clone AF distribution")
    ax.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out_path}")


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-perm", action="store_true",
                        help="Skip permutation null (faster, for testing)")
    parser.add_argument("--n-perm", type=int, default=10_000,
                        help="Number of permutations (default: 10000)")
    args = parser.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)

    # Step 1: Scan for candidates
    print("=== Step 1: Shared clone-specific variant scan ===")
    candidates, n_nonclone = run_scan()

    df_cand = pd.DataFrame(candidates)
    df_cand.to_csv(f"{OUT_DIR}/shared_clone_variants.tsv", sep="\t", index=False)
    print(f"\n  Saved {len(df_cand)} candidates to shared_clone_variants.tsv")

    # Step 2: Permutation null (optional)
    mean_fp_perm = 0.0
    if not args.skip_perm and len(df_cand) > 0:
        mean_fp_perm, fp_lo, fp_hi = permutation_null(n_perm=args.n_perm)

    # Step 3: Age estimation (requires somatic_mutations_neale.py output)
    rate_tsv = f"{OUT_DIR}/mutation_rate_summary.tsv"
    df_age = estimate_founding_age(len(df_cand), mean_fp_perm, rate_tsv)

    # Step 4: Write summary
    summary_path = f"{OUT_DIR}/shared_clone_summary.txt"
    with open(summary_path, "w") as f:
        f.write("=== Shared clone-specific somatic variant summary ===\n\n")
        f.write(f"N_shared candidates:   {len(df_cand)}\n")
        if not args.skip_perm:
            f.write(f"FP (permutation null): {mean_fp_perm:.1f} "
                    f"(95% CI: {fp_lo:.1f}–{fp_hi:.1f})\n")
            f.write(f"N_true_shared:         {max(0, len(df_cand) - mean_fp_perm):.1f}\n")
        f.write(f"\nFilters used:\n")
        f.write(f"  Clone: DP >= {MIN_DP_CLONE}, AO >= {MIN_AO_CLONE}, "
                f"AF {MIN_AF_CLONE}–{MAX_AF_CLONE}\n")
        f.write(f"  Non-clone: AO <= {MAX_AO_NONCLONE} in <= {MAX_AC_NONCLONE} samples\n")
        if df_age is not None:
            f.write("\n=== Founding tree age estimates ===\n")
            f.write(df_age.to_string(index=False))
            f.write("\n")

    print(f"\n  Summary written to {summary_path}")

    # Step 5: Plot
    plot_af_dist(candidates, f"{OUT_DIR}/shared_clone_af_dist.png")

    print(f"\nDone. Outputs in {OUT_DIR}/")


if __name__ == "__main__":
    main()
