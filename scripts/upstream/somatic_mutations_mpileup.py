#!/usr/bin/env python3
# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Search for alleles carried by one stem only (samtools mpileup on the 13 libraries above;
# Supplementary Note S1). Re-run on the current BAMs as scripts/mutation_rate_note/mr02_scan_original.py.
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
"""
BAM-based somatic mutation calling for coast redwood clone ramets.

Uses samtools mpileup directly on the 9 clone BAM files (aligned to SEGI v2.0).
Identifies positions where exactly 1 ramet carries an ALT allele and all others
carry only REF. Avoids multi-sample VCF batch effects (PCR pool artifacts).

False positive calibration uses WITHIN-POOL pairs:
  - Pool 28: 28_WP6 vs 28_WP7 (same pool, same genotype → differences = technical error)
  - Pool 31: 31_c3 vs 31_WP2_1 (same pool, same genotype)

Usage:
    conda run -n ntsynt-310 python3 somatic_mutations_mpileup.py

Runtime: ~30-60 min for exome capture data (~50-200 Mb covered regions)

Outputs in /space/s1/sashanikolaeva/Results/clone/somatic_mutations_mpileup/:
    private_variants.tsv
    per_ramet_counts.tsv
    fp_calibration.tsv
    pairwise_somatic_dist.tsv
    somatic_rate_per_chr.png
    chr3_chr4_enrichment.png
"""

import os
import sys
import re
import subprocess
import bisect
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

try:
    from scipy.stats import binomtest as _binomtest
    def binom_test(k, n, p, alternative="greater"):
        return _binomtest(k, n, p, alternative=alternative).pvalue
except ImportError:
    from scipy.stats import binom_test

# ─── Paths ────────────────────────────────────────────────────────────────────
BAM_DIR  = "/space/s1/sashanikolaeva/Results/clone/nquire/bam"
REF      = "/space/s1/sashanikolaeva/Data/SEGI_genome/SEGI.2.0.fa"
NUPT_BED = "/space/s1/sashanikolaeva/Results/clone/nquire/nupt_mask.bed"
OUT_DIR  = "/space/s1/sashanikolaeva/Results/clone/somatic_mutations_mpileup"

# ─── Samples ─────────────────────────────────────────────────────────────────
CLONE_RAMETS = [
    "13_WP1", "6_WP5",    "31_c3",  "31_WP2_1",
    "28_WP7", "28_WP6",   "14_WP3", "22_WPd2",  "33_c6",
]
CLONE_BAMS = {s: f"{BAM_DIR}/{s}.no_chlor.sorted.bam" for s in CLONE_RAMETS}

# Within-pool FP calibration pairs (same PCR pool, same clone genotype)
FP_PAIRS = {
    "pool28": ["28_WP6", "28_WP7"],
    "pool31": ["31_c3",  "31_WP2_1"],
}

CHROMS = [f"chr{i}" for i in range(1, 12)]

# ─── Filters ─────────────────────────────────────────────────────────────────
MIN_DP          = 10    # minimum depth per sample to count as callable
MIN_ALT_READS   = 2     # minimum ALT-supporting reads in the variant carrier
MIN_ALT_FREQ    = 0.10  # minimum ALT allele frequency in variant carrier
                        # (1/6 dosage at DP=12 → 2 reads = 16.7% > 10%)
MAX_ALT_IN_OTHERS = 0   # strict: zero ALT reads allowed in any other sample
MIN_CALLABLE    = 7     # min samples with DP >= MIN_DP per site
MIN_MAPQ        = 20    # minimum mapping quality (passed to mpileup)
MIN_BASEQ       = 20    # minimum base quality (passed to mpileup)
MAX_STRAND_BIAS = 0.95  # if >95% of ALT reads are from one strand → likely artifact

# ─── Helpers ──────────────────────────────────────────────────────────────────

def load_nupt_mask(path):
    """Returns dict chr -> sorted list of (start, end) intervals."""
    mask = defaultdict(list)
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            mask[parts[0]].append((int(parts[1]), int(parts[2])))
    for chrom in mask:
        mask[chrom].sort()
    return mask


def in_nupt(pos, intervals):
    """Binary search: is pos (int) inside any interval?"""
    if not intervals:
        return False
    starts = [iv[0] for iv in intervals]
    idx = bisect.bisect_right(starts, pos) - 1
    return idx >= 0 and intervals[idx][0] <= pos <= intervals[idx][1]


def parse_pileup_bases(bases_str, depth):
    """
    Parse a mpileup bases string for one sample.
    Returns (n_ref, n_alt_dict, n_fwd_alt, n_rev_alt) where:
      n_ref: count of REF alleles
      n_alt_dict: dict {base -> count} for each distinct ALT base (uppercase)
      n_fwd_alt: forward-strand ALT reads (uppercase letters)
      n_rev_alt: reverse-strand ALT reads (lowercase letters)
    """
    # Remove read-start markers (^X where X is any char)
    bases_str = re.sub(r'\^.', '', bases_str)
    # Remove read-end markers
    bases_str = bases_str.replace('$', '')
    # Remove insertion sequences: +NXX...
    bases_str = re.sub(r'[+-](\d+)', lambda m: ' ' * (int(m.group(1)) + len(m.group(0))), bases_str)
    # Remove the spaces we inserted
    bases_str = re.sub(r'\s', '', bases_str)

    n_ref = bases_str.count('.') + bases_str.count(',')
    n_alt_fwd = defaultdict(int)
    n_alt_rev = defaultdict(int)

    for c in bases_str:
        if c in 'ACGT':
            n_alt_fwd[c] += 1
        elif c in 'acgt':
            n_alt_rev[c.upper()] += 1
        # skip *, <, >, N

    n_alt_dict = {}
    for base in set(list(n_alt_fwd.keys()) + list(n_alt_rev.keys())):
        n_alt_dict[base] = n_alt_fwd[base] + n_alt_rev[base]

    n_fwd_alt = sum(n_alt_fwd.values())
    n_rev_alt = sum(n_alt_rev.values())

    return n_ref, n_alt_dict, n_fwd_alt, n_rev_alt


def stream_mpileup(bam_paths, ref, chroms, min_mapq, min_baseq):
    """
    Stream samtools mpileup output for given BAMs.
    Yields (chrom, pos, ref_base, [(depth, bases_str), ...]) per position.
    Number of samples = len(bam_paths).
    BAMs are already restricted to chr1-chr11 (no_chlor filter), so no -r needed.
    """
    cmd = [
        "samtools", "mpileup",
        "-f", ref,
        "-B",                     # disable BAQ recalibration (faster, fine for short targets)
        "-q", str(min_mapq),
        "-Q", str(min_baseq),
        "--no-output-ins",        # skip insertion sequences (+NUM...)
        "--no-output-del",        # skip deletion sequences (-NUM...)
    ] + list(bam_paths)

    err_log = f"{OUT_DIR}/mpileup_stderr.log"
    with open(err_log, "w") as ferr:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=ferr, text=True)
        for line in proc.stdout:
            parts = line.rstrip("\n").split("\t")
            chrom = parts[0]
            pos = int(parts[1])
            ref_base = parts[2].upper()
            n_samples = len(bam_paths)
            sample_data = []
            for i in range(n_samples):
                depth = int(parts[3 + i * 3])
                bases = parts[4 + i * 3]
                sample_data.append((depth, bases))
            yield chrom, pos, ref_base, sample_data
        proc.wait()


# ─── Core private variant caller ─────────────────────────────────────────────

def find_private_variants_mpileup(
    sample_names, nupt_mask,
    min_dp=MIN_DP, min_called=MIN_CALLABLE,
    min_alt_reads=MIN_ALT_READS, min_alt_freq=MIN_ALT_FREQ,
    max_alt_others=MAX_ALT_IN_OTHERS, max_strand_bias=MAX_STRAND_BIAS,
    label="",
):
    """
    Scan BAMs via mpileup and call private variants.

    Returns:
        private_vars: list of dicts
        callable_sites: dict sample_name -> dict chrom -> count
    """
    bam_paths = [CLONE_BAMS[s] for s in sample_names]
    n = len(sample_names)
    # For FP pairs (n=2), lower min_called to match
    effective_min_called = min(min_called, n)

    private_vars = []
    callable_sites = defaultdict(lambda: defaultdict(int))

    n_processed = 0
    n_nupt = 0
    n_low_depth = 0

    print(f"  [{label}] Running mpileup on {n} BAMs...")

    for chrom, pos, ref_base, sample_data in stream_mpileup(
        bam_paths, REF, None, MIN_MAPQ, MIN_BASEQ
    ):
        n_processed += 1
        if n_processed % 100000 == 0:
            print(f"    {n_processed:,} positions, {len(private_vars)} private so far")

        # NUPT filter
        if in_nupt(pos, nupt_mask.get(chrom, [])):
            n_nupt += 1
            continue

        # Depth filter
        callable_mask = [depth >= min_dp for depth, _ in sample_data]
        n_callable = sum(callable_mask)
        if n_callable < effective_min_called:
            n_low_depth += 1
            continue

        # Update callable site counts
        for i, (is_callable, s) in enumerate(zip(callable_mask, sample_names)):
            if is_callable:
                callable_sites[s][chrom] += 1

        # Parse pileup bases per sample
        parsed = []
        for depth, bases_str in sample_data:
            n_ref, n_alt_dict, n_fwd, n_rev = parse_pileup_bases(bases_str, depth)
            parsed.append((depth, n_ref, n_alt_dict, n_fwd, n_rev))

        # Find private variants
        for i in range(n):
            if not callable_mask[i]:
                continue
            depth_i, n_ref_i, n_alt_i, n_fwd_i, n_rev_i = parsed[i]

            # Get the dominant ALT allele for this sample
            if not n_alt_i:
                continue
            alt_base = max(n_alt_i, key=n_alt_i.get)
            n_alt_reads = n_alt_i[alt_base]

            # ALT read count and frequency thresholds
            if n_alt_reads < min_alt_reads:
                continue
            alt_freq = n_alt_reads / depth_i
            if alt_freq < min_alt_freq:
                continue

            # Strand bias filter
            total_alt = n_fwd_i + n_rev_i
            if total_alt > 0:
                strand_frac = max(n_fwd_i, n_rev_i) / total_alt
                if strand_frac > max_strand_bias:
                    continue

            # Check all other callable samples have zero ALT at this base
            others_clean = True
            for j in range(n):
                if j == i or not callable_mask[j]:
                    continue
                _, _, n_alt_j, _, _ = parsed[j]
                if n_alt_j.get(alt_base, 0) > max_alt_others:
                    others_clean = False
                    break

            if others_clean:
                private_vars.append({
                    "label": label,
                    "sample": sample_names[i],
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref_base,
                    "alt": alt_base,
                    "n_alt_reads": n_alt_reads,
                    "depth": depth_i,
                    "alt_freq": round(alt_freq, 3),
                    "n_callable": n_callable,
                })

    print(f"    Done: {n_processed:,} positions | NUPT filtered: {n_nupt:,} | "
          f"Low depth: {n_low_depth:,} | Private variants: {len(private_vars)}")
    return private_vars, callable_sites


# ─── FP calibration ───────────────────────────────────────────────────────────

def compute_fp_rate_mpileup(nupt_mask):
    """
    Run private variant calling on within-pool pairs.
    Pool 28: WP6 vs WP7 — same PCR pool, same clone genotype.
    Pool 31: c3 vs WP2_1 — same PCR pool, same clone genotype.
    Any "private variant" = technical error within that pool.
    """
    fp_rates = {}
    fp_details = []

    for pair_label, pair_samples in FP_PAIRS.items():
        print(f"\n  FP pair [{pair_label}]: {pair_samples}")
        pvars, callable_sites = find_private_variants_mpileup(
            pair_samples, nupt_mask,
            min_called=2,   # both samples must be callable
            label=f"FP_{pair_label}",
        )
        # Mean callable sites across both samples
        total_callable = np.mean([
            sum(callable_sites[s].values()) for s in pair_samples
        ])
        fp_rate = len(pvars) / total_callable if total_callable > 0 else 0
        print(f"    {len(pvars)} FP variants / {total_callable:.0f} callable sites → rate = {fp_rate:.6f}")
        fp_rates[pair_label] = {"fp_rate": fp_rate, "n_fp": len(pvars), "n_callable": total_callable}
        fp_details.extend(pvars)

    mean_fp_rate = np.mean([v["fp_rate"] for v in fp_rates.values()])
    fp_rates["mean"] = {"fp_rate": mean_fp_rate, "n_fp": None, "n_callable": None}
    print(f"\n  Mean FP rate (within-pool): {mean_fp_rate:.6f}")
    return fp_rates, fp_details


# ─── Per-ramet table & figures ────────────────────────────────────────────────

def make_per_ramet_table(private_vars, callable_sites, fp_rate):
    records = []
    for sample in CLONE_RAMETS:
        for chrom in CHROMS:
            n_raw = sum(
                1 for v in private_vars
                if v["sample"] == sample and v["chrom"] == chrom
            )
            n_callable = callable_sites[sample].get(chrom, 0)
            n_fp_exp = fp_rate * n_callable
            n_corrected = max(0.0, n_raw - n_fp_exp)
            records.append({
                "ramet": sample,
                "chrom": chrom,
                "n_raw": n_raw,
                "n_callable": n_callable,
                "n_fp_expected": round(n_fp_exp, 3),
                "n_corrected": round(n_corrected, 3),
                "rate_raw_per_mb": round(n_raw / n_callable * 1e6, 4) if n_callable else 0,
                "rate_corrected_per_mb": round(n_corrected / n_callable * 1e6, 4) if n_callable else 0,
            })
    return pd.DataFrame(records)


def compute_pairwise_dist(private_vars):
    """Pairwise somatic sharing + Jaccard distance matrix."""
    ramet_positions = defaultdict(set)
    for v in private_vars:
        ramet_positions[v["sample"]].add((v["chrom"], v["pos"]))

    n = len(CLONE_RAMETS)
    dist = np.zeros((n, n))
    sharing = np.zeros((n, n), dtype=int)

    for i in range(n):
        for j in range(n):
            pi = ramet_positions[CLONE_RAMETS[i]]
            pj = ramet_positions[CLONE_RAMETS[j]]
            if i == j:
                sharing[i, j] = len(pi)
            else:
                shared = len(pi & pj)
                sharing[i, j] = shared
                union = len(pi | pj)
                dist[i, j] = 1.0 - shared / union if union > 0 else 0.0

    return (pd.DataFrame(sharing, index=CLONE_RAMETS, columns=CLONE_RAMETS),
            pd.DataFrame(dist,    index=CLONE_RAMETS, columns=CLONE_RAMETS))


def plot_somatic_rate_per_chr(df, out_path):
    fig, axes = plt.subplots(3, 4, figsize=(16, 10), sharey=False)
    axes = axes.flatten()
    for idx, chrom in enumerate(CHROMS):
        ax = axes[idx]
        sub = df[df["chrom"] == chrom].sort_values("ramet")
        color = ["#e74c3c" if chrom in ("chr3", "chr4") else "#2980b9"] * len(sub)
        ax.bar(sub["ramet"], sub["rate_corrected_per_mb"], color=color)
        ax.set_title(chrom, fontweight="bold")
        ax.tick_params(axis="x", rotation=45, labelsize=7)
        if idx % 4 == 0:
            ax.set_ylabel("Somatic variants / Mb")
    for idx in range(len(CHROMS), len(axes)):
        axes[idx].set_visible(False)
    fig.suptitle("BAM-based somatic mutation rate per ramet per chromosome\n"
                 "(red = chr3/chr4; FP-corrected)", fontsize=12)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()


def plot_chr3_chr4_enrichment(df, out_path):
    genome_rate = df.groupby("ramet").apply(
        lambda x: x["n_corrected"].sum() / x["n_callable"].sum()
        if x["n_callable"].sum() > 0 else 0,
        include_groups=False
    )
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, focal_chr in zip(axes, ["chr3", "chr4"]):
        focal = df[df["chrom"] == focal_chr].set_index("ramet")
        ramets = sorted(CLONE_RAMETS)
        x = list(range(len(ramets)))
        focal_rates  = [focal.loc[r, "rate_corrected_per_mb"] if r in focal.index else 0 for r in ramets]
        genome_rates = [genome_rate.get(r, 0) * 1e6 for r in ramets]
        ax.bar(x, focal_rates, color="#e74c3c", alpha=0.8, label=focal_chr)
        ax.plot(x, genome_rates, "ko--", markersize=5, linewidth=1, label="genome-wide")
        for i, r in enumerate(ramets):
            if r not in focal.index:
                continue
            n_obs = int(focal.loc[r, "n_corrected"])
            n_cal = int(focal.loc[r, "n_callable"])
            p_exp = genome_rate.get(r, 0)
            if n_cal > 0 and p_exp > 0 and n_obs > 0:
                try:
                    pval = binom_test(n_obs, n_cal, p_exp, alternative="greater")
                    if pval < 0.05:
                        ax.text(i, focal_rates[i] * 1.05, f"p={pval:.3f}",
                                ha="center", fontsize=7, color="red")
                except Exception:
                    pass
        ax.set_xticks(x)
        ax.set_xticklabels(ramets, rotation=45, fontsize=8)
        ax.set_title(f"{focal_chr} enrichment vs genome-wide")
        ax.set_ylabel("Somatic variants / Mb")
        ax.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    nupt_mask = load_nupt_mask(NUPT_BED)
    print(f"NUPT mask: {sum(len(v) for v in nupt_mask.values())} intervals")

    # ── FP calibration ────────────────────────────────────────────────────────
    print("\n=== Step 1: FP Calibration (within-pool pairs) ===")
    fp_rates, fp_details = compute_fp_rate_mpileup(nupt_mask)
    mean_fp_rate = fp_rates["mean"]["fp_rate"]

    fp_rows = [{"pair": k, **v} for k, v in fp_rates.items()]
    pd.DataFrame(fp_rows).to_csv(f"{OUT_DIR}/fp_calibration.tsv", sep="\t", index=False)

    # ── Clone ramet private variants ──────────────────────────────────────────
    print("\n=== Step 2: Clone Ramet Private Variants ===")
    private_vars, callable_sites = find_private_variants_mpileup(
        CLONE_RAMETS, nupt_mask, label="clones"
    )
    df_private = pd.DataFrame(private_vars)
    if not df_private.empty:
        df_private.to_csv(f"{OUT_DIR}/private_variants.tsv", sep="\t", index=False)

    # ── Per-ramet table ────────────────────────────────────────────────────────
    print("\n=== Step 3: Per-Ramet Counts ===")
    df_counts = make_per_ramet_table(private_vars, callable_sites, mean_fp_rate)
    df_counts.to_csv(f"{OUT_DIR}/per_ramet_counts.tsv", sep="\t", index=False)

    total = df_counts.groupby("ramet")[["n_raw", "n_corrected", "n_callable"]].sum()
    print(f"\n{'Ramet':<14} {'Raw':>6} {'Corrected':>10} {'Rate/Mb':>10}")
    for r in total.index:
        raw = total.loc[r, "n_raw"]
        cor = total.loc[r, "n_corrected"]
        cal = total.loc[r, "n_callable"]
        rate = cor / cal * 1e6 if cal > 0 else 0
        print(f"  {r:<12} {raw:>6} {cor:>10.1f} {rate:>10.4f}")

    # ── Pairwise distance ─────────────────────────────────────────────────────
    print("\n=== Step 4: Pairwise Somatic Distance ===")
    if not df_private.empty:
        df_sharing, df_dist = compute_pairwise_dist(private_vars)
        df_sharing.to_csv(f"{OUT_DIR}/pairwise_somatic_sharing.tsv", sep="\t")
        df_dist.to_csv(f"{OUT_DIR}/pairwise_somatic_dist.tsv", sep="\t")
        print("  Somatic sharing matrix (diagonal = total private per ramet):")
        print(df_sharing.to_string())

    # ── Figures ───────────────────────────────────────────────────────────────
    print("\n=== Step 5: Figures ===")
    plot_somatic_rate_per_chr(df_counts, f"{OUT_DIR}/somatic_rate_per_chr.png")
    plot_chr3_chr4_enrichment(df_counts, f"{OUT_DIR}/chr3_chr4_enrichment.png")

    print(f"\nDone. Outputs in {OUT_DIR}/")


if __name__ == "__main__":
    main()
