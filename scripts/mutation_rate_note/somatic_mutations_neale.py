#!/usr/bin/env python3
# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Stem-private search on the coast redwood assembly (Neale et al. 2022). Returned no sites
# because of the scaffold-name mismatch (Supplementary Note S1).
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
"""
Somatic mutation calling for coast redwood clone ramets — Neale SESE 2.1 reference.

Steps:
  1. FP calibration: mpileup on true technical replicate pairs (MRD, 416)
  2. Private variant calling: mpileup on 9 clone ramets
  3. FP correction and per-ramet counts
  4. Self-calibrated mutation rate (ramet age 165–177 years)

Quality filters match FreeBayes VCF parameters exactly:
  -q 30  (--min-mapping-quality 30)
  -Q 20  (--min-base-quality 20)

Position filter: exome target BED (passed via mpileup -l), replaces SEGI NUPT mask.

Usage:
    conda run -n ntsynt-310 python3 somatic_mutations_neale.py

Runtime: ~1–2 hours (exome capture data, ~152 Mb target on Neale ref)

Outputs in /space/s1/sashanikolaeva/Results/clone/somatic_mutations_neale/:
    fp_calibration.tsv
    private_variants.tsv
    per_ramet_counts.tsv
    mutation_rate_summary.tsv
    somatic_rate_summary.png
"""

import os
import re
import subprocess
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ─── Paths ────────────────────────────────────────────────────────────────────
BAM_DIR  = "/space/s1/sashanikolaeva/Results/clone/somatic_mutations_neale/bams"
REF      = "/space/s1/santang3/prg/results/ref/neale/neale.fa"
EXOME_BED = "/space/s1/santang3/prg/resources/data/Annotation/SESE.2.1.CDS.merge.bed"
OUT_DIR  = "/space/s1/sashanikolaeva/Results/clone/somatic_mutations_neale"

# ─── Samples ──────────────────────────────────────────────────────────────────
CLONE_RAMETS = [
    "13_WP1", "6_WP5", "31_c3", "31_WP2_1",
    "28_WP7", "28_WP6", "14_WP3", "22_WPd2", "33_c6",
]

FP_PAIRS = {
    "MRD": ["22_MRD",   "28_MRD"],
    "416": ["16_416-1", "35_416_1"],
}

def bam(sample):
    return f"{BAM_DIR}/RedwoodCap_{sample}_neale_sorted_dupsMarked.bam"

# ─── Filters (match FreeBayes VCF parameters) ─────────────────────────────────
MIN_DP          = 10
MIN_ALT_READS   = 2
MIN_ALT_FREQ    = 0.08   # 1/6 at DP=12 → 2/12 = 16.7%; lower bound allows DP=25
MAX_ALT_FREQ    = 0.30   # upper bound: >2 copies would be population variant
MAX_ALT_OTHERS  = 0      # strict: zero ALT reads in any other callable sample
MIN_CALLABLE    = 7      # min ramets with DP >= MIN_DP per site (for 9-sample call)
MIN_MAPQ        = 30     # matches FreeBayes --min-mapping-quality 30
MIN_BASEQ       = 20     # matches FreeBayes --min-base-quality 20
MAX_STRAND_BIAS = 0.95   # >95% reads from one strand → artifact

# Ramet age estimates from logging history + 1948 aerial photography
RAMET_AGE_CENTRAL = 177  # years (2026 - 1849 logging end)
RAMET_AGE_LOW     = 165  # conservative (some sprouting delay after logging)

# ─── Pileup parsing ────────────────────────────────────────────────────────────

def parse_pileup_bases(bases_str):
    bases_str = re.sub(r'\^.', '', bases_str)
    bases_str = bases_str.replace('$', '')
    bases_str = re.sub(r'[+-](\d+)', lambda m: ' ' * (int(m.group(1)) + len(m.group(0))), bases_str)
    bases_str = re.sub(r'\s', '', bases_str)

    n_alt_fwd = defaultdict(int)
    n_alt_rev = defaultdict(int)

    for c in bases_str:
        if c in 'ACGT':
            n_alt_fwd[c] += 1
        elif c in 'acgt':
            n_alt_rev[c.upper()] += 1

    n_alt_dict = {}
    for base in set(list(n_alt_fwd) + list(n_alt_rev)):
        n_alt_dict[base] = n_alt_fwd[base] + n_alt_rev[base]

    n_fwd_alt = sum(n_alt_fwd.values())
    n_rev_alt = sum(n_alt_rev.values())
    return n_alt_dict, n_fwd_alt, n_rev_alt


def stream_mpileup(bam_paths, err_log):
    cmd = [
        "samtools", "mpileup",
        "-f", REF,
        "-l", EXOME_BED,
        "-B",
        "-q", str(MIN_MAPQ),
        "-Q", str(MIN_BASEQ),
        "--no-output-ins",
        "--no-output-del",
    ] + list(bam_paths)

    with open(err_log, "a") as ferr:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=ferr, text=True)
        for line in proc.stdout:
            parts = line.rstrip("\n").split("\t")
            chrom = parts[0]
            pos   = int(parts[1])
            n_samp = len(bam_paths)
            sample_data = []
            for i in range(n_samp):
                depth = int(parts[3 + i * 3])
                bases = parts[4 + i * 3]
                sample_data.append((depth, bases))
            yield chrom, pos, sample_data
        proc.wait()


# ─── FP calibration ───────────────────────────────────────────────────────────

def calibrate_fp():
    print("\n=== FP Calibration (true technical replicates) ===")
    fp_results = {}
    err_log = f"{OUT_DIR}/fp_mpileup_stderr.log"

    for pair_label, samples in FP_PAIRS.items():
        bam_paths = [bam(s) for s in samples]
        missing = [p for p in bam_paths if not os.path.exists(p)]
        if missing:
            print(f"  [{pair_label}] SKIP — BAM not found: {missing}")
            continue

        print(f"  [{pair_label}] {samples[0]} vs {samples[1]}")
        n_fp = 0
        n_callable = 0
        n_processed = 0

        for chrom, pos, sample_data in stream_mpileup(bam_paths, err_log):
            n_processed += 1
            if n_processed % 500_000 == 0:
                print(f"    {n_processed:,} pos | {n_fp} FP | {n_callable:,} callable")

            depths = [sd[0] for sd in sample_data]
            if any(d < MIN_DP for d in depths):
                continue
            n_callable += 1

            parsed = [parse_pileup_bases(sd[1]) for sd in sample_data]

            for i in range(2):
                j = 1 - i
                n_alt_i, n_fwd_i, n_rev_i = parsed[i]
                n_alt_j, _, _ = parsed[j]
                depth_i = sample_data[i][0]

                if not n_alt_i:
                    continue
                alt_base = max(n_alt_i, key=n_alt_i.get)
                n_alt_reads = n_alt_i[alt_base]

                if n_alt_reads < MIN_ALT_READS:
                    continue
                af = n_alt_reads / depth_i
                if af < MIN_ALT_FREQ or af > MAX_ALT_FREQ:
                    continue

                total_alt = n_fwd_i + n_rev_i
                if total_alt > 0 and max(n_fwd_i, n_rev_i) / total_alt > MAX_STRAND_BIAS:
                    continue

                if n_alt_j.get(alt_base, 0) > MAX_ALT_OTHERS:
                    continue

                n_fp += 1

        fp_rate = n_fp / n_callable if n_callable > 0 else 0
        fp_results[pair_label] = {
            "pair": f"{samples[0]} vs {samples[1]}",
            "n_fp": n_fp,
            "n_callable": n_callable,
            "fp_rate": fp_rate,
        }
        print(f"    {n_processed:,} pos | callable: {n_callable:,} | FP: {n_fp} | "
              f"rate: {fp_rate:.2e} ({fp_rate*1e6:.2f}/Mb)")

    if not fp_results:
        raise RuntimeError("No FP pairs completed — check BAM paths")

    mean_fp_rate = float(np.mean([v["fp_rate"] for v in fp_results.values()]))
    fp_results["mean"] = {"pair": "mean", "n_fp": None, "n_callable": None, "fp_rate": mean_fp_rate}
    print(f"\n  Mean FP rate: {mean_fp_rate:.2e} ({mean_fp_rate*1e6:.2f}/Mb)")

    pd.DataFrame(fp_results).T.reset_index(drop=True).to_csv(
        f"{OUT_DIR}/fp_calibration.tsv", sep="\t", index=False
    )
    return mean_fp_rate


# ─── Private variant calling ──────────────────────────────────────────────────

def call_private_variants():
    print("\n=== Private Variant Calling (9 clone ramets) ===")
    bam_paths = [bam(s) for s in CLONE_RAMETS]
    missing = [p for p, s in zip(bam_paths, CLONE_RAMETS) if not os.path.exists(p)]
    if missing:
        raise FileNotFoundError(f"BAMs not found: {missing}")

    err_log = f"{OUT_DIR}/mpileup_stderr.log"
    private_vars = []
    callable_sites = defaultdict(lambda: defaultdict(int))  # sample -> scaffold -> count
    n_processed = 0

    for chrom, pos, sample_data in stream_mpileup(bam_paths, err_log):
        n_processed += 1
        if n_processed % 500_000 == 0:
            print(f"  {n_processed:,} pos | {len(private_vars)} private so far")

        callable_mask = [sd[0] >= MIN_DP for sd in sample_data]
        if sum(callable_mask) < MIN_CALLABLE:
            continue

        for i, s in enumerate(CLONE_RAMETS):
            if callable_mask[i]:
                callable_sites[s][chrom] += 1

        parsed = [parse_pileup_bases(sd[1]) for sd in sample_data]

        for i, s in enumerate(CLONE_RAMETS):
            if not callable_mask[i]:
                continue
            n_alt_i, n_fwd_i, n_rev_i = parsed[i]
            depth_i = sample_data[i][0]

            if not n_alt_i:
                continue
            alt_base = max(n_alt_i, key=n_alt_i.get)
            n_alt_reads = n_alt_i[alt_base]

            if n_alt_reads < MIN_ALT_READS:
                continue
            af = n_alt_reads / depth_i
            if af < MIN_ALT_FREQ or af > MAX_ALT_FREQ:
                continue

            total_alt = n_fwd_i + n_rev_i
            if total_alt > 0 and max(n_fwd_i, n_rev_i) / total_alt > MAX_STRAND_BIAS:
                continue

            # All other callable ramets must have zero ALT reads
            others_clean = all(
                parsed[j][0].get(alt_base, 0) <= MAX_ALT_OTHERS
                for j in range(len(CLONE_RAMETS))
                if j != i and callable_mask[j]
            )
            if not others_clean:
                continue

            private_vars.append({
                "sample":      s,
                "chrom":       chrom,
                "pos":         pos,
                "alt":         alt_base,
                "n_alt_reads": n_alt_reads,
                "depth":       depth_i,
                "alt_freq":    round(af, 3),
                "n_callable":  sum(callable_mask),
            })

    print(f"  Done: {n_processed:,} positions | {len(private_vars)} private variants")
    return private_vars, callable_sites


# ─── Shared variant calling ───────────────────────────────────────────────────

def call_shared_variants():
    """
    Find sites where ALL 9 clone ramets carry the same ALT allele at ~1/6 dosage.
    These are candidates for somatic mutations arising in the founding tree.

    The VCF-based approach misses these because batch variant calling across 284
    samples filters out alleles at ~0.5% population AF. BAM-level detection is
    required.

    Returns list of candidate dicts with (chrom, pos, alt, mean_af_across_ramets).
    A subsequent permutation step or VCF AF lookup is needed to estimate FP rate.
    """
    print("\n=== Shared Variant Calling (all 9 clone ramets carry ALT) ===")
    bam_paths = [bam(s) for s in CLONE_RAMETS]
    missing = [p for p, s in zip(bam_paths, CLONE_RAMETS) if not os.path.exists(p)]
    if missing:
        raise FileNotFoundError(f"BAMs not found: {missing}")

    err_log = f"{OUT_DIR}/shared_mpileup_stderr.log"
    shared_vars = []
    n_processed = 0
    n_all_callable = 0

    for chrom, pos, sample_data in stream_mpileup(bam_paths, err_log):
        n_processed += 1
        if n_processed % 500_000 == 0:
            print(f"  {n_processed:,} pos | {len(shared_vars)} shared so far")

        # All 9 must be callable
        if any(sd[0] < MIN_DP for sd in sample_data):
            continue
        n_all_callable += 1

        parsed = [parse_pileup_bases(sd[1]) for sd in sample_data]

        # Candidate ALT bases: union of ALTs seen across all samples
        all_alts = set()
        for n_alt_i, _, _ in parsed:
            all_alts.update(n_alt_i.keys())

        for alt_base in all_alts:
            # Every ramet must carry this ALT at ~1/6 freq
            all_carry = True
            ao_list = []
            af_list = []
            for i, (n_alt_i, n_fwd_i, n_rev_i) in enumerate(parsed):
                depth_i = sample_data[i][0]
                n_alt_reads = n_alt_i.get(alt_base, 0)
                af = n_alt_reads / depth_i if depth_i > 0 else 0

                if n_alt_reads < MIN_ALT_READS or af < MIN_ALT_FREQ or af > MAX_ALT_FREQ:
                    all_carry = False
                    break

                # Strand bias per ramet
                total_alt = n_fwd_i + n_rev_i
                if total_alt > 0 and max(n_fwd_i, n_rev_i) / total_alt > MAX_STRAND_BIAS:
                    all_carry = False
                    break

                ao_list.append(n_alt_reads)
                af_list.append(af)

            if not all_carry or len(ao_list) < 9:
                continue

            mean_af = float(np.mean(af_list))
            shared_vars.append({
                "chrom":       chrom,
                "pos":         pos,
                "alt":         alt_base,
                "mean_af":     round(mean_af, 3),
                "min_ao":      min(ao_list),
                "n_callable_all9": n_all_callable,
            })

    print(f"  Done: {n_processed:,} positions | {n_all_callable:,} all-callable | "
          f"{len(shared_vars)} shared variant candidates")
    return shared_vars, n_all_callable


def annotate_shared_with_vcf_af(shared_vars):
    """
    For each shared candidate, look up its population AF in the VCF INFO field.
    Sites absent from the VCF (AF not found) are stronger somatic candidates.
    Sites in the VCF with low AF may still be population variants.
    """
    if not shared_vars:
        return []

    print(f"\n  Cross-checking {len(shared_vars)} shared candidates against VCF...")
    vcf = "/space/s1/santang3/prg/results/variants/vcf/allScaffs_allRegions_filtered.vcf.gz"

    # Build a lookup: (chrom, pos) -> AF from VCF
    # Use bcftools query for speed
    positions = "\n".join(f"{v['chrom']}\t{v['pos']}" for v in shared_vars)
    import tempfile, subprocess
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write(positions)
        pos_file = f.name

    cmd = ["bcftools", "query", "-R", pos_file, "-f", "%CHROM\t%POS\t%INFO/AF\n", vcf]
    result = subprocess.run(cmd, capture_output=True, text=True)
    vcf_af = {}
    for line in result.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) >= 3:
            try:
                vcf_af[(parts[0], int(parts[1]))] = float(parts[2])
            except ValueError:
                pass
    os.unlink(pos_file)

    for v in shared_vars:
        key = (v["chrom"], v["pos"])
        v["vcf_af"] = vcf_af.get(key, None)   # None = not in VCF = stronger somatic candidate
        v["in_vcf"] = key in vcf_af

    n_not_in_vcf = sum(1 for v in shared_vars if not v["in_vcf"])
    print(f"  {n_not_in_vcf}/{len(shared_vars)} candidates not found in population VCF "
          f"(stronger somatic evidence)")
    return shared_vars


# ─── Per-ramet table ──────────────────────────────────────────────────────────

def make_per_ramet_table(private_vars, callable_sites, fp_rate):
    records = []
    for s in CLONE_RAMETS:
        n_raw      = sum(1 for v in private_vars if v["sample"] == s)
        n_callable = sum(callable_sites[s].values())
        n_fp_exp   = fp_rate * n_callable
        n_corrected = max(0.0, n_raw - n_fp_exp)
        records.append({
            "ramet":       s,
            "n_raw":       n_raw,
            "n_callable":  n_callable,
            "n_fp_expected": round(n_fp_exp, 2),
            "n_corrected": round(n_corrected, 2),
            "rate_raw_per_mb":  round(n_raw / n_callable * 1e6, 4) if n_callable else 0,
            "rate_corr_per_mb": round(n_corrected / n_callable * 1e6, 4) if n_callable else 0,
        })
    return pd.DataFrame(records)


# ─── Mutation rate calibration ────────────────────────────────────────────────

def compute_mutation_rate(df_counts, fp_rate):
    """
    Self-calibrated somatic mutation rate using known ramet age.

    μ = mean(n_corrected) / (ramet_age × callable_bases × ploidy)

    Ploidy = 6 (hexaploid) — somatic mutations arise in one of 6 haploid copies.
    """
    # Exclude c3 and WP2_1 if they are outliers (same pool 31, may share errors)
    # Use all 9 ramets by default; flag outliers separately
    PLOIDY = 6

    # Use median callable bases (robust to one ramet with unusually low coverage)
    callable_all = df_counts["n_callable"].values
    mean_callable = float(np.mean(callable_all))
    n_corrected_all = df_counts["n_corrected"].values

    records = []
    for age in [RAMET_AGE_LOW, RAMET_AGE_CENTRAL]:
        for excl_pool31 in [False, True]:
            mask = df_counts["ramet"].str.contains("31_") if excl_pool31 else [True] * len(df_counts)
            subset = df_counts[~mask] if excl_pool31 else df_counts
            mean_n = float(np.mean(subset["n_corrected"].values))
            mean_cal = float(np.mean(subset["n_callable"].values))
            mu = mean_n / (age * mean_cal * PLOIDY) if mean_cal > 0 else 0
            records.append({
                "ramet_age_yr": age,
                "exclude_pool31": excl_pool31,
                "n_ramets": len(subset),
                "mean_n_corrected": round(mean_n, 2),
                "mean_callable":    int(mean_cal),
                "mu_per_site_per_yr": f"{mu:.3e}",
            })
    return pd.DataFrame(records)


# ─── Figures ──────────────────────────────────────────────────────────────────

def plot_summary(df_counts, out_path):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    x = np.arange(len(CLONE_RAMETS))
    labels = [r.replace("_", "\n") for r in CLONE_RAMETS]
    pool31 = [r.startswith("31_") for r in CLONE_RAMETS]

    ax = axes[0]
    colors = ["#e74c3c" if p else "#2980b9" for p in pool31]
    ax.bar(x - 0.2, df_counts["n_raw"],       0.4, label="Raw",       color="#95a5a6")
    ax.bar(x + 0.2, df_counts["n_corrected"], 0.4, label="Corrected", color=colors)
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("Private variants")
    ax.set_title("Raw vs FP-corrected private variants\n(red = pool31 outlier)")
    ax.legend()

    ax = axes[1]
    ax.bar(x, df_counts["rate_corr_per_mb"], color=colors)
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("Corrected variants / Mb")
    ax.set_title("Per-ramet somatic mutation rate\n(Neale reference, exome-callable positions)")

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out_path}")


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    # Step 1: FP calibration
    fp_rate = calibrate_fp()

    # Step 2: Private variant calling
    private_vars, callable_sites = call_private_variants()

    pd.DataFrame(private_vars).to_csv(f"{OUT_DIR}/private_variants.tsv", sep="\t", index=False)

    # Step 3: Per-ramet counts
    df_counts = make_per_ramet_table(private_vars, callable_sites, fp_rate)
    df_counts.to_csv(f"{OUT_DIR}/per_ramet_counts.tsv", sep="\t", index=False)

    print("\n=== Per-ramet summary ===")
    print(f"{'Ramet':<14} {'Raw':>6} {'FP_exp':>8} {'Corrected':>10} {'Callable(Mb)':>13} {'Rate/Mb':>9}")
    for _, row in df_counts.iterrows():
        print(f"  {row['ramet']:<12} {row['n_raw']:>6} {row['n_fp_expected']:>8.1f} "
              f"{row['n_corrected']:>10.1f} {row['n_callable']/1e6:>13.1f} {row['rate_corr_per_mb']:>9.4f}")

    # Step 4: Mutation rate calibration
    df_rate = compute_mutation_rate(df_counts, fp_rate)
    df_rate.to_csv(f"{OUT_DIR}/mutation_rate_summary.tsv", sep="\t", index=False)

    print("\n=== Self-calibrated somatic mutation rate ===")
    print(df_rate.to_string(index=False))

    print(f"\n  Published comparison:")
    print(f"  Moore 2023 (redwood): 6.9×10⁻⁹ to 3.4×10⁻⁸/site/yr")
    print(f"  Hanlon 2019 (Picea):  1–2×10⁻⁹/site/yr")

    # Step 5: Shared variant calling (founding tree age)
    print("\n=== Step 5: Shared Variant Calling ===")
    shared_vars, n_all_callable = call_shared_variants()
    shared_vars = annotate_shared_with_vcf_af(shared_vars)

    df_shared = pd.DataFrame(shared_vars)
    df_shared.to_csv(f"{OUT_DIR}/shared_variants_raw.tsv", sep="\t", index=False)

    print(f"\n  Raw shared candidates: {len(df_shared)}")
    if not df_shared.empty:
        not_in_vcf = df_shared[df_shared["in_vcf"] == False]
        print(f"  Not in population VCF: {len(not_in_vcf)} (strongest somatic candidates)")
        print(f"  Run shared_somatic_variants.py --permutation-mode for FP correction")

    # Figures
    plot_summary(df_counts, f"{OUT_DIR}/somatic_rate_summary.png")

    print(f"\nDone. Outputs in {OUT_DIR}/")


if __name__ == "__main__":
    main()
