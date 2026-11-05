#!/usr/bin/env python3
# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Read depth per giant sequoia chromosome from the survey VCF. Its output defines the VCF
# sample set read by scripts/finalize/common.py.
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
"""
Mean read depth per chromosome from VCF — full population.

Uses DP field from VCF at called sites (same SNP map + sample set as het_density_vcf_fullpop.py).
Tests whether clones have elevated depth on chr6/7 (misalignment artifact from chr4 pentasomy?)
or whether depth is flat (suggesting het density signals are real).

Outputs:
  Results/clone/het_density_fullpop/depth_per_chr.tsv
  Results/clone/het_density_fullpop/depth_per_chr_normalized_boxplot.png
"""

import subprocess
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy import stats
from pathlib import Path

VCF     = Path("/space/s1/santang3/prg/results/variants/vcf/allScaffs_allRegions_filtered.vcf.gz")
SNP_MAP = Path("/space/s1/sashanikolaeva/Results/clone/snp_segi_alignment/snp_to_segi_chr.tsv")
OUT_DIR = Path("/space/s1/sashanikolaeva/Results/clone/het_density_fullpop")

CHROMOSOMES = [f"chr{i}" for i in range(1, 12)]

CLONE_SAMPLES = {
    'RedwoodCap_13_WP1', 'RedwoodCap_6_WP5',  'RedwoodCap_31_c3',
    'RedwoodCap_31_WP2_1', 'RedwoodCap_28_WP7', 'RedwoodCap_28_WP6',
    'RedwoodCap_14_WP3', 'RedwoodCap_22_WPd2', 'RedwoodCap_33_c6',
}

SAMPLES_TO_REMOVE = {
    'RedwoodCap_35_150', 'RedwoodCap_36_JSM4',
    'RedwoodCap_22_MRD', 'RedwoodCap_13_286_1', 'RedwoodCap_16_416-1',
    'RedwoodCap_17_111_2', 'RedwoodCap_36_M_spl', 'RedwoodCap_8_374',
    'RedwoodCap_14_86', 'RedwoodCap_26_46H_sg',
    'RedwoodCap_24_57_4', 'RedwoodCap_12_57', 'RedwoodCap_23_57-3',
    'RedwoodCap_10_NQ2', 'RedwoodCap_27_NQH1',
    'RedwoodCap_35_97_extra', 'RedwoodCap_10_368-2', 'RedwoodCap_21_63_5',
    'RedwoodCap_21_45-2', 'RedwoodCap_19_382_3', 'RedwoodCap_30_157-0X',
    'RedwoodCap_1_263_1', 'RedwoodCap_8_63-8', 'RedwoodCap_4_BH4',
    'RedwoodCap_15_BH2', 'RedwoodCap_36_NP7', 'RedwoodCap_20_RW3',
    'RedwoodCap_29_RW1', 'RedwoodCap_30_RW5', 'RedwoodCap_36_197',
    'RedwoodCap_33_256', 'RedwoodCap_23_FR2', 'RedwoodCap_17_380',
    'RedwoodCap_17_380-2', 'RedwoodCap_24_380-3', 'RedwoodCap_20_365',
    'RedwoodCap_16_257-1', 'RedwoodCap_13_340-1', 'RedwoodCap_30_PG2',
    'RedwoodCap_23_182', 'RedwoodCap_30_500', 'RedwoodCap_28_95',
    'RedwoodCap_29_fl3', 'RedwoodCap_23_17', 'RedwoodCap_31_J5',
    'RedwoodCap_32_J4',
}


def load_snp_map(path):
    snp_to_chr = {}
    chr_set = set(CHROMOSOMES)
    with open(path) as f:
        next(f)
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            scaffold, pos, segi_chr = parts[0], parts[1], parts[2]
            if segi_chr in chr_set:
                snp_to_chr[(scaffold, int(pos))] = segi_chr
    return snp_to_chr


def main():
    print("Reading VCF sample list ...")
    result = subprocess.run(['bcftools', 'query', '-l', str(VCF)],
                            capture_output=True, text=True, check=True)
    all_vcf_samples = result.stdout.strip().split('\n')

    kept_samples     = [s for s in all_vcf_samples if s not in SAMPLES_TO_REMOVE]
    clone_samples    = [s for s in kept_samples if s in CLONE_SAMPLES]
    nonclone_samples = [s for s in kept_samples if s not in CLONE_SAMPLES]
    print(f"  Clone: {len(clone_samples)}  Non-clone: {len(nonclone_samples)}  Total: {len(kept_samples)}")
    sys.stdout.flush()

    print(f"Loading SNP map ...")
    snp_to_chr = load_snp_map(SNP_MAP)
    print(f"  {len(snp_to_chr):,} SNPs mapped to chr1-chr11")
    sys.stdout.flush()

    # Accumulate DP values per sample per chromosome
    dp_sum   = {s: {c: 0.0 for c in CHROMOSOMES} for s in kept_samples}
    dp_count = {s: {c: 0   for c in CHROMOSOMES} for s in kept_samples}

    print(f"Streaming VCF ...")
    sys.stdout.flush()

    samples_arg = ','.join(kept_samples)
    cmd = ['bcftools', 'view', '-s', samples_arg, str(VCF)]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)

    sample_indices = None
    n_processed = 0

    for line in proc.stdout:
        if line.startswith('##'):
            continue
        if line.startswith('#CHROM'):
            fields = line.rstrip('\n').split('\t')
            sample_indices = {}
            for i, col in enumerate(fields[9:], start=9):
                if col in kept_samples:
                    sample_indices[i] = col
            print(f"  Found {len(sample_indices)} samples in header")
            sys.stdout.flush()
            continue

        fields   = line.rstrip('\n').split('\t')
        scaffold = fields[0]
        pos      = int(fields[1])

        chrom = snp_to_chr.get((scaffold, pos))
        if chrom is None:
            continue

        n_processed += 1
        fmt    = fields[8].split(':')
        dp_idx = fmt.index('DP') if 'DP' in fmt else None
        if dp_idx is None:
            continue

        for col_idx, sample in sample_indices.items():
            sample_fields = fields[col_idx].split(':')
            if dp_idx >= len(sample_fields):
                continue
            dp_str = sample_fields[dp_idx]
            if dp_str == '.':
                continue
            dp = int(dp_str)
            if dp < 10:
                continue
            dp_sum[sample][chrom]   += dp
            dp_count[sample][chrom] += 1

        if n_processed % 50000 == 0:
            print(f"  {n_processed:,} SNPs processed ...")
            sys.stdout.flush()

    proc.wait()
    print(f"Done. {n_processed:,} mapped SNPs processed.")

    # Mean DP per sample per chromosome
    mean_dp = {}
    for s in kept_samples:
        mean_dp[s] = {}
        for c in CHROMOSOMES:
            n = dp_count[s][c]
            mean_dp[s][c] = dp_sum[s][c] / n if n > 0 else float('nan')

    # Save TSV
    out_tsv = OUT_DIR / "depth_per_chr.tsv"
    with open(out_tsv, 'w') as f:
        f.write("sample\tgroup\tchrom\tmean_dp\tn_sites\n")
        for s in kept_samples:
            group = "clone" if s in CLONE_SAMPLES else "non-clone"
            for c in CHROMOSOMES:
                f.write(f"{s}\t{group}\t{c}\t{mean_dp[s][c]:.3f}\t{dp_count[s][c]}\n")
    print(f"Saved: {out_tsv}")

    # Within-sample normalization (divide each chr by sample genome-wide mean)
    clone_norm    = {c: [] for c in CHROMOSOMES}
    nonclone_norm = {c: [] for c in CHROMOSOMES}
    for s in kept_samples:
        vals = [mean_dp[s][c] for c in CHROMOSOMES]
        valid = [v for v in vals if not np.isnan(v)]
        if not valid:
            continue
        gm = np.mean(valid)
        for c in CHROMOSOMES:
            if not np.isnan(mean_dp[s][c]):
                norm = mean_dp[s][c] / gm
                if s in CLONE_SAMPLES:
                    clone_norm[c].append(norm)
                else:
                    nonclone_norm[c].append(norm)

    # Print stats
    print(f"\nNormalized mean depth (per-sample / genome-wide mean):")
    print(f"{'chr':>6}  {'clone':>8}  {'non-clone':>10}  {'diff':>8}  {'p':>8}")
    print("-" * 52)
    norm_pvals = []
    for c in CHROMOSOMES:
        cm = np.mean(clone_norm[c])
        nm = np.mean(nonclone_norm[c])
        _, p = stats.mannwhitneyu(clone_norm[c], nonclone_norm[c], alternative='two-sided')
        norm_pvals.append(p)
        sig = '*' if p < 0.05 else ''
        print(f"{c:>6}  {cm:.4f}    {nm:.4f}      {cm-nm:+.4f}  {p:>8.4f} {sig}")

    nc = len(clone_samples)
    nn = len(nonclone_samples)

    # Plot: normalized depth boxplot
    x = np.arange(len(CHROMOSOMES))
    fig, ax = plt.subplots(figsize=(14, 6))

    positions_c = x - 0.2
    positions_n = x + 0.2

    bp_c = ax.boxplot([clone_norm[c]    for c in CHROMOSOMES], positions=positions_c, widths=0.3,
                      patch_artist=True, showfliers=True, medianprops=dict(color='black'))
    bp_n = ax.boxplot([nonclone_norm[c] for c in CHROMOSOMES], positions=positions_n, widths=0.3,
                      patch_artist=True, showfliers=True, medianprops=dict(color='black'))

    for patch in bp_c['boxes']:
        patch.set_facecolor('#d62728'); patch.set_alpha(0.6)
    for patch in bp_n['boxes']:
        patch.set_facecolor('#1f77b4'); patch.set_alpha(0.6)

    ax.axhline(1.0, color='gray', ls='--', lw=1, alpha=0.5)

    for i, c in enumerate(CHROMOSOMES):
        if norm_pvals[i] < 0.05:
            y_max = max(max(clone_norm[c]) if clone_norm[c] else 1,
                        max(nonclone_norm[c]) if nonclone_norm[c] else 1)
            ax.text(i, y_max + 0.01, '*', ha='center', fontsize=14, fontweight='bold')

    ax.set_xticks(x); ax.set_xticklabels(CHROMOSOMES, fontsize=10)
    ax.set_xlabel("Chromosome (SEGI)", fontsize=12)
    ax.set_ylabel("Normalized mean DP\n(relative to sample genome-wide mean)", fontsize=11)
    ax.set_title(
        f"Mean read depth per chromosome: Clones (red, n={nc}) vs Non-clones (blue, n={nn})\n"
        f"VCF DP field at SNP sites; depth ≥ 10; * = p<0.05 (Mann-Whitney)\n"
        f"Test: do chr6/7 depth elevations mirror chr4 pentasomy (misalignment artifact)?",
        fontsize=11, fontweight="bold"
    )
    ax.legend(handles=[Patch(facecolor='#d62728', alpha=0.6, label=f'Clones (n={nc})'),
                        Patch(facecolor='#1f77b4', alpha=0.6, label=f'Non-clones (n={nn})')],
              fontsize=10)

    plt.tight_layout()
    out_fig = OUT_DIR / "depth_per_chr_normalized_boxplot.png"
    plt.savefig(out_fig, dpi=200, bbox_inches='tight')
    print(f"Saved: {out_fig}")
    plt.close()

    print(f"\nDone. Results in {OUT_DIR}")


if __name__ == "__main__":
    main()
