#!/usr/bin/env python3
# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Heterozygous-site density per chromosome from VCF allele depths. Read by
# scripts/finalize/02_R1ab.py as a secondary column; not reported in the paper.
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
"""
Het density per chromosome from VCF — full population, AD-based het calling.

Replaces het_density_vcf_fullpop.py which used GT field (FreeBayes hexaploid model).
GT-based calling is ploidy-dependent: for pentasomic chr4 in clones (5 copies),
FreeBayes GT calls are misspecified (model assumes 6 copies). This is not a fair
comparison with hexaploid non-clones.

This script uses the AD (allelic depth) field instead:
  ref_dp = AD[0]
  alt_dp = sum(AD[1:])   (handles multiallelic sites)
  total  = ref_dp + alt_dp

  A site is heterozygous if:
    total >= MIN_DEPTH (10)         — adequate coverage
    min(ref_dp, alt_dp) >= MIN_ALT (3)  — both alleles genuinely present

MIN_ALT=3 safely captures 1/6 dosage at ~30x (~5 alt reads expected) while
filtering sequencing errors (~1% error rate = ~0.3 alt reads at 30x).

This approach is completely ploidy-agnostic: works equally for 5-copy or 6-copy
chromosomes. Het density = n_het / n_total per sample per chromosome.
Within-sample normalization and Mann-Whitney test as before.

Outputs:
  Results/clone/het_density_fullpop/het_density_ad_fullpop.tsv
  Results/clone/het_density_fullpop/het_density_ad_barplot.png
  Results/clone/het_density_fullpop/het_density_ad_normalized_boxplot.png
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

MIN_DEPTH = 10   # minimum total read depth at site
MIN_ALT   = 3    # minimum alt reads to call het (filters sequencing errors)

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
    all_vcf_samples  = result.stdout.strip().split('\n')
    kept_samples     = [s for s in all_vcf_samples if s not in SAMPLES_TO_REMOVE]
    clone_samples    = [s for s in kept_samples if s in CLONE_SAMPLES]
    nonclone_samples = [s for s in kept_samples if s not in CLONE_SAMPLES]
    print(f"  Clone: {len(clone_samples)}  Non-clone: {len(nonclone_samples)}  Total: {len(kept_samples)}")
    sys.stdout.flush()

    print("Loading SNP map ...")
    snp_to_chr = load_snp_map(SNP_MAP)
    print(f"  {len(snp_to_chr):,} SNPs mapped to chr1-chr11")
    sys.stdout.flush()

    n_het   = {s: {c: 0 for c in CHROMOSOMES} for s in kept_samples}
    n_total = {s: {c: 0 for c in CHROMOSOMES} for s in kept_samples}

    print(f"\nStreaming VCF ({len(kept_samples)} samples) ...")
    sys.stdout.flush()

    samples_arg = ','.join(kept_samples)
    cmd  = ['bcftools', 'view', '-s', samples_arg, str(VCF)]
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
        ad_idx = fmt.index('AD') if 'AD' in fmt else None
        if ad_idx is None:
            continue

        for col_idx, sample in sample_indices.items():
            sf = fields[col_idx].split(':')
            if ad_idx >= len(sf) or sf[ad_idx] == '.':
                continue

            ad_parts = sf[ad_idx].split(',')
            if len(ad_parts) < 2:
                continue

            try:
                ref_dp = int(ad_parts[0])
                alt_dp = sum(int(x) for x in ad_parts[1:])
            except ValueError:
                continue

            total = ref_dp + alt_dp
            if total < MIN_DEPTH:
                continue

            n_total[sample][chrom] += 1
            if min(ref_dp, alt_dp) >= MIN_ALT:
                n_het[sample][chrom] += 1

        if n_processed % 50000 == 0:
            print(f"  {n_processed:,} SNPs processed ...")
            sys.stdout.flush()

    proc.wait()
    print(f"\nDone. {n_processed:,} mapped SNPs processed.")

    # ── Het rate per sample per chromosome ──
    het_rate = {}
    for s in kept_samples:
        het_rate[s] = {}
        for c in CHROMOSOMES:
            nt = n_total[s][c]
            het_rate[s][c] = n_het[s][c] / nt if nt > 0 else float('nan')

    # ── Save TSV ──
    out_tsv = OUT_DIR / "het_density_ad_fullpop.tsv"
    print(f"\nSaving TSV to {out_tsv} ...")
    with open(out_tsv, 'w') as f:
        f.write("sample\tgroup\tchrom\tn_het\tn_total\thet_rate\n")
        for s in kept_samples:
            group = "clone" if s in CLONE_SAMPLES else "non-clone"
            for c in CHROMOSOMES:
                f.write(f"{s}\t{group}\t{c}\t{n_het[s][c]}\t{n_total[s][c]}\t{het_rate[s][c]:.6f}\n")

    # ── Per-chromosome stats (absolute) ──
    clone_rates    = {c: [het_rate[s][c] for s in clone_samples    if not np.isnan(het_rate[s][c])] for c in CHROMOSOMES}
    nonclone_rates = {c: [het_rate[s][c] for s in nonclone_samples if not np.isnan(het_rate[s][c])] for c in CHROMOSOMES}

    print(f"\nAbsolute het rate (AD-based, min_alt={MIN_ALT}):")
    print(f"{'chr':>6}  {'clone':>10}  {'non-clone':>10}  {'diff':>8}  {'p':>10}")
    print("-" * 55)
    abs_pvals = []
    for c in CHROMOSOMES:
        cm = np.mean(clone_rates[c])
        nm = np.mean(nonclone_rates[c])
        _, p = stats.mannwhitneyu(clone_rates[c], nonclone_rates[c], alternative='two-sided')
        abs_pvals.append(p)
        sig = '*' if p < 0.05 else ''
        print(f"{c:>6}  {cm:.6f}    {nm:.6f}    {cm-nm:+.6f}  {p:>10.4f} {sig}")

    # ── Within-sample normalization ──
    clone_norm    = {c: [] for c in CHROMOSOMES}
    nonclone_norm = {c: [] for c in CHROMOSOMES}
    for s in kept_samples:
        rates = [het_rate[s][c] for c in CHROMOSOMES]
        valid = [r for r in rates if not np.isnan(r)]
        if not valid:
            continue
        gm = np.mean(valid)
        for c in CHROMOSOMES:
            if not np.isnan(het_rate[s][c]):
                norm = het_rate[s][c] / gm
                if s in CLONE_SAMPLES:
                    clone_norm[c].append(norm)
                else:
                    nonclone_norm[c].append(norm)

    print(f"\nNormalized het rate:")
    print(f"{'chr':>6}  {'clone':>10}  {'non-clone':>10}  {'diff':>8}  {'p':>10}")
    print("-" * 55)
    norm_pvals = []
    for c in CHROMOSOMES:
        cm = np.mean(clone_norm[c])
        nm = np.mean(nonclone_norm[c])
        _, p = stats.mannwhitneyu(clone_norm[c], nonclone_norm[c], alternative='two-sided')
        norm_pvals.append(p)
        sig = '*' if p < 0.05 else ''
        print(f"{c:>6}  {cm:.6f}    {nm:.6f}    {cm-nm:+.6f}  {p:>10.4f} {sig}")

    nc = len(clone_samples)
    nn = len(nonclone_samples)

    # ── Plot 1: bar chart ──
    x     = np.arange(len(CHROMOSOMES))
    width = 0.35
    clone_means    = [np.mean(clone_rates[c])    for c in CHROMOSOMES]
    nonclone_means = [np.mean(nonclone_rates[c]) for c in CHROMOSOMES]
    clone_ses      = [np.std(clone_rates[c])    / np.sqrt(len(clone_rates[c]))    for c in CHROMOSOMES]
    nonclone_ses   = [np.std(nonclone_rates[c]) / np.sqrt(len(nonclone_rates[c])) for c in CHROMOSOMES]
    ratios         = [cm / nm for cm, nm in zip(clone_means, nonclone_means)]

    fig, axes = plt.subplots(1, 2, figsize=(18, 6))

    ax = axes[0]
    ax.bar(x - width/2, clone_means, width, yerr=clone_ses, color="#d62728", alpha=0.7,
           label=f"Clones (n={nc})", capsize=4, edgecolor="white")
    ax.bar(x + width/2, nonclone_means, width, yerr=nonclone_ses, color="#1f77b4", alpha=0.7,
           label=f"Non-clones (n={nn})", capsize=4, edgecolor="white")
    ax.set_xticks(x); ax.set_xticklabels(CHROMOSOMES, fontsize=9)
    ax.set_xlabel("Chromosome", fontsize=12)
    ax.set_ylabel("Het site density (AD-based: min_alt ≥ 3)", fontsize=11)
    ax.set_title(f"Het site density per chromosome — AD-based calling\n"
                 f"(ploidy-agnostic; min total DP={MIN_DEPTH}, min alt reads={MIN_ALT})",
                 fontsize=12, fontweight="bold")
    ax.legend(fontsize=10)

    ax = axes[1]
    colors = ["#d62728" if abs_pvals[i] < 0.05 else "steelblue" for i in range(len(CHROMOSOMES))]
    ax.bar(x, ratios, color=colors, alpha=0.8, edgecolor="white")
    ax.axhline(1.0, color="black", ls="--", lw=1)
    ax.set_xticks(x); ax.set_xticklabels(CHROMOSOMES, fontsize=9)
    ax.set_xlabel("Chromosome", fontsize=12)
    ax.set_ylabel("Clone / Non-clone ratio", fontsize=12)
    ax.set_title("Het density ratio (red = p<0.05)", fontsize=12, fontweight="bold")

    plt.tight_layout()
    out1 = OUT_DIR / "het_density_ad_barplot.png"
    plt.savefig(out1, dpi=200, bbox_inches="tight")
    print(f"\nSaved: {out1}")
    plt.close()

    # ── Plot 2: normalized boxplot ──
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

    ax.axhline(1.0, color="gray", ls="--", lw=1, alpha=0.5)

    for i, c in enumerate(CHROMOSOMES):
        if norm_pvals[i] < 0.05:
            y_max = max(max(clone_norm[c]) if clone_norm[c] else 1,
                        max(nonclone_norm[c]) if nonclone_norm[c] else 1)
            ax.text(i, y_max + 0.01, '*', ha='center', fontsize=14, fontweight='bold')

    ax.set_xticks(x); ax.set_xticklabels(CHROMOSOMES, fontsize=10)
    ax.set_xlabel("Chromosome", fontsize=12)
    ax.set_ylabel("Normalized het density\n(relative to sample genome-wide mean)", fontsize=11)
    ax.set_title(
        f"Normalized het density: Clones (red, n={nc}) vs Non-clones (blue, n={nn})\n"
        f"AD-based calling (ploidy-agnostic): min DP={MIN_DEPTH}, min alt reads={MIN_ALT}",
        fontsize=12, fontweight="bold"
    )
    ax.legend(handles=[Patch(facecolor='#d62728', alpha=0.6, label=f'Clones (n={nc})'),
                        Patch(facecolor='#1f77b4', alpha=0.6, label=f'Non-clones (n={nn})')],
              fontsize=10)

    plt.tight_layout()
    out2 = OUT_DIR / "het_density_ad_normalized_boxplot.png"
    plt.savefig(out2, dpi=200, bbox_inches="tight")
    print(f"Saved: {out2}")
    plt.close()

    print(f"\nAll done. Results in {OUT_DIR}")


if __name__ == "__main__":
    main()
