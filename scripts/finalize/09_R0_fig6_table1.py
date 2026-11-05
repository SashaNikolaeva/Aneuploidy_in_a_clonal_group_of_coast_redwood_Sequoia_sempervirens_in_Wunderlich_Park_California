# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Computes alpha_ij (read share of a sample on a chromosome divided by the median share of
# the 224 survey trees) and the expected one-copy levels, saved as fig6_source.tsv and fig6_expected_levels.tsv
# (Methods and Results: Chromosome copy number; Supplementary Table S1). Also writes the stem coordinates of
# Table 1. Its figure was replaced by scripts/figures/fig3_alpha_segi.py.
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
"""R0b Fig 6 regenerated from the SEGI MQ60 read counts, R0c Table 1 additions, and a supplementary chr4 positional
profile figure from R1c (decision 2.1, 4.2).

Fig 6: alpha(s, c) = share(s, c) / median_PANEL_ALL share(c), share over chr1-11 (as the published figure; d04).
One grey line per sample (all 282 samples with > 1 M reads); the 9 ramets coloured; WPR, C1, 14_364 labelled.
Expected one-copy levels alpha5 = (5/6)/(1 - p_c/6), alpha7 = (7/6)/(1 + p_c/6), p_c = PANEL_ALL median share of c;
per-chromosome bands 1 +- z_B * SD_PANEL_ALL(alpha), z_B from R2 (Bonferroni, printed in 04_R2_screen.log).
Outputs: clone_paper_figures/aneuploidy/fig6_alpha_segi.{png,pdf}, figS_chr4_positional_profile.{png,pdf};
analysis/finalize_run/fig6_source.tsv, fig6_expected_levels.tsv, R0_table1_additions.tsv
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *  # noqa

FD = f"{FIGDIR}/aneuploidy"
os.makedirs(FD, exist_ok=True)
INK, INK2, GRID = "#0b0b0b", "#52514e", "#d6d5d0"
C_CLONE, C_WPR, C_C1, C_364 = "#2a78d6", "#eb6834", "#1baf7a", "#e87ba4"   # reference palette slots 1-3, 5
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "font.family": "DejaVu Sans"})

S = sample_sets()
rc, P = S["rc"], S["PANEL_ALL"]
samples = S["plot"]
share = pd.DataFrame({s: {c: rc[s][c] / sum(rc[s].values()) for c in CH} for s in samples}).T[CH]
med = share.loc[P].median()
alpha = share / med
sd = alpha.loc[P].std(ddof=1)
nscreen = sum(1 for s in samples if s not in CLONE)
zB = stats.norm.isf(0.05 / (2 * nscreen * 11))
a5 = (5 / 6) / (1 - med / 6)
a7 = (7 / 6) / (1 + med / 6)
cells = pd.read_csv(f"{OUT}/R2_screen_cells.tsv", sep="\t")
chat = cells.pivot(index="sample", columns="chrom", values="c_hat")[CH]
src = alpha.add_prefix("alpha_").join(chat.add_prefix("c_hat_"))
src.insert(0, "group", ["clone_ramet" if s in CLONE else ("PANEL_ALL" if s in P else "other_sample") for s in src.index])
src.to_csv(f"{OUT}/fig6_source.tsv", sep="\t", float_format="%.4f", index_label="sample")
lv = pd.DataFrame({"panel_median_share_p_c": med, "alpha5_one_copy_loss": a5, "alpha7_one_copy_gain": a7,
                   "sd_alpha_PANEL_ALL": sd, "band_lo": 1 - zB * sd, "band_hi": 1 + zB * sd})
lv.to_csv(f"{OUT}/fig6_expected_levels.tsv", sep="\t", float_format="%.4f", index_label="chrom")
print(f"z_B {zB:.3f}; chr4 p_c {med['chr4']:.4f}, alpha5 {a5['chr4']:.4f}; clone chr4 alpha "
      f"{alpha.loc[CLONE, 'chr4'].min():.3f}-{alpha.loc[CLONE, 'chr4'].max():.3f}")

# ---------------- Fig 6
x = np.arange(1, 12)
fig, ax = plt.subplots(figsize=(9.5, 6.0))
for i, c in enumerate(CH):
    ax.add_patch(plt.Rectangle((x[i] - 0.32, 1 - zB * sd[c]), 0.64, 2 * zB * sd[c], color="#ecebe7", lw=0, zorder=0))
    for lvl in (a5[c], a7[c]):
        ax.hlines(lvl, x[i] - 0.34, x[i] + 0.34, color=INK2, lw=1.1, ls=(0, (3, 2)), zorder=1)
hl = set(CLONE) | {WPR, C1, S364}
for s in samples:
    if s not in hl:
        ax.plot(x, alpha.loc[s].values, color="#9a9994", lw=0.5, alpha=0.35, zorder=2)
for s in CLONE:
    ax.plot(x, alpha.loc[s].values, color=C_CLONE, lw=1.4, alpha=0.9, zorder=4)
for s, col, lab in [(WPR, C_WPR, "WPR"), (C1, C_C1, "C1"), (S364, C_364, "14_364")]:
    ax.plot(x, alpha.loc[s].values, color=col, lw=1.6, zorder=5)
ax.annotate("9 clone ramets", (4, alpha.loc[CLONE, "chr4"].min()), xytext=(4.45, 0.835), fontsize=8.5, color=INK,
            arrowprops=dict(arrowstyle="-", color=INK2, lw=0.6))
ax.annotate("14_364", (4, alpha.loc[S364, "chr4"]), xytext=(2.9, 0.915), fontsize=8.5, color=INK,
            arrowprops=dict(arrowstyle="-", color=INK2, lw=0.6))
ax.annotate("14_364", (9, alpha.loc[S364, "chr9"]), xytext=(9.25, 1.10), fontsize=8.5, color=INK,
            arrowprops=dict(arrowstyle="-", color=INK2, lw=0.6))
ax.text(4.36, a5["chr4"] + 0.003, "5 copies", fontsize=7.5, color=INK2, va="bottom")
ax.text(4.36, a7["chr4"] + 0.004, "7 copies", fontsize=7.5, color=INK2, va="bottom")
ax.set_xticks(x)
ax.set_xticklabels([str(i) for i in x])
ax.set_xlim(0.5, 11.5)
ax.set_xlabel("giant sequoia (SEGI v2.0) chromosome")
ax.set_ylabel(r"$\alpha$ = read share / population median share")
ax.yaxis.grid(True, color=GRID, lw=0.6)
ax.set_axisbelow(True)
ax.spines[["top", "right"]].set_visible(False)
h = [plt.Line2D([], [], color=C_CLONE, lw=1.6, label="Wunderlich clone (9 ramets)"),
     plt.Line2D([], [], color=C_WPR, lw=1.6, label="WPR"), plt.Line2D([], [], color=C_C1, lw=1.6, label="C1"),
     plt.Line2D([], [], color=C_364, lw=1.6, label="14_364"),
     plt.Line2D([], [], color="#9a9994", lw=0.8, label=f"other survey {'trees' if SET == 'popgen224' else 'samples'} (n = {len(samples) - 12})"),
     plt.Line2D([], [], color=INK2, lw=1.1, ls=(0, (3, 2)), label="expected for one copy lost / gained"),
     plt.Rectangle((0, 0), 1, 1, color="#ecebe7", label=f"population range (±{zB:.1f} SD)")]
ax.legend(handles=h, frameon=False, fontsize=7.8, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=4)
ax.set_ylim(0.80, 1.20)
fig.tight_layout()
fig.savefig(f"{FD}/fig6_alpha_segi.png", dpi=250)
fig.savefig(f"{FD}/fig6_alpha_segi.pdf")
plt.close(fig)

# ---------------- Table 1 additions
inv = S["inv"]
cen = S["centroid"]
rows = []
for s in CLONE + [WPR, C1]:
    r = inv[s]
    near = min(hav((float(r["lat"]), float(r["lon"])), (float(inv[t]["lat"]), float(inv[t]["lon"]))) * 1000
               for t in CLONE if t != s)
    rows.append(dict(sample=s, genet_id=r["genet_id"], tissue="leaf", pcr_pool=int(float(r["pcr_pool"])),
                     pcr_group=int(float(r["pcr_group"])), plate_well=r["plate_well"], pcr_date=r["pcr_date"],
                     segi_reads_chr1_11=int(sum(rc[s].values())), lat=float(r["lat"]), lon=float(r["lon"]),
                     m_to_nearest_ramet=round(near, 1), m_to_clone_centroid=round(hav(cen, (float(r["lat"]), float(r["lon"]))) * 1000, 1),
                     qc_flags=r["qc_flags"]))
t1 = pd.DataFrame(rows)
t1.to_csv(f"{OUT}/R0_table1_additions.tsv", sep="\t", index=False)
print(t1.to_string(index=False))

# ---------------- supplementary: chr4 positional profile (R1c) + 14_364 chr4 / chr9 (R2-364)
b = pd.read_csv(f"{OUT}/R1c_bins.tsv", sep="\t")
fig, axs = plt.subplots(1, 3, figsize=(12, 3.6), sharey=True, gridspec_kw=dict(width_ratios=[37, 37, 32]))
for ax, c, s, col, ttl in [(axs[0], "chr4", "clone_mean_log2", C_CLONE, "A  Wunderlich clone, chr4 (mean of 9 ramets)"),
                            (axs[1], "chr4", S364, C_364, "B  14_364, chr4"), (axs[2], "chr9", S364, C_364, "C  14_364, chr9")]:
    d = b[b.chrom == c]
    xm = d.bin_start_Mb + 10
    ax.fill_between(xm, d.panel_q05, d.panel_q95, color="#ecebe7", lw=0, label="population 5-95%")
    ax.axhline(0, color=INK2, lw=0.7)
    ax.axhline(np.log2(5 / 6), color=INK2, lw=1, ls=(0, (3, 2)))
    ax.axhline(np.log2(7 / 6), color=INK2, lw=1, ls=(0, (3, 2)))
    if s == "clone_mean_log2":
        for r_ in CLONE:
            ax.plot(xm, d[r_], color=C_CLONE, lw=0.5, alpha=0.3)
    ax.plot(xm, d[s], "o-", color=col, lw=1.6, ms=3.5)
    ax.set_title(ttl, loc="left", fontsize=9.5, color=INK)
    ax.set_xlabel(f"SEGI {c} position (Mb, 20-Mb bins)")
    ax.spines[["top", "right"]].set_visible(False)
    ax.yaxis.grid(True, color=GRID, lw=0.5)
axs[0].set_ylabel("log2 depth vs population median")
axs[0].text(5, np.log2(5 / 6) - 0.06, "5 copies", fontsize=7.5, color=INK2)
axs[0].text(5, np.log2(7 / 6) + 0.02, "7 copies", fontsize=7.5, color=INK2)
axs[0].set_ylim(-1.0, 1.0)
fig.tight_layout()
fig.savefig(f"{FD}/figS_chr4_positional_profile.png", dpi=250)
fig.savefig(f"{FD}/figS_chr4_positional_profile.pdf")
print("figures written to", FD)
