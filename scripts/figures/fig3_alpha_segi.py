# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Figure 3 (alpha per chromosome for the stems, WPR, C1 and the survey trees), plotted from
# the saved output of scripts/finalize/09_R0_fig6_table1.py.
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
"""Fig 3 (aneuploidy alpha plot) re-plotted from the saved outputs of finalize/09_R0_fig6_table1.py, without the
14_364 highlight (Sasha, 2026-11-05: 14_364 is no longer discussed in the text). Nothing is recomputed: alpha and
the expected levels/bands are read from analysis/finalize_run_224/{fig6_source,fig6_expected_levels}.tsv; 14_364 is
drawn as an ordinary survey tree. Style follows 09_R0_fig6_table1.py.
Output: clone_paper_figures/aneuploidy/fig3_alpha_segi.{png,pdf}
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

REPO = "/space/s1/sashanikolaeva/Scripts/GitHub/redwood-clone-paper"
RUN = f"{REPO}/analysis/finalize_run_224"
FD = "/space/s1/sashanikolaeva/Scripts/GitHub/clone_paper_figures/aneuploidy"
WPR, C1 = "22_WPR", "16_C1"
CH = [f"chr{i}" for i in range(1, 12)]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#d6d5d0"
C_CLONE, C_WPR, C_C1 = "#2a78d6", "#eb6834", "#1baf7a"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "font.family": "DejaVu Sans"})

src = pd.read_csv(f"{RUN}/fig6_source.tsv", sep="\t", index_col=0)
lv = pd.read_csv(f"{RUN}/fig6_expected_levels.tsv", sep="\t", index_col=0).loc[CH]
alpha = src[[f"alpha_{c}" for c in CH]]
alpha.columns = CH
CLONE = list(src.index[src.group == "clone_ramet"])
assert len(CLONE) == 9 and WPR in alpha.index and C1 in alpha.index
zB = float(((lv.band_hi - 1) / lv.sd_alpha_PANEL_ALL).median())
a5, a7 = lv.alpha5_one_copy_loss, lv.alpha7_one_copy_gain

x = np.arange(1, 12)
fig, ax = plt.subplots(figsize=(9.5, 6.0))
for i, c in enumerate(CH):
    ax.add_patch(plt.Rectangle((x[i] - 0.32, lv.band_lo[c]), 0.64, lv.band_hi[c] - lv.band_lo[c],
                               color="#ecebe7", lw=0, zorder=0))
    for lvl in (a5[c], a7[c]):
        ax.hlines(lvl, x[i] - 0.34, x[i] + 0.34, color=INK2, lw=1.1, ls=(0, (3, 2)), zorder=1)
hl = set(CLONE) | {WPR, C1}
others = [s for s in alpha.index if s not in hl]
for s in others:
    ax.plot(x, alpha.loc[s].values, color="#9a9994", lw=0.5, alpha=0.35, zorder=2)
for s in CLONE:
    ax.plot(x, alpha.loc[s].values, color=C_CLONE, lw=1.4, alpha=0.9, zorder=4)
for s, col in [(WPR, C_WPR), (C1, C_C1)]:
    ax.plot(x, alpha.loc[s].values, color=col, lw=1.6, zorder=5)
ax.annotate("9 clone ramets", (4, alpha.loc[CLONE, "chr4"].min()), xytext=(4.45, 0.835), fontsize=8.5, color=INK,
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
     plt.Line2D([], [], color="#9a9994", lw=0.8, label=f"other survey trees (n = {len(others)})"),
     plt.Line2D([], [], color=INK2, lw=1.1, ls=(0, (3, 2)), label="expected for one copy lost / gained"),
     plt.Rectangle((0, 0), 1, 1, color="#ecebe7", label=f"population range (±{zB:.1f} SD)")]
ax.legend(handles=h, frameon=False, fontsize=7.8, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=3)
ax.set_ylim(0.80, 1.20)
fig.tight_layout()
fig.savefig(f"{FD}/fig3_alpha_segi.png", dpi=250)
fig.savefig(f"{FD}/fig3_alpha_segi.pdf")
print(f"n other = {len(others)}; zB = {zB:.2f}; written to {FD}/fig3_alpha_segi.{{png,pdf}}")
