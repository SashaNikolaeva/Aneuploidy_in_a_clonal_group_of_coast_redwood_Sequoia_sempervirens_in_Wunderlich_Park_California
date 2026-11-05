#!/usr/bin/env python3
# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Combines the allele-dosage and clonality results into the calls set in advance and draws a
# summary figure (not in the paper).
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
"""06: combine A9, lattice, B1/B4 and B2 into the pre-registered calls (PREREG_v2 A9d, L2/L3, Bx) and draw the summary
figure (clone_paper_figures/dosage_clonality/dosage_clonality_summary.{png,pdf}; colours as in Fig 2).
Outputs: RUN/06_calls.json; stdout -> RUN/06_summary.log
"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dc_common as c  # noqa: E402

A = json.load(open(f"{c.RUN}/02_A9_summary.json"))
L = json.load(open(f"{c.RUN}/03_lattice_gate.json")) if os.path.exists(f"{c.RUN}/03_lattice_gate.json") else None
B1 = json.load(open(f"{c.RUN}/04_B1_summary.json"))
B2 = json.load(open(f"{c.RUN}/05_B2_summary.json"))
b1, b4 = B1["corrected"], B1["uncorrected"]
contra = [s for s in b1["contradicted_candidates_B1"] if s in B2["B2_flagged_stems"]]
if contra:
    bcall = "CONTRADICTED"
elif b1["call_B1_part"] == "ONE GENET" and B2["call_B2"] == "ONE GENET":
    bcall = "ONE GENET"
else:
    bcall = "AMBIGUOUS"
robust_b = b4["call_B1_part"] == b1["call_B1_part"]
calls = dict(A9=dict(call=A["call"], reason_code=A["reason_code"], T=A["T_chr4"], p=A["p"], max_T_c=A["max_T_c"],
                     power=A["power_pi033_delta1_15"], p_depth200=A["p200"]),
             lattice=None if L is None else dict(gate=L["gate"], max_boot_sd=L["max_boot_sd"],
                                                 flagged=L["n_flagged_lo_gt_0p03"]),
             B=dict(call=bcall, B1=b1["call_B1_part"], B2=B2["call_B2"], B4_same_as_B1=robust_b,
                    contradicted_stems=contra, B1_U=b1["U"], B2_U=B2["U"], B2_n_clone_out=B2["n_clone_out"]))
json.dump(calls, open(f"{c.RUN}/06_calls.json", "w"), indent=1, default=float)
print(json.dumps(calls, indent=1, default=float))

# ---------------------------------------------------------------- figure
GREEN, ORANGE, LORANGE, PURPLE, GREY = "#1b7837", "#e08214", "#fdb863", "#8073ac", "#9a9a9a"
null = np.load(f"{c.RUN}/02_A9_null.npy")
pc = pd.read_csv(f"{c.RUN}/02_A9_per_chromosome.tsv", sep="\t")
P1 = pd.read_csv(f"{c.RUN}/04_B1_pairs_corrected.tsv", sep="\t")
P2 = pd.read_csv(f"{c.RUN}/05_B2_pairs.tsv", sep="\t")
fig, ax = plt.subplots(1, 3, figsize=(16, 4.8), gridspec_kw=dict(wspace=0.32))

a = ax[0]
cnt, _, _ = a.hist(null, bins=60, color="#d9d9d9", edgecolor="white", linewidth=0.5, label="10,000 random sets of 9 survey trees")
a.set_ylim(0, cnt.max() * 1.55)  # headroom for the legend
a.axvline(A["T_chr4"], color=GREEN, lw=2.2, label=f"clone, chromosome 4 (p = {A['p']:.3f})")
for t in pc.T_c_control.dropna():
    a.axvline(t, color=GREEN, lw=0.8, alpha=0.45, ymax=0.25)
a.axvline(A["null_q999"], color="0.3", ls="--", lw=1, label="99.9th percentile (decision threshold)")
a.plot([], [], color=GREEN, lw=0.8, alpha=0.45, label="clone, each other chromosome")
a.set_xlabel("standardized chr4 excess deviation from the survey (T)")
a.set_ylabel("number of random sets")
a.legend(frameon=False, fontsize=7.5, loc="upper left")
a.spines[["top", "right"]].set_visible(False)
a.set_title("A", loc="left", fontweight="bold", fontsize=14)


def strip(axx, P, col, U, ylab, logy=False):
    order = [("identity", "same tree or tissue (7)", PURPLE), ("clone", "within the clone (36)", GREEN),
             ("anchor_WPR", "clone vs WPR (9)", ORANGE), ("anchor_C1", "clone vs C1 (9)", LORANGE),
             ("anchor_south", "clone vs southern trees (396)", GREY)]
    rng = np.random.default_rng(1)
    for i, (k, lab, cc) in enumerate(order):
        v = P.loc[P.cls == k, col].values
        axx.scatter(i + rng.uniform(-0.25, 0.25, len(v)), v, s=14 if len(v) > 50 else 26, color=cc,
                    alpha=0.35 if len(v) > 50 else 0.9, edgecolor="none")
    axx.axhline(U, color=PURPLE, ls="--", lw=1)
    axx.text(4.45, U, "99% bound of\nsame-tree pairs", color=PURPLE, fontsize=7.5, va="bottom", ha="right")
    axx.set_xticks(range(5))
    axx.set_xticklabels([o[1] for o in order], rotation=30, ha="right", fontsize=8)
    axx.set_ylabel(ylab)
    if logy:
        axx.set_yscale("log")
    axx.spines[["top", "right"]].set_visible(False)


strip(ax[1], P1, "d_c", b1["U"], "batch-corrected genetic distance (IBD)")
ax[1].set_title("B", loc="left", fontweight="bold", fontsize=14)
strip(ax[2], P2, "rate", B2["U"], "allele-presence discordance (thinned to 30 reads; log)", logy=True)
ax[2].set_title("C", loc="left", fontweight="bold", fontsize=14)
fig.savefig(f"{c.FIG}/dosage_clonality_summary.png", dpi=200, bbox_inches="tight")
fig.savefig(f"{c.FIG}/dosage_clonality_summary.pdf", bbox_inches="tight")
print("figure written")
