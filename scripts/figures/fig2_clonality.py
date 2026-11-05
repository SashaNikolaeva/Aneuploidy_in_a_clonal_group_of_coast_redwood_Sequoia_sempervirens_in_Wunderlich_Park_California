#!/usr/bin/env python3
# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Figure 2 (distance heatmap and pairwise distance classes) and the all-pairs comparison with
# the survey (Results: 24,976 pairs of distinct trees).
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
"""New Figure 2 (merges old Figs 2 and 3; 2x2 layout per Sasha 2026-10-26): is the Wunderlich Park group one clone?

Distance: PCR-corrected IBD distance, the matrix the popstructure paper uses (PCoA_from_IBD.py, mantel_test_corrected.py):
Results/Population_structure/April2025files/genetic_distance_matrix_IBD_corrected.csv (dated 2026-05-14; supersedes the
September 2024 one_of_six_matrix/normalized_recalculated_matrix_september_IBD.csv used in the first draft).
Sasha chose this matrix on 2026-10-26; the old figures used the April 2025 normalized Manhattan matrix.

Panels
  A  heatmap: 9 clone ramets, WPR (6.7 m), C1 (~240 m), the two same-tissue technical repeats (MRD, 416_1) and
     the 8 geographically nearest other southern trees (one library per genet). Diagonal masked.
  B  every pairwise distance, by comparison class (points + median):
       technical repeat (same tissue, sequenced twice): MRD, 416_1
       same tree, re-sampled: 286, 256, 17 (same tree label, separate sample)
       (other natural clones BH2/BH3, 279, 292 removed from the panel at Sasha's request, 2026-10-26)
       within the Wunderlich clone: 36 pairs
       clone vs WPR (6.7 m from the nearest ramet), clone vs C1 (173-309 m), clone vs other southern trees (lat < 37.82,
       south of San Francisco Bay), southern background (lat < 37.82, distinct
       genets, clone excluded)
     Dashed line: largest within-clone distance.
  (C distance vs km and D all-pairs histogram were removed 2026-10-26 at Sasha's request; earlier D versions:
   chr4-excluded distances, distance vs stem spacing, southern PCoA. The all-pairs background statistics are still
   written to fig2_allpairs_background.tsv for the text.)
Same-genet pairs 9.6-153 km apart (G001, G003, G004, G012 in het_cause T5) are left out of B as unresolved
sample-identity questions.
Output: clone_paper_figures/clonality/fig2_clonality.{png,pdf}, fig2_clonality_pairs.tsv, fig2_clonality_summary.tsv,
        fig2_allpairs_background.tsv
Run: ~/miniconda3/envs/quast/bin/python3 -I scripts/figures/fig2_clonality.py
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402
from scipy.cluster.hierarchy import leaves_list, linkage  # noqa: E402
from scipy.spatial.distance import squareform  # noqa: E402

S1 = "/space/s1/sashanikolaeva"
MAT = f"{S1}/Results/Population_structure/April2025files/genetic_distance_matrix_IBD_corrected.csv"
INV = f"{S1}/Scripts/GitHub/redwood-reference-genome/results/popmap/redwood_tree_inventory.tsv"
OUT = f"{S1}/Scripts/GitHub/clone_paper_figures/clonality"
os.makedirs(OUT, exist_ok=True)

CLONE = ["13_WP1", "6_WP5", "31_c3", "31_WP2_1", "28_WP7", "28_WP6", "14_WP3", "22_WPd2", "33_c6"]
TECH = [("22_MRD", "28_MRD"), ("16_416-1", "35_416_1")]
# tree 286 dropped as a replicate pair (popgen keeps only 3_286-1; Sasha 2026-10-26, DV14)
SAME_TREE = [("1_256S", "33_256"), ("21_17-1", "23_17")]
SOUTH = 37.82  # south of San Francisco Bay (population break); 38.0 took in Muir Woods and 19_FL



def readmat(path):
    m = pd.read_csv(path, index_col=0)
    m.index = [s.replace("RedwoodCap_", "") for s in m.index]
    m.columns = [s.replace("RedwoodCap_", "") for s in m.columns]
    return m


D = readmat(MAT)
inv = pd.read_csv(INV, sep="\t", low_memory=False).set_index("sample_id")
for s in CLONE + ["22_WPR", "16_C1"] + [x for p in TECH + SAME_TREE for x in p]:
    assert s in D.index, s

lat, lon = inv["lat"], inv["lon"]
c_lat, c_lon = lat[CLONE].mean(), lon[CLONE].mean()


def km(a_lat, a_lon, b_lat, b_lon):
    r = np.pi / 180
    h = np.sin((b_lat - a_lat) * r / 2) ** 2 + np.cos(a_lat * r) * np.cos(b_lat * r) * np.sin((b_lon - a_lon) * r / 2) ** 2
    return 2 * 6371.0 * np.arcsin(np.sqrt(h))


# population: canonical 224, one library per genet, clone neighbours handled separately
pop = inv[(inv.in_canonical_224 == True) & inv.index.isin(D.index)]  # noqa: E712
# population = exactly the popgen paper's 224 trees, incl. 7_241-4 (Sasha 2026-10-26, DV14)
pop = pop[~pop.index.isin(CLONE)]
assert len(pop) == 224 and not pop.genet_id.duplicated().any(), len(pop)
pop = pop.assign(dist_km=km(c_lat, c_lon, pop.lat.values, pop.lon.values),
                 d_clone=D.loc[pop.index, CLONE].mean(1).values)

# ---- B: pairwise classes
rows = []


def add(cls, a, b):
    rows.append(dict(cls=cls, a=a, b=b, d=float(D.loc[a, b])))


for a, b in TECH:
    add("technical repeat", a, b)
for a, b in SAME_TREE:
    add("same tree,\ncollected twice", a, b)
for i, a in enumerate(CLONE):
    for b in CLONE[i + 1:]:
        add("within\nWunderlich clone", a, b)
for a in CLONE:
    add("clone vs WPR\n(7 m)", a, "22_WPR")
    add("clone vs C1\n(173 m)", a, "16_C1")
south = pop[(pop.lat < SOUTH) & ~pop.index.isin(["22_WPR", "16_C1"])]
for a in CLONE:
    for b in south.index:
        add("clone vs other\nsouthern trees", a, b)
sl = list(south.index)
for i, a in enumerate(sl):
    for b in sl[i + 1:]:
        add("southern\nbackground", a, b)
pairs = pd.DataFrame(rows)
pairs.to_csv(f"{OUT}/fig2_clonality_pairs.tsv", sep="\t", index=False, float_format="%.5f")
order = list(dict.fromkeys(pairs.cls))
summ = pairs.groupby("cls", sort=False).d.agg(["count", "min", "median", "max"]).reindex(order)
summ.index = [i.replace("\n", " ") for i in summ.index]
summ.to_csv(f"{OUT}/fig2_clonality_summary.tsv", sep="\t", float_format="%.4f")
print(summ.to_string())

# ---- A: heatmap set
near = south.sort_values("dist_km").index[:8].tolist()
cl = D.loc[CLONE, CLONE].values.copy()
np.fill_diagonal(cl, 0)
clone_ord = [CLONE[i] for i in leaves_list(linkage(squareform(cl, checks=False), "average"))]
groups = ([("Wunderlich clone", s) for s in clone_ord] + [("neighbour", "22_WPR"), ("neighbour", "16_C1")]
          + [("technical repeat", s) for p in TECH for s in p] + [("other southern", s) for s in near])
hs = [s for _, s in groups]
H = D.loc[hs, hs].values.astype(float)
np.fill_diagonal(H, np.nan)

GCOL = {"Wunderlich clone": "#1b7837", "neighbour": "#e08214", "technical repeat": "#8073ac", "other southern": "#878787"}
CCOL = {"technical repeat": "#8073ac", "same tree,\ncollected twice": "#b2abd2",
        "within\nWunderlich clone": "#1b7837", "clone vs WPR\n(7 m)": "#e08214", "clone vs C1\n(173 m)": "#fdb863",
        "clone vs other\nsouthern trees": "#878787", "southern\nbackground": "#bababa"}

fig = plt.figure(figsize=(16, 8))
gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1], wspace=0.3)
axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[0, 1])

cmap = ListedColormap(plt.get_cmap("viridis_r")(np.linspace(0, 1, 256)))
cmap.set_bad("#ffffff")
vmin, vmax = np.nanmin(H), np.nanmax(H)
im = axA.imshow(H, cmap=cmap, vmin=vmin, vmax=vmax)
n = len(hs)
axA.set_xticks(range(n))
axA.set_yticks(range(n))
lab = [s if g == "technical repeat" else (s.split("_", 1)[1] if s.split("_", 1)[0].isdigit() else s) for g, s in groups]
axA.set_xticklabels(lab, rotation=90, fontsize=9)
axA.set_yticklabels(lab, fontsize=9)
for t, (g, _) in zip(axA.get_xticklabels(), groups):
    t.set_color(GCOL[g])
for t, (g, _) in zip(axA.get_yticklabels(), groups):
    t.set_color(GCOL[g])
for k in np.cumsum([9, 2, 4])[:]:
    axA.axhline(k - 0.5, color="white", lw=2.5)
    axA.axvline(k - 0.5, color="white", lw=2.5)
for i in range(n):
    for j in range(n):
        if i != j and (i < 11 and j < 11):
            axA.text(j, i, f"{H[i, j]:.3f}"[1:], ha="center", va="center", fontsize=6.2,
                     color="white" if H[i, j] > (vmin + vmax) / 2 else "black")
cb = fig.colorbar(im, ax=axA, fraction=0.035, pad=0.02)
cb.set_label("genetic distance (PCR-corrected IBD)")
for g, c in GCOL.items():
    axA.plot([], [], "s", color=c, label=g, ms=9)
axA.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=4, frameon=False, fontsize=9)
axA.set_title("A", loc="left", fontweight="bold", fontsize=15)

rng = np.random.default_rng(1)
for x, c in enumerate(order):
    v = pairs.loc[pairs.cls == c, "d"].values
    jit = rng.uniform(-0.28, 0.28, len(v)) if len(v) > 1 else np.zeros(1)
    a = 0.25 if len(v) > 200 else (0.6 if len(v) > 40 else 0.9)
    axB.scatter(x + jit, v, s=10 if len(v) > 40 else 26, color=CCOL[c], alpha=a, edgecolor="none", zorder=2)
    axB.hlines(np.median(v), x - 0.36, x + 0.36, color="black", lw=1.6, zorder=3)
wmax = pairs.loc[pairs.cls == "within\nWunderlich clone", "d"].max()
axB.axhline(wmax, ls="--", color="#1b7837", lw=1, zorder=1)
axB.set_xticks(range(len(order)))
axB.set_xticklabels([o.replace("\n", " ") for o in order], fontsize=8.5, rotation=35, ha="right", rotation_mode="anchor")
axB.set_ylabel("genetic distance (PCR-corrected IBD)")
axB.spines[["top", "right"]].set_visible(False)
lo = pairs.d.min()
axB.set_ylim(lo - 0.012, pairs.d.max() + 0.006)
for x, c in enumerate(order):
    axB.text(x, lo - 0.009, f"n={int((pairs.cls == c).sum())}", ha="center", fontsize=7.5, color="0.35")
axB.set_title("B", loc="left", fontweight="bold", fontsize=15)

# panels C (distance vs km) and D (all-pairs histogram) removed 2026-10-26 at Sasha's request; the all-pairs background
# statistics are still written because the text cites them (24,753 pairs, min 0.136, median 0.152)
pn = list(pop.index)
bg = D.loc[pn, pn].values[np.triu_indices(len(pn), 1)]
pd.DataFrame({"n_pairs": [len(bg)], "min": [bg.min()], "q001": [np.quantile(bg, 0.001)], "median": [np.median(bg)]}).to_csv(
    f"{OUT}/fig2_allpairs_background.tsv", sep="\t", index=False, float_format="%.5f")

fig.savefig(f"{OUT}/fig2_clonality.png", dpi=200, bbox_inches="tight")
fig.savefig(f"{OUT}/fig2_clonality.pdf", bbox_inches="tight")
print("heatmap neighbours:", near)
print(f"background: {len(bg)} pairs of distinct survey trees, min {bg.min():.4f}, median {np.median(bg):.4f}")
print("WPR km", round(pop.dist_km["22_WPR"], 3), "C1 km", round(pop.dist_km["16_C1"], 3), "n pop", len(pop))
