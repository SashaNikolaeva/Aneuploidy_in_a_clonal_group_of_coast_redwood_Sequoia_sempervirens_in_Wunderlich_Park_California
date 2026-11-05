#!/usr/bin/env python3
# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Within-clone genetic distance versus distance between stems (Spearman rho, 9,999
# permutations; Methods and Results).
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
"""Within-clone genetic distance vs metres between stems (exploratory; cited in the Results, no figure).

Restores the test that was Fig 2D for a short time (redwood-clone-paper e5eb577) so the numbers in the text have a saved
script and output. 36 stem pairs; corrected IBD matrix; V5 GPS from the tree inventory in the Fig 1B local projection;
Spearman rho, significance from 9,999 permutations of stem labels (Mantel-type; pairs sharing a stem are not independent).
Output: clone_paper_figures/exploratory/within_clone_spacing.tsv (pairs) and stdout (rho, p).
Run: ~/miniconda3/envs/quast/bin/python3 -I scripts/exploratory/within_clone_spacing.py
"""
import numpy as np
import pandas as pd

S1 = "/space/s1/sashanikolaeva"
MAT = f"{S1}/Results/Population_structure/April2025files/genetic_distance_matrix_IBD_corrected.csv"
INV = f"{S1}/Scripts/GitHub/redwood-reference-genome/results/popmap/redwood_tree_inventory.tsv"
OUT = f"{S1}/Scripts/GitHub/clone_paper_figures/exploratory/within_clone_spacing.tsv"
CLONE = ["13_WP1", "6_WP5", "31_c3", "31_WP2_1", "28_WP7", "28_WP6", "14_WP3", "22_WPd2", "33_c6"]

D = pd.read_csv(MAT, index_col=0)
D.index = [s.replace("RedwoodCap_", "") for s in D.index]
D.columns = list(D.index)
inv = pd.read_csv(INV, sep="\t", low_memory=False).set_index("sample_id")
lat, lon = inv.loc[CLONE, "lat"], inv.loc[CLONE, "lon"]
c_lat, c_lon = lat.mean(), lon.mean()
x = (lon - c_lon) * 111320 * np.cos(np.radians(c_lat))
y = (lat - c_lat) * 110574
G = D.loc[CLONE, CLONE].values
M = np.hypot(x.values[:, None] - x.values[None, :], y.values[:, None] - y.values[None, :])
iu = np.triu_indices(len(CLONE), 1)


def srho(g, m):
    return pd.Series(g[iu]).corr(pd.Series(m[iu]), method="spearman")


rho = srho(G, M)
rng = np.random.default_rng(7)
null = np.array([srho(G[np.ix_(q, q)], M) for q in (rng.permutation(len(CLONE)) for _ in range(9999))])
p = (1 + np.sum(np.abs(null) >= abs(rho))) / (1 + len(null))
pd.DataFrame({"a": [CLONE[i] for i in iu[0]], "b": [CLONE[j] for j in iu[1]], "metres": M[iu], "d": G[iu]}).to_csv(
    OUT, sep="\t", index=False, float_format="%.5f")
print(f"36 pairs; stems {M[iu].min():.0f}-{M[iu].max():.0f} m apart; Spearman rho {rho:.3f}; permutation p {p:.3f}")
