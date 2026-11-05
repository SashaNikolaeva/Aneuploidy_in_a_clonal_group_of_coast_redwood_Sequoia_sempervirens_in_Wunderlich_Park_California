# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Clonality tables: classes of pairwise distances, within-clone distance versus library
# depth, and the comparison of distances with and without chromosome 4 (Results: the nine stems are one genet).
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
"""R5 clonality tables, consistent with (not a duplicate of) the parent's scripts/figures/fig2_clonality.py
(decision 2.8). Classes come ONLY from the parent's fig2_clonality_pairs.tsv; this script never re-defines them.

Gate R5-g: class counts and min / median / max recomputed from genetic_distance_matrix_IBD_corrected.csv for the
parent's pairs equal the parent's fig2_clonality_summary.tsv to 4 decimals (and the decision's quoted values).
On mismatch: write R5_GATE_FAIL.md, print, exit 3 (no substitution).
T1 replicate ceiling + gap; T2 the 36 within-clone pairs with pools / groups / depths, Spearman distance ~ min depth;
T3 chr4-excluded sensitivity on the uncorrected full vs no-chr4 matrices; T4 southern-set definition check.
Outputs: R5_gate.tsv, R5_T1_classes.tsv, R5_T2_within_clone.tsv, R5_T3_no_chr4.tsv, R5_T4_south.tsv; log 08_R5_clonality.log
"""
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *  # noqa
import vcfcache  # noqa


def readmat(p):
    D = pd.read_csv(p, index_col=0)
    D.index = [short(s) for s in D.index]
    D.columns = [short(s) for s in D.columns]
    return D


pairs = pd.read_csv(INPUTS["fig2_pairs"], sep="\t")
pairs["cls1"] = pairs.cls.str.replace("\n", " ")
summ = pd.read_csv(INPUTS["fig2_summary"], sep="\t", index_col=0)
D = readmat(INPUTS["ibd_corrected"])
pairs["d_mat"] = [float(D.loc[a, b]) for a, b in zip(pairs.a, pairs.b)]
maxdiff = (pairs.d_mat - pairs.d).abs().max()
re_ = pairs.groupby("cls1", sort=False).d_mat.agg(["count", "min", "median", "max"])
print(f"pairs {len(pairs)}; max |pairs.d - matrix| = {maxdiff:.2e}")
QUOTED = {"within Wunderlich clone": (36, 0.0966, None, 0.1075), "technical repeat": (2, 0.1021, None, 0.1035),
          "clone vs WPR (7 m)": (9, 0.1278, None, 0.1341), "clone vs C1 (~240 m)": (9, 0.1396, None, 0.1440),
          "clone vs C1 (173 m)": (9, 0.1396, None, 0.1440),  # Fig 2 label renamed after the 2026-10-26 run (2026-10-29)
          "southern background": (None, 0.1401, None, 0.1577)}
ok = maxdiff <= 5.5e-6 and list(re_.index) == list(summ.index)  # pairs.tsv is rounded to 5 dp
rows = []
for c in summ.index:
    r, s = re_.loc[c], summ.loc[c]
    m = int(r["count"]) == int(s["count"]) and all(round(r[k], 4) == round(s[k], 4) for k in ["min", "median", "max"])
    q = QUOTED.get(c)
    mq = True if q is None else all(v is None or (round(float(r[k]), 4) == v if k != "count" else int(r[k]) == v)
                                    for k, v in zip(["count", "min", "median", "max"], q))
    ok &= m and mq
    rows.append((c, int(r["count"]), r["min"], r["median"], r["max"], int(s["count"]), s["min"], s["median"], s["max"], m, mq))
with open(f"{OUT}/R5_gate.tsv", "w") as f:
    f.write("class\tn_recomputed\tmin\tmedian\tmax\tn_parent\tmin_parent\tmedian_parent\tmax_parent\tmatch_4dp\tmatch_decision_quote\n")
    for r in rows:
        f.write("\t".join(f"{v:.4f}" if isinstance(v, float) else str(v) for v in r) + "\n")
    f.write(f"GATE\t{'PASS' if ok else 'FAIL'}\n")
print(open(f"{OUT}/R5_gate.tsv").read())
if not ok:
    open(f"{OUT}/R5_GATE_FAIL.md", "w").write("# R5-g failed\n\nRecomputed class summaries differ from the parent's "
                                              "fig2_clonality_summary.tsv; R5 not run (S4). See R5_gate.tsv.\n")
    sys.exit(3)

src = open(INPUTS["fig2_script"]).read()
south_parent = float(re.search(r"^SOUTH\s*=\s*([0-9.]+)", src, re.M).group(1))
print(f"parent fig2_clonality.py SOUTH = {south_parent}; excludes 7_241-4: {'7_241-4' in src}")

# ---- T1
cls = lambda c: pairs.loc[pairs.cls1 == c, "d_mat"]
w = cls("within Wunderlich clone")
wpr = cls("clone vs WPR (7 m)")
same_tree = "same tree, re-sampled" if (pairs.cls1 == "same tree, re-sampled").any() else "same tree, collected twice"  # label renamed in Fig 2 after the run (2026-10-29)
rep_max = max(cls("technical repeat").max(), cls(same_tree).max())
gap = wpr.min() - w.max()
one_genet = (w.max() <= rep_max + 0.005) and (w.max() < wpr.min())
with open(f"{OUT}/R5_T1_classes.tsv", "w") as f:
    f.write("class\tn\tmin\tmedian\tmax\tmean\tsd\n")
    for c in summ.index:
        v = cls(c)
        f.write(f"{c}\t{len(v)}\t{v.min():.4f}\t{v.median():.4f}\t{v.max():.4f}\t{v.mean():.4f}\t"
                f"{v.std(ddof=1) if len(v) > 1 else float('nan'):.4f}\n")
    f.write(f"gap_minWPR_minus_maxWithin\t\t\t\t{gap:.4f}\t\tin within-clone SD units: {gap / w.std(ddof=1):.2f}\n")
    f.write(f"one_genet_rule\t\t\t\t\t\tmax within {w.max():.4f} <= max(technical, same-tree) {rep_max:.4f} + 0.005 "
            f"and < min clone-WPR {wpr.min():.4f}: {one_genet}\n")
print(f"T1: within-clone {w.min():.4f}-{w.max():.4f} (mean {w.mean():.4f}, SD {w.std(ddof=1):.4f}); replicate max {rep_max:.4f}; "
      f"min clone-WPR {wpr.min():.4f}; gap {gap:.4f} = {gap / w.std(ddof=1):.1f} within-clone SD; one genet: {one_genet}")

# ---- T2
S = sample_sets()
inv, rc = S["inv"], S["rc"]
V = vcfcache.load()
vix = {s: i for i, s in enumerate(V["samples"])}
mdp = V["DP"].mean(0)
lib = pd.DataFrame({s: dict(pool=inv[s]["pcr_pool"], group=inv[s]["pcr_group"], well=inv[s]["plate_well"],
                            segi_reads=int(sum(rc[s].values())), vcf_mean_dp=float(mdp[vix[s]])) for s in CLONE}).T
wp = pairs[pairs.cls1 == "within Wunderlich clone"].copy()
wp["pool_a"] = lib.loc[wp.a, "pool"].values; wp["pool_b"] = lib.loc[wp.b, "pool"].values
wp["group_a"] = lib.loc[wp.a, "group"].values; wp["group_b"] = lib.loc[wp.b, "group"].values
wp["same_pool"] = wp.pool_a == wp.pool_b
wp["same_group"] = wp.group_a == wp.group_b
wp["min_vcf_mean_dp"] = np.minimum(lib.loc[wp.a, "vcf_mean_dp"].values.astype(float), lib.loc[wp.b, "vcf_mean_dp"].values.astype(float))
wp["min_segi_reads"] = np.minimum(lib.loc[wp.a, "segi_reads"].values.astype(float), lib.loc[wp.b, "segi_reads"].values.astype(float))
rho_dp = stats.spearmanr(wp.d_mat, wp.min_vcf_mean_dp)
rho_rd = stats.spearmanr(wp.d_mat, wp.min_segi_reads)
wp[["a", "b", "d_mat", "pool_a", "pool_b", "group_a", "group_b", "same_pool", "same_group", "min_vcf_mean_dp",
    "min_segi_reads"]].sort_values("d_mat").to_csv(f"{OUT}/R5_T2_within_clone.tsv", sep="\t", index=False, float_format="%.4f")
with open(f"{OUT}/R5_T2_within_clone.tsv", "a") as f:
    f.write(f"# Spearman distance ~ min VCF mean DP: rho {rho_dp.statistic:.3f} (p {rho_dp.pvalue:.3g}); ~ min SEGI reads: "
            f"rho {rho_rd.statistic:.3f} (p {rho_rd.pvalue:.3g}); n 36 pairs (not independent: 9 libraries)\n")
    f.write("# libraries: " + "; ".join(f"{s} pool {lib.loc[s,'pool']} group {lib.loc[s,'group']} SEGI {lib.loc[s,'segi_reads']} "
                                         f"VCF meanDP {float(lib.loc[s,'vcf_mean_dp']):.1f}" for s in CLONE) + "\n")
sp = wp[wp.same_pool]
print(f"T2: Spearman d ~ min VCF mean DP rho {rho_dp.statistic:.3f} (p {rho_dp.pvalue:.2g}); ~ min SEGI reads rho "
      f"{rho_rd.statistic:.3f}; same-pool pairs {len(sp)}: " + ", ".join(f"{a}-{b} {d:.4f}" for a, b, d in zip(sp.a, sp.b, sp.d_mat))
      + f"; max pair {wp.loc[wp.d_mat.idxmax(), 'a']}-{wp.loc[wp.d_mat.idxmax(), 'b']}")

# ---- T3 chr4-excluded sensitivity (uncorrected full vs uncorrected no-chr4)
Du, Dn = readmat(INPUTS["ibd_uncorrected"]), readmat(INPUTS["ibd_no_chr4"])
south = sorted(set(pairs.loc[pairs.cls1 == "clone vs other southern trees", "b"]))
sets = {"within clone": [(a, b) for i, a in enumerate(CLONE) for b in CLONE[i + 1:]],
        "clone vs WPR": [(a, WPR) for a in CLONE], "clone vs C1": [(a, C1) for a in CLONE],
        "clone vs other southern trees": [(a, b) for a in CLONE for b in south],
        "technical repeat": [("22_MRD", "28_MRD"), ("16_416-1", "35_416_1")]}
with open(f"{OUT}/R5_T3_no_chr4.tsv", "w") as f:
    f.write("class\tn\tfull_min\tfull_median\tfull_max\tnochr4_min\tnochr4_median\tnochr4_max\tcorrected_min\tcorrected_max\n")
    for k, ps in sets.items():
        a = np.array([Du.loc[x, y] for x, y in ps]); b = np.array([Dn.loc[x, y] for x, y in ps])
        c = np.array([D.loc[x, y] for x, y in ps])
        f.write(f"{k}\t{len(ps)}\t{a.min():.4f}\t{np.median(a):.4f}\t{a.max():.4f}\t{b.min():.4f}\t{np.median(b):.4f}\t"
                f"{b.max():.4f}\t{c.min():.4f}\t{c.max():.4f}\n")
print(open(f"{OUT}/R5_T3_no_chr4.tsv").read())

# ---- T4 southern set
lat = {s: float(inv[s]["lat"]) for s in south}
with open(f"{OUT}/R5_T4_south.tsv", "w") as f:
    f.write("item\tvalue\n")
    f.write(f"parent SOUTH threshold\t{south_parent}\n")
    f.write(f"southern trees in parent's 'clone vs other southern trees'\t{len(south)}\n")
    f.write(f"max lat in set\t{max(lat.values()):.4f}\n")
    f.write(f"switched to 37.82 (no T4 contrast needed)\t{south_parent == 37.82}\n")
print(f"T4: parent SOUTH {south_parent}; {len(south)} southern trees, max lat {max(lat.values()):.4f}; switched: {south_parent == 37.82}")
