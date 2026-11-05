# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Replication of stem-private alleles (kappa) in the other stems versus same-tree comparison
# pairs (Methods: Replication of stem-private alleles; Results; Supplementary Note S1).
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
"""R4 G000 somatic non-replication test, no rate (decision 2.7).

Input: het_cause 05_replicate_alleles.tsv.gz = single-library alleles (lib1 >= 3 alt reads, no other genet >= 2 alt
reads, lib2 depth >= 8); observed = alt_reads_lib2 >= 2; expected = expected_P_ge2_alt_lib2.
Groups: G000 (all ordered ramet pairs); comparator = same-tree / same-tissue genets G002 G005 G006 G007 G008 G009 G010
G013; sensitivity row = all 13 non-G000 genets.
kappa = sum obs / sum exp; Delta = kappa_comp - kappa_G000; E = (1 - kappa_G000/kappa_comp) x single-library alleles
per G000 library. 95% CI: two-level bootstrap, 2,000 reps (G000 resamples ramets keeping pairs among distinct ramets;
comparator resamples genets), seed 20261026. Verdict order: DEVIATIONS DV07.
Outputs: R4_results.tsv, R4_decomposition.tsv; stdout -> 07_R4_G000.log
"""
import gzip
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *  # noqa

COMP = ["G002", "G005", "G006", "G007", "G008", "G009", "G010", "G013"]
if SET == "popgen224":
    COMP = [g for g in COMP if g != "G007"]   # G007 = tree 286 (13_286_1 / 3_286-1), not a replicate pair under DV14
EXCL = ["G001", "G003", "G004", "G012", "G011"]
POOL31 = {"31_c3", "31_WP2_1"}
rng = np.random.default_rng(SEED)
d = pd.read_csv(INPUTS["rep_alleles"], sep="\t")
d["obs"] = (d.alt_reads_lib2 >= 2).astype(int)
d["exp"] = d.expected_P_ge2_alt_lib2
d["af1"] = d.alt_reads_lib1 / d.depth_lib1
g0 = d[d.genet == "G000"].copy()
cp = d[d.genet.isin(COMP)].copy()
al = d[(d.genet != "G000") & ~((SET == "popgen224") & (d.genet == "G007"))].copy()
ram = sorted(set(g0.lib1) | set(g0.lib2))
assert set(ram) == set(CLONE), ram
npairs = g0[["lib1", "lib2"]].drop_duplicates().shape[0]
print(f"G000: {len(g0)} allele tests, {npairs} ordered pairs among {len(ram)} ramets; comparator {cp.genet.nunique()} genets, "
      f"{len(cp)} tests; all-13 {al.genet.nunique()} genets, {len(al)} tests")
conc = pd.read_csv(INPUTS["rep_concord"], sep="\t")
geo = conc[conc.level == "per_pair"].drop_duplicates("genet").set_index("genet").geo_km
print("comparator geo distance (km):", ", ".join(f"{g} {geo.get(g, np.nan):.3f}" for g in COMP))
print("excluded genets geo distance (km):", ", ".join(f"{g} {geo.get(g, np.nan)}" for g in EXCL))

kap = lambda x: x.obs.sum() / x.exp.sum()
nlib = g0.drop_duplicates(["lib1", "chrom", "pos", "alt"]).groupby("lib1").size()
per_lib = nlib.reindex(ram).fillna(0)
A = per_lib.mean()
k0, kc, ka = kap(g0), kap(cp), kap(al)
D = kc - k0
E = (1 - k0 / kc) * A
print(f"\nkappa G000 {k0:.4f} (obs {g0.obs.sum()} / exp {g0.exp.sum():.1f}); comparator {kc:.4f} (obs {cp.obs.sum()} / exp "
      f"{cp.exp.sum():.1f}); all-13 {ka:.4f}")
print(f"single-library alleles per G000 library: mean {A:.1f} (range {per_lib.min():.0f}-{per_lib.max():.0f})")

# bootstrap
pair_o = g0.groupby(["lib1", "lib2"]).obs.sum()
pair_e = g0.groupby(["lib1", "lib2"]).exp.sum()
ri = {r: i for i, r in enumerate(ram)}
O = np.zeros((9, 9)); Ex = np.zeros((9, 9))
for (a, b), v in pair_o.items():
    O[ri[a], ri[b]] = v
for (a, b), v in pair_e.items():
    Ex[ri[a], ri[b]] = v
go = cp.groupby("genet").obs.sum().reindex(COMP).values
ge = cp.groupby("genet").exp.sum().reindex(COMP).values
nl = per_lib.values
B = 2000
bs = []
for _ in range(B):
    w = np.bincount(rng.integers(0, 9, 9), minlength=9).astype(float)
    W = np.outer(w, w)
    np.fill_diagonal(W, 0)
    kg = (W * O).sum() / (W * Ex).sum()
    wc = np.bincount(rng.integers(0, len(COMP), len(COMP)), minlength=len(COMP))
    kcb = (wc * go).sum() / (wc * ge).sum()
    Ab = (w * nl).sum() / w.sum()
    bs.append((kg, kcb, kcb - kg, (1 - kg / kcb) * Ab))
bs = np.array(bs)
ci = np.percentile(bs, [2.5, 97.5], axis=0)
print(f"Delta kappa {D:+.4f} (95% CI {ci[0,2]:+.4f} to {ci[1,2]:+.4f}); E {E:.1f} alleles per library "
      f"(95% CI {ci[0,3]:.1f} to {ci[1,3]:.1f}); kappa G000 CI {ci[0,0]:.4f}-{ci[1,0]:.4f}; comparator CI {ci[0,1]:.4f}-{ci[1,1]:.4f}")

# decompositions
dec = []
g0n = g0[~g0.lib1.isin(POOL31) & ~g0.lib2.isin(POOL31)]
k0n = kap(g0n)
Dn = kc - k0n
dec.append(("G000 without pool-31 libraries", len(g0n), k0n, kc, Dn, Dn / D if D else np.nan))
both = pd.concat([g0, cp])
cuts = np.quantile(both.depth_lib2, [1 / 3, 2 / 3])
both["tert"] = np.digitize(both.depth_lib2, cuts)
for t in range(3):
    x0 = both[(both.tert == t) & (both.genet == "G000")]
    xc = both[(both.tert == t) & (both.genet != "G000")]
    k0t, kct = kap(x0), kap(xc)
    dec.append((f"lib2 depth tertile {t+1} ({'<' if t == 0 else '>='}{cuts[0] if t < 2 else cuts[1]:.0f}"
                f"{'' if t != 1 else f' to <{cuts[1]:.0f}'})", len(x0), k0t, kct, kct - k0t, (kct - k0t) / D if D else np.nan))
dec.append(("sensitivity: comparator = all 13 genets", len(g0), k0, ka, ka - k0, (ka - k0) / D if D else np.nan))
print("\n## decompositions (Delta = kappa_comp - kappa_G000; retained = Delta_sub / Delta_all)")
for r in dec:
    print(f"{r[0]}: n_G000 {r[1]}, kappa_G000 {r[2]:.4f}, kappa_comp {r[3]:.4f}, Delta {r[4]:+.4f}, retained {r[5]:.2f}")
print("\n## lib1 alt-fraction spectrum (share with af1 <= 0.20)")
spec = []
for lab, x in [("G000", g0), ("comparator", cp)]:
    for o, nm in [(0, "non-replicating"), (1, "replicating")]:
        y = x[x.obs == o]
        spec.append((lab, nm, len(y), (y.af1 <= 0.20).mean()))
        print(f"{lab} {nm}: n {len(y)}, share af1 <= 0.20 = {(y.af1 <= 0.20).mean():.3f}")

ret31 = dec[0][5]
rettop = dec[3][5]
if ci[0, 2] <= 0 <= ci[1, 2] or ci[0, 2] < 0:
    verdict = f"no excess beyond library replicates; upper bound E_hi = {ci[1,3]:.1f} alleles per library"
elif ret31 < 0.5:
    verdict = "attributable to pool-31 libraries"
elif rettop >= 0.5:
    verdict = f"small stem-specific excess: E = {E:.1f} ({ci[0,3]:.1f}-{ci[1,3]:.1f}) alleles per library at survey SNP sites"
else:
    verdict = "excess not robust to depth (no PREREG label)"
print(f"\nR4 verdict: {verdict}")
with open(f"{OUT}/R4_results.tsv", "w") as f:
    f.write("quantity\tpoint\tci95_lo\tci95_hi\tnote\n")
    f.write(f"kappa_G000\t{k0:.4f}\t{ci[0,0]:.4f}\t{ci[1,0]:.4f}\t{len(g0)} tests, {npairs} ordered pairs\n")
    f.write(f"kappa_comparator_8\t{kc:.4f}\t{ci[0,1]:.4f}\t{ci[1,1]:.4f}\t{','.join(COMP)}; {len(cp)} tests\n")
    f.write(f"kappa_all13\t{ka:.4f}\t\t\tsensitivity\n")
    f.write(f"delta_kappa\t{D:+.4f}\t{ci[0,2]:+.4f}\t{ci[1,2]:+.4f}\tcomparator - G000\n")
    f.write(f"E_alleles_per_library\t{E:.2f}\t{ci[0,3]:.2f}\t{ci[1,3]:.2f}\tsingle-library alleles per G000 library {A:.1f}\n")
    f.write(f"retained_without_pool31\t{ret31:.3f}\t\t\t\n")
    f.write(f"retained_top_depth_tertile\t{rettop:.3f}\t\t\t\n")
    f.write(f"VERDICT\t\t\t\t{verdict}\n")
with open(f"{OUT}/R4_decomposition.tsv", "w") as f:
    f.write("subset\tn_G000_tests\tkappa_G000\tkappa_comp\tdelta\tretained_share\n")
    for r in dec:
        f.write(f"{r[0]}\t{r[1]}\t{r[2]:.4f}\t{r[3]:.4f}\t{r[4]:+.4f}\t{r[5]:.3f}\n")
    for r in spec:
        f.write(f"spectrum {r[0]} {r[1]}\t{r[2]}\tshare_af1<=0.20={r[3]:.3f}\t\t\t\n")
