# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Allele balance on chromosome 4 relative to the other chromosomes, against depth-matched
# survey trees (Methods and Results: allele fractions; allele-share ratio in Supplementary Table S1).
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
"""R3 allele balance, difference-in-differences (decision 2.6; d03 definitions).

Per sample x SEGI chromosome, from FORMAT/RO,AO at the 207,546 SEGI-mapped survey SNPs (vcfcache):
  MAF (primary) = mean min(ref, alt)/DP over sites with DP >= 20 and min(ref, alt) >= 3
  HET (secondary) = share of DP >= 10 sites with min(ref, alt) >= 3
  ratio = chr4 / mean(chr1-3, chr5-11)
Reference: depth-matched panel = PANEL_ALL with chr4 RO and AO binomially thinned to 5/6, mean ratio over 20 draws
(seed 20261026). Un-thinned PANEL_ALL and PANEL_LOCAL alongside.
Test: clone genet ratio (mean of 9 ramets) vs depth-matched panel, prediction t, one-sided (MAF up, HET down);
Holm over m = 2.
Outputs: R3_per_sample.tsv, R3_results.tsv; stdout -> 06_R3_allelic.log
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *  # noqa
import vcfcache  # noqa

S = sample_sets()
V = vcfcache.load()
samples = list(V["samples"])
ix = {s: i for i, s in enumerate(samples)}
P = [s for s in S["PANEL_ALL"] if s in ix]
LOC = [s for s in S["PANEL_LOCAL"] if s in ix]
chrom = V["chrom"]
RO, AO = V["RO"], V["AO"]
rng = np.random.default_rng(SEED)


def stats_chr(ro, ao):
    d = ro + ao
    mn = np.minimum(ro, ao)
    m10 = d >= 10
    h = m10 & (mn >= 3)
    m20 = h & (d >= 20)
    het = h.sum(0) / np.maximum(m10.sum(0), 1)
    maf = np.where(m20, mn / np.maximum(d, 1), 0).sum(0) / np.maximum(m20.sum(0), 1)
    return het, maf


per = {}
for c in CH:
    sel = chrom == c
    per[c] = stats_chr(RO[sel], AO[sel])


def ratios(per_):
    out = {}
    for j, name in enumerate(["HET", "MAF"]):
        other = np.mean([per_[c][j] for c in CH if c != "chr4"], axis=0)
        out[name] = per_["chr4"][j] / other
    return out


R0 = ratios(per)
# thinned chr4 for panel samples, 20 draws
pidx = np.array([ix[s] for s in P])
sel4 = chrom == "chr4"
ro4, ao4 = RO[sel4][:, pidx], AO[sel4][:, pidx]
thin = {"HET": [], "MAF": []}
for b in range(20):
    tr = rng.binomial(ro4, 5 / 6)
    ta = rng.binomial(ao4, 5 / 6)
    h4, m4 = stats_chr(tr, ta)
    for j, name, v4 in [(0, "HET", h4), (1, "MAF", m4)]:
        other = np.mean([per[c][j][pidx] for c in CH if c != "chr4"], axis=0)
        thin[name].append(v4 / other)
thin = {k: np.mean(v, axis=0) for k, v in thin.items()}

with open(f"{OUT}/R3_per_sample.tsv", "w") as f:
    f.write("sample\tset\tHET_ratio\tMAF_ratio\tHET_ratio_thinned_mean20\tMAF_ratio_thinned_mean20\n")
    for s in CLONE:
        f.write(f"{s}\tCLONE\t{R0['HET'][ix[s]]:.5f}\t{R0['MAF'][ix[s]]:.5f}\t\t\n")
    for k, s in enumerate(P):
        f.write(f"{s}\tPANEL_ALL\t{R0['HET'][ix[s]]:.5f}\t{R0['MAF'][ix[s]]:.5f}\t{thin['HET'][k]:.5f}\t{thin['MAF'][k]:.5f}\n")

EXPECT = {"MAF": "+3.7% (k=6) / +11.9% (k=3) / +0.0% (k=2)", "HET": "-4.5% (k=6) / -16.0% (k=3) / -33.6% (k=2)"}
rows = []
print(f"panel n {len(P)} (depth-matched, 20 thinning draws); LOCAL n {len(LOC)}")
for name, direction in [("MAF", "up"), ("HET", "down")]:
    g = np.mean([R0[name][ix[s]] for s in CLONE])
    ram = [R0[name][ix[s]] for s in CLONE]
    refs = {"depth-matched": thin[name], "unthinned": np.array([R0[name][ix[s]] for s in P]),
            "LOCAL_unthinned": np.array([R0[name][ix[s]] for s in LOC])}
    for rn, ref in refs.items():
        t, df, pl, pu = pred_t(g, ref)
        p1 = pu if direction == "up" else pl
        eff = 100 * (g / np.median(ref) - 1)
        rows.append([name, rn, len(ref), g, np.median(ref), eff, t, df, p1, direction])
        print(f"{name} vs {rn}: genet ratio {g:.4f}, panel median {np.median(ref):.4f}, effect {eff:+.2f}%, t {t:+.2f} "
              f"(df {df}), one-sided p ({direction}) {p1:.3g}")
    print(f"  ramet ratios {name}: " + ", ".join(f"{v:.3f}" for v in ram) + f"; expected allelic-only effect {EXPECT[name]}")
prim = [r for r in rows if r[1] == "depth-matched"]
adj = holm([r[8] for r in prim])
indir = [(r[5] > 0) if r[0] == "MAF" else (r[5] < 0) for r in prim]
opp_sig = [(not d) and (r[8] > 0.5) and (1 - r[8]) < 0.05 for r, d in zip(prim, indir)]
# opposite direction with nominal p < 0.05: one-sided p in the other direction = 1 - p1
if any(opp_sig):
    verdict = "not supported"
elif any(a < 0.05 for a in adj) and all(indir):
    verdict = "allelic data independently support a copy loss"
elif all(indir):
    verdict = "consistent, not independently significant"
else:
    verdict = "mixed directions, neither significant (no PREREG label)"
print(f"\nHolm (m=2) on depth-matched: MAF {adj[0]:.3g}, HET {adj[1]:.3g}; both in predicted direction: {all(indir)}")
print(f"R3 verdict: {verdict}")
with open(f"{OUT}/R3_results.tsv", "w") as f:
    f.write("statistic\treference\tn_ref\tclone_genet_ratio\tref_median\teffect_pct\tpred_t\tdf\tp_one_sided\tpredicted_direction\t"
            "holm_m2\texpected_allelic_only_effect\n")
    for r in rows:
        h = f"{adj[0]:.3g}" if (r[1] == "depth-matched" and r[0] == "MAF") else (f"{adj[1]:.3g}" if r[1] == "depth-matched" else "")
        f.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]:.5f}\t{r[4]:.5f}\t{r[5]:+.2f}\t{r[6]:.2f}\t{r[7]}\t{r[8]:.3g}\t{r[9]}\t{h}\t"
                f"{EXPECT[r[0]]}\n")
    f.write(f"VERDICT\t\t\t\t\t\t\t\t\t\t\t{verdict}\n")
