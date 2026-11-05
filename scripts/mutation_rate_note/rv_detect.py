# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Probability of detecting a stem-private mutation under the search filters at the stems'
# depths (Supplementary Note S1: about 85% at fraction 1/6, 28% at 1/12).
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
"""Reviewer check: detection probability of a stem-private mutation under the scan's real carrier filters
(alt reads >=2, alt fraction >=0.10, not >95% of alt reads on one strand), at fixed depth and averaged over the
chr1 depth distribution of callable stems (rv_calib_chr1_depth.tsv). Strands assumed 50:50 and binomial sampling;
ignores the 'zero alt reads in other stems' condition (true mutations pass it unless an error coincides).
Also recomputes the illustrative expected counts (draft Table S2) and adds the published per-year rates."""
import csv, math
from collections import defaultdict
OUT = "/space/s1/sashanikolaeva/Scripts/GitHub/redwood-clone-paper/analysis/mutation_rate_note"
def binom(n, k, p): return math.comb(n, k) * p**k * (1 - p)**(n - k)
def p_detect(dp, f, strand=True):
    tot = 0.0
    for k in range(2, dp + 1):
        if k / dp < 0.10: continue
        ps = 1.0
        if strand:  # fail if max(fwd, rev)/k > 0.95
            ps = sum(binom(k, j, 0.5) for j in range(k + 1) if max(j, k - j) / k <= 0.95)
        tot += binom(dp, k, f) * ps
    return tot
print("fixed depth: dp, f, P(>=2 alt only), P(all carrier filters)")
for dp in (10, 15, 20, 30, 40):
    for f in (1/6, 1/12):
        print(dp, f"{f:.3f}", round(1 - binom(dp, 0, f) - binom(dp, 1, f), 3), round(p_detect(dp, f), 3))
hist = defaultdict(dict)
for r in csv.DictReader(open(f"{OUT}/rv_calib_chr1_depth.tsv"), delimiter="\t"):
    hist[r["stem"]][int(r["depth_capped200"])] = int(r["n_sites"])
print("\nchr1 callable sites: per stem median depth and depth-averaged detection (all carrier filters)")
allh = defaultdict(int)
for s, h in hist.items():
    n = sum(h.values()); cum = 0
    for d in sorted(h):
        cum += h[d]
        if cum >= n / 2: med = d; break
    for d, c in h.items(): allh[d] += c
    pa = {f: sum(c * p_detect(d, f) for d, c in h.items()) / n for f in (1/6, 1/12)}
    print(s, "median DP", med, "share DP<15", round(sum(c for d, c in h.items() if d < 15) / n, 3), "P(1/6)", round(pa[1/6], 3), "P(1/12)", round(pa[1/12], 3))
n = sum(allh.values())
print("ALL stems pooled: P(1/6)", round(sum(c * p_detect(d, 1/6) for d, c in allh.items()) / n, 3),
      "P(1/12)", round(sum(c * p_detect(d, 1/12) for d, c in allh.items()) / n, 3))
L = 20.4e6
print("\nexpected new mutations per stem lineage = rate x years x 6 copies x 20.4 Mb (every one detected)")
for lab, mu in (("Hanlon2019 7.4e-11", 7.4e-11), ("1e-10", 1e-10), ("Hofmeister2020 1.33e-10", 1.33e-10), ("Satake2024 8e-10", 8e-10),
                ("1e-9", 1e-9), ("Orr2020 upper 1.12e-9", 1.12e-9), ("Moore2023 min 6.9e-9", 6.9e-9), ("1e-8", 1e-8),
                ("Moore2023 3.4e-8", 3.4e-8), ("Moore2023 max 8.29e-7", 8.29e-7)):
    print(f"{lab:26s} 78 yr: {mu*78*6*L:9.1f}   177 yr: {mu*177*6*L:9.1f}   (per-site rate, no x6: {mu*78*L:8.1f} / {mu*177*L:8.1f})")
