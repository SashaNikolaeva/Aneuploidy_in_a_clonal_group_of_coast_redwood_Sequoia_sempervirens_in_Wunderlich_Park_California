# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Chromosome-4 copy-number statistic per stem with prediction-interval t, bootstrap CI and
# the bound on the fraction of chromosome 4 lost. Part of the frozen run; these numbers are NOT reported in the
# current version of the paper, which reports alpha_ij (scripts/finalize/09_R0_fig6_table1.py).
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
"""R1a (chr4 copy number per ramet and genet) and R1b (other chromosomes), SEGI direct read counts (decision 2.2-2.3).

Primary: R(s,c) = reads(c) / mean(reads chr1-11 excluding c and chr4); c_hat = 6 R / median(PANEL_ALL).
Genet = mean of the 9 ramets. Prediction-interval t vs PANEL_ALL / LOCAL / SOUTH.
Genet CI: joint bootstrap (resample ramets and panel genets), 10,000 reps, seed 20261026.
f = 6 - c_hat_genet; f_lo = f - 1.96 * 6 * sd_P / median_P.
Secondary: untruncated VCF-DP c_hat (summed per-SNP FORMAT/DP, d02 method) and AD-het ratio (rev05 definition,
het_density_ad_fullpop.tsv).
Outputs: R1a_ramets.tsv, R1a_genet.tsv, R1b_chromosomes.tsv, R1b_WPR_C1.tsv; stdout -> 02_R1ab.log
"""
import csv
import os
import statistics as st
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *  # noqa
import vcfcache  # noqa

S = sample_sets()
rc, inv = S["rc"], S["inv"]
PANELS = {"ALL": S["PANEL_ALL"], "LOCAL": S["PANEL_LOCAL"], "SOUTH": S["PANEL_SOUTH"]}
print("panel sizes:", {k: len(v) for k, v in PANELS.items()})
rng = np.random.default_rng(SEED)

# ---------- secondary sources
V = vcfcache.load()
vs = list(V["samples"])
vix = {s: i for i, s in enumerate(vs)}
tot = {c: V["DP"][V["chrom"] == c].sum(0).astype(float) for c in CH}


def vcf_ratio(s, c):
    i = vix[s]
    return tot[c][i] / sum(tot[k][i] for k in CH if k not in (c, "chr4"))


adhet = defaultdict(dict)
with open(INPUTS["ad_het_per_chr"]) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        adhet[short(r["sample"])][r["chrom"]] = float(r["het_rate"])


def ad_ratio(s, c):
    return adhet[s][c] / st.mean(adhet[s][k] for k in CH if k not in (c, "chr4"))


def genet_stats(fun, c, panel):
    ref = [fun(s, c) for s in panel]
    g = st.mean(fun(s, c) for s in CLONE)
    t, df, pl, pu = pred_t(g, ref)
    return dict(g=g, med=st.median(ref), rel=g / st.median(ref), chat=6 * g / st.median(ref), t=t, df=df,
                p_lower=pl, p_two=2 * min(pl, pu), n=len(ref), mean=st.mean(ref), sd=st.stdev(ref))


vpanel = [s for s in S["PANEL_ALL"] if s in vix]
apanel = [s for s in S["PANEL_ALL"] if s in adhet]
print(f"VCF-DP panel n = {len(vpanel)}; AD-het panel n = {len(apanel)}")

# ---------- R1a per ramet
P = S["PANEL_ALL"]
R4 = {s: segi_ratio(rc, s, "chr4") for s in S["big"]}
pv = np.array([R4[s] for s in P])
med, mu, sd = np.median(pv), pv.mean(), pv.std(ddof=1)
vmed = np.median([vcf_ratio(s, "chr4") for s in vpanel])
amed = np.median([ad_ratio(s, "chr4") for s in apanel])
pool = lambda s: inv[s]["pcr_pool"]
with open(f"{OUT}/R1a_ramets.tsv", "w") as f:
    f.write("sample\tpcr_pool\tpcr_group\tplate_well\tsegi_reads_chr1_11\tR_chr4\tc_hat\tz_vs_ALL\tn_PANEL_ALL_le\t"
            "pool_panel_n\tc_hat_pool_matched\tvcf_dp_c_hat_untrunc\tad_het_rel\tin_[4.5,5.5]\tbelow_all_panel\n")
    print("\n## R1a per ramet (SEGI counts)")
    ok_ramets = True
    for s in CLONE + [WPR, C1]:
        pp = [R4[x] for x in P if pool(x) == pool(s) and x != s]
        refP = pv if s in CLONE else np.array([R4[x] for x in P if x != s])
        m_ = np.median(refP)
        ch = 6 * R4[s] / m_
        z = (R4[s] - refP.mean()) / refP.std(ddof=1)
        nle = int((refP <= R4[s]).sum())
        vch = 6 * vcf_ratio(s, "chr4") / vmed
        arel = ad_ratio(s, "chr4") / amed
        inr = 4.5 <= ch <= 5.5
        below = nle == 0
        if s in CLONE:
            ok_ramets &= inr and below
        f.write(f"{s}\t{inv[s]['pcr_pool']}\t{inv[s]['pcr_group']}\t{inv[s]['plate_well']}\t{int(sum(rc[s].values()))}\t"
                f"{R4[s]:.5f}\t{ch:.3f}\t{z:+.2f}\t{nle}\t{len(pp)}\t{6*R4[s]/np.median(pp):.3f}\t{vch:.3f}\t{arel:.4f}\t"
                f"{inr}\t{below}\n")
        print(f"{s}\tpool {pool(s)}\tc_hat {ch:.3f}\tz {z:+.2f}\tn<= {nle}\tpool-matched {6*R4[s]/np.median(pp):.3f} "
              f"(n {len(pp)})\tVCF-DP {vch:.3f}\tAD-het rel {arel:.3f}")

# ---------- R1a genet
res = {k: genet_stats(lambda s, c: segi_ratio(rc, s, c), "chr4", v) for k, v in PANELS.items()}
g = res["ALL"]["g"]
cl = np.array([R4[s] for s in CLONE])
B = 10000
bi = rng.integers(0, len(cl), (B, len(cl)))
bp = rng.integers(0, len(pv), (B, len(pv)))
boot = 6 * cl[bi].mean(1) / np.median(pv[bp], axis=1)
lo, hi = np.percentile(boot, [2.5, 97.5])
chat = 6 * g / med
fval = 6 - chat
f_lo = fval - 1.96 * 6 * sd / med
vres = genet_stats(vcf_ratio, "chr4", vpanel)
ares = genet_stats(ad_ratio, "chr4", apanel)
crit = {k: res[k]["p_lower"] < 0.001 for k in PANELS}
met = all(crit.values()) and ok_ramets and hi < 5.5
with open(f"{OUT}/R1a_genet.tsv", "w") as f:
    f.write("metric\tpanel\tn_panel\tgenet_value\tpanel_median\trel\tc_hat\tci95_lo\tci95_hi\tpred_t\tdf\tp_one_sided_lower\n")
    for k in PANELS:
        r = res[k]
        cil, cih = (f"{lo:.3f}", f"{hi:.3f}") if k == "ALL" else ("", "")
        f.write(f"SEGI_counts\t{k}\t{r['n']}\t{r['g']:.5f}\t{r['med']:.5f}\t{r['rel']:.4f}\t{r['chat']:.3f}\t{cil}\t{cih}\t"
                f"{r['t']:.2f}\t{r['df']}\t{r['p_lower']:.3g}\n")
    f.write(f"VCF_DP_untruncated\tALL\t{vres['n']}\t{vres['g']:.5f}\t{vres['med']:.5f}\t{vres['rel']:.4f}\t{vres['chat']:.3f}"
            f"\t\t\t{vres['t']:.2f}\t{vres['df']}\t{vres['p_lower']:.3g}\n")
    f.write(f"AD_het_ratio\tALL\t{ares['n']}\t{ares['g']:.5f}\t{ares['med']:.5f}\t{ares['rel']:.4f}\t\t\t\t{ares['t']:.2f}"
            f"\t{ares['df']}\t{ares['p_lower']:.3g}\n")
    f.write(f"f_bound\tALL\t{len(pv)}\tf={fval:.3f}\tsd_P={sd:.5f}\t\tf_lo={f_lo:.3f}\t\t\t\t\t\n")
    f.write(f"R1a_threshold\t\t\tp<0.001 all panels: {crit}; ramets in [4.5,5.5] and below every panel genet: "
            f"{ok_ramets}; CI upper < 5.5: {hi < 5.5}\tMET={met}\t\t\t\t\t\t\t\n")
print(f"\n## R1a genet: c_hat {chat:.3f} (95% CI {lo:.3f}-{hi:.3f}, joint bootstrap B={B}); "
      + "; ".join(f"{k}: t {res[k]['t']:.2f} df {res[k]['df']} p {res[k]['p_lower']:.2g}" for k in PANELS))
print(f"f = {fval:.3f}; f_lo = {f_lo:.3f} (sd_P {sd:.5f}, median_P {med:.5f})")
print(f"secondary: VCF-DP untruncated genet c_hat {vres['chat']:.3f} (t {vres['t']:.2f}, n {vres['n']}); "
      f"AD-het rel {ares['rel']:.4f} (t {ares['t']:.2f})")
print(f"R1a threshold met: {met}")

# ---------- R1b
rows = []
for c in [x for x in CH if x != "chr4"]:
    r = {k: genet_stats(lambda s, cc: segi_ratio(rc, s, cc), c, v) for k, v in PANELS.items()}
    rv = genet_stats(vcf_ratio, c, vpanel)
    ra = genet_stats(ad_ratio, c, apanel)
    rows.append((c, r, rv, ra))
adj = holm([r["ALL"]["p_two"] for _, r, _, _ in rows])
with open(f"{OUT}/R1b_chromosomes.tsv", "w") as f:
    f.write("chrom\tc_hat_ALL\tt_ALL\tp_two_ALL\tholm_ALL_m10\tc_hat_LOCAL\tt_LOCAL\tp_two_LOCAL\tc_hat_SOUTH\tt_SOUTH\t"
            "p_two_SOUTH\tdeviates\tvcf_dp_c_hat\tvcf_dp_t\tad_het_rel\tad_het_t\n")
    print("\n## R1b genet, other chromosomes (SEGI counts; two-sided; Holm m=10 vs ALL)")
    for (c, r, rv, ra), h in zip(rows, adj):
        dev = (h < 0.05 and np.sign(r["LOCAL"]["t"]) == np.sign(r["ALL"]["t"]) and r["LOCAL"]["p_two"] < 0.05
               and abs(r["ALL"]["chat"] - 6) >= 0.5)
        f.write(f"{c}\t{r['ALL']['chat']:.3f}\t{r['ALL']['t']:.2f}\t{r['ALL']['p_two']:.3g}\t{h:.3g}\t{r['LOCAL']['chat']:.3f}\t"
                f"{r['LOCAL']['t']:.2f}\t{r['LOCAL']['p_two']:.3g}\t{r['SOUTH']['chat']:.3f}\t{r['SOUTH']['t']:.2f}\t"
                f"{r['SOUTH']['p_two']:.3g}\t{dev}\t{rv['chat']:.3f}\t{rv['t']:.2f}\t{ra['rel']:.4f}\t{ra['t']:.2f}\n")
        print(f"{c}\tc_hat {r['ALL']['chat']:.3f}\tt {r['ALL']['t']:+.2f}\tp {r['ALL']['p_two']:.2g}\tHolm {h:.2g}\t"
              f"LOCAL t {r['LOCAL']['t']:+.2f}\tdeviates {dev}\tVCF-DP {rv['chat']:.3f}\tAD-het {ra['rel']:.3f} (t {ra['t']:+.2f})")
    print(f"chromosomes deviating: {sum(1 for (c, r, rv, ra), h in zip(rows, adj) if h < 0.05)} with Holm < 0.05")

# ---------- WPR and C1 rows (not tested; leave-one-out panel)
with open(f"{OUT}/R1b_WPR_C1.tsv", "w") as f:
    f.write("sample\tchrom\tc_hat\tz_vs_ALL_loo\n")
    for s in [WPR, C1]:
        out = []
        for c in CH:
            ref = np.array([segi_ratio(rc, x, c) for x in P if x != s])
            v = segi_ratio(rc, s, c)
            ch, z = 6 * v / np.median(ref), (v - ref.mean()) / ref.std(ddof=1)
            f.write(f"{s}\t{c}\t{ch:.3f}\t{z:+.2f}\n")
            out.append(f"{c} {ch:.2f} ({z:+.1f})")
        print(f"{s}: " + "; ".join(out))
