# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Read depth along chromosome 4 in 20-Mb bins. Part of the frozen run; NOT reported in the
# current version of the paper (dropped because it had little power).
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
"""R1c-g, R1c (whole vs partial chr4, clone) and R2-364 (14_364 chr4 / chr9) positional scans (decision 2.4, 2.5).

Profile: SEGI 20-Mb bins with >= 50 SNPs (SEGI position = best PAF hit, vcfcache.py). Per sample: summed FORMAT/DP in
the bin / summed DP on chromosomes other than chr4 and the bin's own chromosome; log2 vs the PANEL_ALL bin median.
Clone profile = mean of the 9 ramets' log2.
Noise: sigma_total per bin = 1.4826 MAD over the panel; sigma_tech from TECH (sqrt(mean d^2 / 2));
sigma_genet^2 = sigma_total^2 - sigma_tech^2; clone sigma_b = sqrt(sigma_genet^2 + sigma_tech^2 / 9).
Scan: every contiguous segment of k = 2..n-1 bins: z = (mean_in - mean_out) / (sigma_bar sqrt(1/k + 1/(n-k))),
sigma_bar = median per-bin sigma (clone: sigma_b; single samples: sigma_total). T = max |z|.
Null: T of each PANEL_ALL genet's profile on that chromosome (sigma_total), leaving out the tested sample (DV02).
Outputs: R1c_bins.tsv, R1c_null.tsv, R1c_results.tsv, R2_364_results.tsv; stdout -> 05_R1c_positional.log
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

BIN = 20_000_000
S = sample_sets()
V = vcfcache.load()
samples = list(V["samples"])
ix = {s: i for i, s in enumerate(samples)}
panel = [s for s in S["PANEL_ALL"] if s in ix]
pidx = np.array([ix[s] for s in panel])
print(f"VCF samples {len(samples)}; PANEL_ALL in VCF {len(panel)}; clone in VCF {sum(s in ix for s in CLONE)}")

chrom, pos, DP = V["chrom"], V["pos"], V["DP"].astype(float)
bkey = defaultdict(list)
for j, (c, p) in enumerate(zip(chrom, pos)):
    bkey[(c, int(p) // BIN)].append(j)
keys = sorted((k for k in bkey if len(bkey[k]) >= 50), key=lambda k: (CH.index(k[0]), k[1]))
M = np.array([DP[bkey[k]].sum(0) for k in keys])
nsnp = np.array([len(bkey[k]) for k in keys])
cof = np.array([k[0] for k in keys])
print(f"bins (20 Mb, >= 50 SNPs): {len(keys)}; chr4 {int((cof == 'chr4').sum())}; chr9 {int((cof == 'chr9').sum())}")
tot = {c: M[cof == c].sum(0) for c in CH}
L = np.zeros_like(M)
with np.errstate(divide="ignore"):
    for c in CH:
        r = M[cof == c] / sum(tot[x] for x in CH if x not in (c, "chr4"))
        L[cof == c] = np.log2(r / np.median(r[:, pidx], axis=1, keepdims=True))
mad = lambda x: 1.4826 * np.median(np.abs(x - np.median(x)))
sig_tot = np.array([mad(L[j, pidx]) for j in range(len(keys))])
dd = np.array([[L[j, ix[a]] - L[j, ix[b]] for a, b in TECH] for j in range(len(keys))])
sig_tech = np.sqrt(np.mean(dd ** 2, axis=1) / 2)
sig_gen = np.sqrt(np.clip(sig_tot ** 2 - sig_tech ** 2, 0, None))
sig_b = np.sqrt(sig_gen ** 2 + sig_tech ** 2 / 9)


def scan(x, sig):
    n = len(x)
    cs = np.concatenate([[0], np.cumsum(x)])
    tot_ = cs[-1]
    best, segs = -1.0, []
    for k in range(2, n):
        for s0 in range(0, n - k + 1):
            si = cs[s0 + k] - cs[s0]
            mi, mo = si / k, (tot_ - si) / (n - k)
            z = (mi - mo) / (sig * np.sqrt(1 / k + 1 / (n - k)))
            if abs(z) > best + 1e-9:
                best, segs = abs(z), [(s0, k, mi, mo, z)]
            elif abs(abs(z) - best) <= 1e-9:
                segs.append((s0, k, mi, mo, z))
    return best, segs


def null(c, exclude=()):
    sel = cof == c
    sb = np.median(sig_tot[sel])
    Ts = np.array([scan(L[sel, ix[s]], sb)[0] for s in panel if s not in exclude])
    return Ts, sb


nulls = {}
with open(f"{OUT}/R1c_null.tsv", "w") as f:
    f.write("chrom\tn_bins\tn_null\tsigma_total_med\tsigma_tech_med\tsigma_genet_med\tT95\tT99\tT95_no14_364\tT99_no14_364\n")
    for c in CH:
        Ts, sb = null(c)
        Ts2, _ = null(c, exclude=(S364,))
        nulls[c] = (Ts, sb)
        sel = cof == c
        f.write(f"{c}\t{sel.sum()}\t{len(Ts)}\t{sb:.4f}\t{np.median(sig_tech[sel]):.4f}\t{np.median(sig_gen[sel]):.4f}\t"
                f"{np.quantile(Ts,0.95):.3f}\t{np.quantile(Ts,0.99):.3f}\t{np.quantile(Ts2,0.95):.3f}\t{np.quantile(Ts2,0.99):.3f}\n")
print(open(f"{OUT}/R1c_null.tsv").read())


def q(c, exclude):
    sel = cof == c
    sb = np.median(sig_tot[sel])
    Ts = np.array([scan(L[sel, ix[s]], sb)[0] for s in panel if s not in exclude])
    return np.quantile(Ts, 0.95), np.quantile(Ts, 0.99), len(Ts)


def segtxt(segs, c):
    starts = [k[1] for k in keys if k[0] == c]
    out = []
    for s0, k, mi, mo, z in segs:
        out.append(f"bins {s0}-{s0+k-1} (SEGI {starts[s0]*BIN/1e6:.0f}-{(starts[s0+k-1]+1)*BIN/1e6:.0f} Mb, k={k}): "
                   f"in {mi:+.3f} out {mo:+.3f} delta {mi-mo:+.3f} z {z:+.2f}")
    return "; ".join(out)


res = []
# ---- R1c-g
print("## R1c-g: whole-copy gain controls on their gained chromosome (sigma_total; null leaves the tree out)")
gfail = 0
for t, c in GAINCTRL.items():
    sel = cof == c
    T, segs = scan(L[sel, ix[t]], np.median(sig_tot[sel]))
    t95, t99, nn = q(c, (t,))
    ok = T < t99
    gfail += not ok
    lvl = L[sel, ix[t]].mean()
    print(f"{t} {c}: T {T:.2f} (T95 {t95:.2f}, T99 {t99:.2f}, n {nn}) -> {'flat (T<T99)' if ok else 'FAIL'}; "
          f"mean log2 {lvl:+.3f} (c_hat {6*2**lvl:.2f}); best {segtxt(segs[:1], c)}")
    res.append(("R1c-g", t, c, T, t95, t99, nn, lvl, 6 * 2 ** lvl, "flat" if ok else "fail", segtxt(segs, c)))
gate = gfail < 2
print(f"R1c-g: {3-gfail} of 3 flat -> gate {'PASS' if gate else 'FAIL (R1c and R2-364 descriptive only)'}")

# ---- R1c clone
sel4 = cof == "chr4"
n4 = int(sel4.sum())
prof = np.mean([L[sel4, ix[s]] for s in CLONE], axis=0)
sb4 = np.median(sig_b[sel4])
T, segs = scan(prof, sb4)
Ts, _ = nulls["chr4"]
t95, t99 = np.quantile(Ts, 0.95), np.quantile(Ts, 0.99)
lvl = prof.mean()
f_lo = float([l for l in open(f"{OUT}/R1a_genet.tsv") if l.startswith("f_bound")][0].split("f_lo=")[1].split("\t")[0])
s0, k, mi, mo, z = segs[0]
if T < t95 and f_lo >= 0.80:
    cls = "WHOLE"
elif T >= t99 and k >= 2 and abs(mi - mo) >= 0.13:
    cls = "PARTIAL/STRUCTURED"
else:
    cls = "INCONCLUSIVE"
if not gate:
    cls = f"descriptive only (R1c-g failed); scan pattern would read {cls}"
delta = abs(np.log2(5 / 6))
mds = [kk for kk in range(1, n4) if delta / (sb4 * np.sqrt(1 / kk + 1 / (n4 - kk))) >= 4]
mds_k = min(mds) if mds else None
print(f"\n## R1c clone chr4: {n4} bins; sigma_b median {sb4:.4f}; mean log2 {lvl:+.3f} (c_hat {6*2**lvl:.3f}; expected "
      f"log2(5/6) = -0.263); bins below 0: {(prof < 0).sum()}/{n4}")
print(f"T {T:.2f} vs T95 {t95:.2f} / T99 {t99:.2f} (null n {len(Ts)}); f_lo {f_lo:.3f}; best {segtxt(segs, 'chr4')}")
print(f"MDS (z >= 4, one-copy step 0.263, sigma_b {sb4:.4f}): {mds_k} bins = {mds_k*BIN/1e6:.0f} Mb of {n4*BIN/1e6:.0f} Mb")
print(f"R1c classification: {cls}")
# descriptive power diagnostics (design calculation, not a test): null spread, and the largest expected |z| that a
# one-copy step over the best possible segment (k = n/2) would give, for the clone (sigma_b) and a single sample
kh = n4 // 2
ez_clone = delta / (sb4 * np.sqrt(1 / kh + 1 / (n4 - kh)))
ez_single = delta / (np.median(sig_tot[sel4]) * np.sqrt(1 / kh + 1 / (n4 - kh)))
print(f"chr4 null T: median {np.median(Ts):.2f}, P25 {np.quantile(Ts,0.25):.2f}, P75 {np.quantile(Ts,0.75):.2f}; clone T "
      f"{T:.2f} is at null percentile {100*(Ts < T).mean():.1f}")
print(f"POWER NOTE: max expected z for a one-copy step on half of chr4 = {ez_clone:.2f} (clone, sigma_b) / "
      f"{ez_single:.2f} (single sample, sigma_total), vs T95 {t95:.2f} / T99 {t99:.2f}: the PARTIAL and SEGMENTAL rules "
      f"cannot be met by a one-copy segmental change of any length against this null")
res.append(("R1c-power(descr)", "CLONE/single", "chr4", np.nan, t95, t99, len(Ts), np.nan, np.nan,
            f"null median T {np.median(Ts):.2f}; clone T at null pct {100*(Ts < T).mean():.1f}; max expected z for a "
            f"half-chromosome one-copy step: clone {ez_clone:.2f}, single sample {ez_single:.2f}", ""))
res.append(("R1c", "CLONE(9-ramet mean)", "chr4", T, t95, t99, len(Ts), lvl, 6 * 2 ** lvl, cls, segtxt(segs, "chr4")))
res.append(("R1c-MDS", "CLONE", "chr4", np.nan, np.nan, np.nan, n4, np.nan, np.nan,
            f"MDS {mds_k} bins = {mds_k*BIN/1e6:.0f} Mb", f"sigma_b {sb4:.4f}"))
# per-ramet descriptive T (sigma_total), not a test
for s in CLONE:
    Tr, sg = scan(L[sel4, ix[s]], np.median(sig_tot[sel4]))
    res.append(("R1c-ramet(descr)", s, "chr4", Tr, t95, t99, len(Ts), L[sel4, ix[s]].mean(),
                6 * 2 ** L[sel4, ix[s]].mean(), "descriptive", segtxt(sg[:1], "chr4")))
print("per-ramet T (descriptive): " + ", ".join(f"{r[1]} {r[3]:.2f}" for r in res if r[0] == "R1c-ramet(descr)"))

# ---- R2-364
print("\n## R2-364: 14_364 on chr4 and chr9 (sigma_total; null leaves 14_364 out)")
cls364 = {}
for c, target in [("chr4", 5), ("chr9", 7)]:
    sel = cof == c
    x = L[sel, ix[S364]]
    T, segs = scan(x, np.median(sig_tot[sel]))
    t95, t99, nn = q(c, (S364,))
    lvl = x.mean()
    seg_ok = False
    for s0, k, mi, mo, z in segs:
        for a, b in [(mi, mo), (mo, mi)]:     # DV06: both assignments of a prefix/suffix tie
            if abs(6 * 2 ** a - target) <= 0.5 and abs(6 * 2 ** b - 6) <= 0.5:
                seg_ok = True
    if T >= t99 and seg_ok:
        cc = "SEGMENTAL"
    elif T < t95:
        cc = "UNIFORM-FRACTIONAL"
    else:
        cc = "INCONCLUSIVE"
    if not gate:
        cc = f"descriptive only (R1c-g failed); pattern {cc}"
    cls364[c] = cc
    s0, k, mi, mo, z = segs[0]
    print(f"{c}: {sel.sum()} bins; mean log2 {lvl:+.3f} (c_hat {6*2**lvl:.2f}); T {T:.2f} (T95 {t95:.2f}, T99 {t99:.2f}); "
          f"best {segtxt(segs, c)}; segment c_hat {6*2**mi:.2f}, complement c_hat {6*2**mo:.2f} -> {cc}")
    res.append(("R2-364", S364, c, T, t95, t99, nn, lvl, 6 * 2 ** lvl, cc, segtxt(segs, c)))

with open(f"{OUT}/R1c_results.tsv", "w") as f:
    f.write("item\tsample\tchrom\tT\tT95\tT99\tn_null\tmean_log2\tc_hat_from_mean\tverdict\tbest_segment(s)\n")
    for r in res:
        if r[0] != "R2-364":
            f.write("\t".join(f"{v:.3f}" if isinstance(v, (float, np.floating)) else str(v) for v in r) + "\n")

# ---- 14_364 mixture criteria and library QC
adh = defaultdict(lambda: [0, 0])
with open(INPUTS["ad_het_per_chr"]) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["chrom"] in CH:
            a = adh[short(r["sample"])]; a[0] += int(r["n_het"]); a[1] += int(r["n_total"])
gw = {s: v[0] / v[1] for s, v in adh.items()}
RO, AO = V["RO"], V["AO"]
D = RO + AO
mn = np.minimum(RO, AO)
both = (D >= 24) & (mn >= 1)
sub = both & (mn < D / 12)
subshare = sub.sum(0) / np.maximum(both.sum(0), 1)
pan_gw = np.array([gw[s] for s in panel if s != S364 and s in gw])
pan_sub = np.array([subshare[ix[s]] for s in panel if s != S364])
p90, p95 = np.quantile(pan_gw, 0.90), np.quantile(pan_sub, 0.95)
mix_dis = gw[S364] <= p90 and subshare[ix[S364]] <= p95
print(f"\n14_364 genome-wide AD-het {gw[S364]:.4f} (PANEL_ALL P90 {p90:.4f}, median {np.median(pan_gw):.4f}); "
      f"sub-lattice share {subshare[ix[S364]]:.4f} (P95 {p95:.4f}, median {np.median(pan_sub):.4f}) -> mixture "
      f"{'disfavoured' if mix_dis else 'NOT disfavoured'}")
dup = {}
with open(INPUTS["neale_samtools_stats"]) as fh:
    rr = csv.DictReader(fh, delimiter="\t")
    for r in rr:
        s = short(r["Sample"].split("|")[-1].strip().replace("_neale_samtools", ""))
        dup[s] = float(r["reads_duplicated_percent"])
pd_ = np.array([dup[s] for s in panel if s in dup and s != S364])
segi_reads = {s: sum(S["rc"][s].values()) for s in S["big"]}
pr = np.array([segi_reads[s] for s in panel if s != S364])
inv = S["inv"]
qc = (f"SEGI reads chr1-11 {int(segi_reads[S364])} (panel percentile {100*(pr < segi_reads[S364]).mean():.0f}); "
      f"Neale dup% {dup.get(S364, float('nan')):.2f} (panel median {np.median(pd_):.2f}, percentile "
      f"{100*(pd_ < dup.get(S364, np.nan)).mean():.0f}); pool {inv[S364]['pcr_pool']}; group {inv[S364]['pcr_group']}; "
      f"well {inv[S364]['plate_well']}; PCR date {inv[S364]['pcr_date']}")
print("14_364 library QC: " + qc)
with open(f"{OUT}/R2_364_results.tsv", "w") as f:
    f.write("item\tsample\tchrom\tT\tT95\tT99\tn_null\tmean_log2\tc_hat_from_mean\tverdict\tbest_segment(s)\n")
    for r in res:
        if r[0] == "R2-364":
            f.write("\t".join(f"{v:.3f}" if isinstance(v, (float, np.floating)) else str(v) for v in r) + "\n")
    f.write(f"mixture_check\t{S364}\tgenome\t\t\t\t{len(pan_gw)}\t\t\tAD-het {gw[S364]:.4f} vs P90 {p90:.4f}; sub-lattice "
            f"{subshare[ix[S364]]:.4f} vs P95 {p95:.4f}; mixture {'disfavoured' if mix_dis else 'not disfavoured'}"
            f" (applies to the UNIFORM reading)\t\n")
    f.write(f"library_QC\t{S364}\t\t\t\t\t\t\t\t{qc}\t\n")

# ---- bin table (profiles) for figures
with open(f"{OUT}/R1c_bins.tsv", "w") as f:
    f.write("chrom\tbin_start_Mb\tn_snps\tsigma_total\tsigma_tech\tsigma_genet\tsigma_b\tclone_mean_log2\t"
            + "\t".join(CLONE) + f"\t{S364}\t{WPR}\t{C1}\t" + "\t".join(GAINCTRL) + "\tpanel_q05\tpanel_q95\n")
    for j, k in enumerate(keys):
        vals = [np.mean([L[j, ix[s]] for s in CLONE])] + [L[j, ix[s]] for s in CLONE + [S364, WPR, C1] + list(GAINCTRL)]
        f.write(f"{k[0]}\t{k[1]*BIN/1e6:.0f}\t{nsnp[j]}\t{sig_tot[j]:.4f}\t{sig_tech[j]:.4f}\t{sig_gen[j]:.4f}\t{sig_b[j]:.4f}\t"
                + "\t".join(f"{v:.4f}" for v in vals)
                + f"\t{np.quantile(L[j, pidx], 0.05):.4f}\t{np.quantile(L[j, pidx], 0.95):.4f}\n")
