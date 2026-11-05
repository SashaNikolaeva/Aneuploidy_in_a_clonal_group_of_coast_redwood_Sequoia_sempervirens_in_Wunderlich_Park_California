#!/usr/bin/env python3
# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Allele-dosage test: chromosome 4 versus the other chromosomes, block jackknife, 10,000 PCR-
# group-matched random sets of nine trees, and its sensitivity (Methods and Results: allele fractions).
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
"""02: A9 primary allele-dosage test + A4c control + A2p power + A6 positional + A8 sensitivity (PREREG_v2).

D2 per site = (f - p_hat)^2 / (p_hat (1 - p_hat)) - 1/n (f, n: pool; p_hat: reference R minus the pool's own trees);
T = (mean D2 chr4 - mean D2 other ten) / SE, SE by delete-one-block jackknife over 25-Mb SEGI blocks.
Null: 10,000 PCR-group-matched pseudo-clones of 9 R trees. p = (1 + #{T_pseudo >= T_clone}) / 10,001.
Outputs: RUN/02_A9_*.tsv, RUN/02_A9_null.npy, RUN/02_A9_summary.json; stdout -> RUN/02_A9_dosage.log
"""
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dc_common as c  # noqa: E402

arr = c.load_array()
inv = c.load_inventory()
ix, P224, R = c.sets(arr["samples"], inv)
S = np.load(f"{c.RUN}/sites.npz")
meta = c.site_meta(arr)
COMP = {1: 2, 2: 1, 3: 3, 4: 1, 5: 2}
assert sorted(c.pcr_group(s) for s in c.CLONE) == sorted(g for g, k in COMP.items() for _ in range(k))


class Data:
    def __init__(self, mask):
        self.chrom = arr["chrom"][mask]
        self.pos = arr["pos"][mask]
        self.nloci = meta["nloci"][mask]
        self.RO = arr["RO"][mask].astype(np.float64)
        self.AO = arr["AO"][mask].astype(np.float64)
        ri = np.array([ix[s] for s in R])
        self.ri = ri
        self.AOr = self.AO[:, ri].sum(1)
        self.Dr = (self.AO[:, ri] + self.RO[:, ri]).sum(1)
        blk = np.char.add(self.chrom, (self.pos // 25_000_000).astype(str))
        self.blocks, self.bidx = np.unique(blk, return_inverse=True)
        self.Bmat = np.zeros((len(self.chrom), len(self.blocks)))
        self.Bmat[np.arange(len(self.chrom)), self.bidx] = 1.0

    def d2(self, ao, d, ao_own, d_own):
        """ao, d: pool sums (pools x sites); own: the pool's contribution to R (0 for the clone)."""
        ph = (self.AOr - ao_own) / np.maximum(self.Dr - d_own, 1)
        ph = np.clip(ph, 1e-3, 1 - 1e-3)
        f = ao / np.maximum(d, 1)
        x = (f - ph) ** 2 / (ph * (1 - ph)) - 1.0 / np.maximum(d, 1)
        use = d >= 30
        return np.where(use, x, 0.0), use.astype(float)


def contrast(D2, U, ma, mb, Bmat):
    """T for groups a vs b (site masks ma, mb) per pool row; jackknife over blocks."""
    Sa, Na = (D2 * ma) @ Bmat, (U * ma) @ Bmat
    Sb, Nb = (D2 * mb) @ Bmat, (U * mb) @ Bmat
    sa, na, sb, nb = Sa.sum(1, keepdims=True), Na.sum(1, keepdims=True), Sb.sum(1, keepdims=True), Nb.sum(1, keepdims=True)
    theta = (sa / na - sb / nb)[:, 0]
    with np.errstate(invalid="ignore", divide="ignore"):
        th_b = (sa - Sa) / (na - Na) - (sb - Sb) / (nb - Nb)
    occ = (Na + Nb) > 0
    th_b = np.where(occ, th_b, np.nan)
    B = occ.sum(1)
    mean_b = np.nanmean(th_b, axis=1)
    se = np.sqrt((B - 1) / B * np.nansum((th_b - mean_b[:, None]) ** 2, axis=1))
    return theta / se, theta, se


def clone_pool(D, stems):
    ci = [ix[s] for s in stems]
    ao = D.AO[:, ci].sum(1)[None, :]
    d = (D.AO[:, ci] + D.RO[:, ci]).sum(1)[None, :]
    return ao, d


def pseudo_T(D, n, seed, ma, mb):
    rng = np.random.default_rng(seed)
    byg = {g: [ix[s] for s in R if c.pcr_group(s) == g] for g in COMP}
    out = []
    for start in range(0, n, 1000):
        k = min(1000, n - start)
        idx = np.array([np.concatenate([rng.choice(byg[g], COMP[g], replace=False) for g in COMP]) for _ in range(k)])
        ao = D.AO[:, idx].sum(2).T
        d = (D.AO[:, idx] + D.RO[:, idx]).sum(2).T
        x, u = D.d2(ao, d, ao, d)
        out.append(contrast(x, u, ma, mb, D.Bmat)[0])
    return np.concatenate(out)


# ------------------------------------------------------------------ primary
D = Data(S["A9"])
is4, oth = D.chrom == "chr4", D.chrom != "chr4"
ao, d = clone_pool(D, c.CLONE)
x, u = D.d2(ao, d, 0, 0)
T4, th4, se4 = [v[0] for v in contrast(x, u, is4, oth, D.Bmat)]
used = u[0].astype(bool)
print(f"A9 sites {len(D.chrom)} ({is4.sum()} chr4); used by the clone (pool depth >= 30): {used.sum()} ({(used & is4).sum()} chr4); "
      f"jackknife blocks {len(D.blocks)}")
per_chr = pd.DataFrame([dict(chrom=ch, sites_used=int((used & (D.chrom == ch)).sum()),
                             mean_D2=float(x[0][used & (D.chrom == ch)].mean()),
                             median_pool_depth=float(np.median(d[0][used & (D.chrom == ch)])),
                             mean_signed_dev=float(np.mean((ao[0] / d[0] - np.clip(D.AOr / D.Dr, 1e-3, 1 - 1e-3))[used & (D.chrom == ch)])))
                        for ch in c.CH])
null = pseudo_T(D, 10_000, c.SEED, is4, oth)
np.save(f"{c.RUN}/02_A9_null.npy", null)
p = (1 + np.sum(null >= T4)) / (1 + len(null))
q999 = float(np.quantile(null, 0.999))
print(f"clone: theta {th4:+.4f}, SE {se4:.4f}, T {T4:+.3f}; pseudo-clone null: mean {null.mean():+.3f}, SD {null.std():.3f}, "
      f"99.9th pct {q999:.3f}; one-sided p = {p:.5f}")

# A4c: genotype-matched control on the clone's other chromosomes
Tc = {}
for ch in c.CH:
    if ch == "chr4":
        continue
    mc = D.chrom == ch
    mo = (D.chrom != ch) & (D.chrom != "chr4")
    Tc[ch] = float(contrast(x, u, mc, mo, D.Bmat)[0][0])
per_chr["T_c_control"] = per_chr.chrom.map(Tc)
per_chr.to_csv(f"{c.RUN}/02_A9_per_chromosome.tsv", sep="\t", index=False, float_format="%.5f")
print(per_chr.to_string(index=False))
maxTc = max(Tc.values())
q90Tc = float(np.quantile(list(Tc.values()), 0.9))
print(f"A4c: max T_c over the 10 other chromosomes = {maxTc:+.3f}; 90th pct {q90Tc:+.3f}; T_chr4 > max: {T4 > maxTc}")

# A2p power: spikes into the clone's own other chromosomes
rows = []
for pi_ in (0.25, 0.33, 0.50):
    for de in (1 / 15, 2 / 15, 1 / 6):
        hits, vals = 0, []
        for k_c, ch in enumerate([h for h in c.CH if h != "chr4"]):
            mc = D.chrom == ch
            mo = (D.chrom != ch) & (D.chrom != "chr4")
            site_ix = np.where(mc & used)[0]
            rng = np.random.default_rng(c.SEED + 1000 * k_c + int(round(pi_ * 100)) * 10 + int(round(de * 30)))
            AOs = np.repeat(ao, 100, axis=0)
            for r_ in range(100):
                pick = rng.choice(site_ix, int(round(pi_ * len(site_ix))), replace=False)
                sign = rng.choice([-1.0, 1.0], len(pick))
                AOs[r_, pick] = np.clip(np.round(AOs[r_, pick] + sign * de * d[0, pick]), 0, d[0, pick])
            xs, us = D.d2(AOs, np.repeat(d, 100, axis=0), 0, 0)
            t = contrast(xs, us, mc, mo, D.Bmat)[0]
            vals.append(t)
        vals = np.concatenate(vals)
        rows.append(dict(pi=pi_, delta=round(de, 4), power=float(np.mean(vals >= q999)), median_T=float(np.median(vals))))
PW = pd.DataFrame(rows)
PW.to_csv(f"{c.RUN}/02_A2p_power.tsv", sep="\t", index=False, float_format="%.4f")
print(PW.to_string(index=False))
power_gate = float(PW[(PW.pi == 0.33) & (np.isclose(PW.delta, 1 / 15, atol=1e-3))].power.iloc[0])

# A6 positional (descriptive)
rows = []
for w in range(0, int(D.pos[is4].max() // 100_000_000) + 1):
    mw = is4 & (D.pos // 100_000_000 == w)
    if (mw & used).sum() == 0:
        continue
    t, th, se = [v[0] for v in contrast(x, u, mw, oth, D.Bmat)]
    rows.append(dict(window_Mb=f"{w * 100}-{(w + 1) * 100}", sites=int((mw & used).sum()), diff_mean_D2=th, se=se,
                     lo95=th - 1.96 * se, hi95=th + 1.96 * se))
A6 = pd.DataFrame(rows)
A6.to_csv(f"{c.RUN}/02_A6_positional.tsv", sep="\t", index=False, float_format="%.5f")
print(A6.to_string(index=False))

# A8 (i) depth >= 200
D200 = Data(S["A9_200"])
i4, o4 = D200.chrom == "chr4", D200.chrom != "chr4"
a2, d2_ = clone_pool(D200, c.CLONE)
x2, u2 = D200.d2(a2, d2_, 0, 0)
T200 = float(contrast(x2, u2, i4, o4, D200.Bmat)[0][0])
n200 = pseudo_T(D200, 2000, c.SEED + 1, i4, o4)
p200 = (1 + np.sum(n200 >= T200)) / (1 + len(n200))
robust = bool(p200 < 0.01)
print(f"A8(i) depth>=200: sites {len(D200.chrom)}; T {T200:+.3f}; p {p200:.4f} (2,000 pseudo-clones); robust (p<0.01): {robust}")

# A8 (ii) leave-one-stem-out and without pilot stems (descriptive, vs the primary null)
rows = []
for lab, stems in [(f"without {s}", [t for t in c.CLONE if t != s]) for s in c.CLONE] + \
                  [("without pilot stems", [t for t in c.CLONE if t not in c.PILOT_STEMS])]:
    a_, d_ = clone_pool(D, stems)
    xx, uu = D.d2(a_, d_, 0, 0)
    t = float(contrast(xx, uu, is4, oth, D.Bmat)[0][0])
    rows.append(dict(pool=lab, T=t, p_vs_primary_null=(1 + np.sum(null >= t)) / (1 + len(null))))
LO = pd.DataFrame(rows)
LO.to_csv(f"{c.RUN}/02_A8_leave_out.tsv", sep="\t", index=False, float_format="%.5f")
print(LO.to_string(index=False))

# A8 (iii) single- vs multi-Neale-locus strata (descriptive; own 2,000-pseudo-clone null for context)
rows = []
for lab, ms in [("single Neale locus", D.nloci == 1), ("multiple Neale loci", D.nloci >= 2)]:
    t = float(contrast(x, u, is4 & ms, oth & ms, D.Bmat)[0][0])
    nn = pseudo_T(D, 2000, c.SEED + 2, is4 & ms, oth & ms)
    rows.append(dict(stratum=lab, sites_chr4=int((is4 & ms & used).sum()), sites_other=int((oth & ms & used).sum()),
                     T=t, p=(1 + np.sum(nn >= t)) / (1 + len(nn))))
ST = pd.DataFrame(rows)
ST.to_csv(f"{c.RUN}/02_A8_strata.tsv", sep="\t", index=False, float_format="%.5f")
print(ST.to_string(index=False))

# decision (PREREG_v2 A9d)
if p < 0.001 and T4 > maxTc and robust:
    call, code = "SUPPORTS", ""
elif power_gate >= 0.80 and p > 0.20 and T4 <= q90Tc:
    call, code = "CONTRADICTS", ""
else:
    call = "INCONCLUSIVE"
    if p < 0.001 and not (T4 > maxTc):
        code = "M"
    elif p < 0.001 and not robust:
        code = "R"
    elif power_gate < 0.80:
        code = "U"
    else:
        code = "N"
summ = dict(T_chr4=T4, theta=th4, se=se4, p=p, null_q999=q999, max_T_c=maxTc, q90_T_c=q90Tc,
            power_pi033_delta1_15=power_gate, T200=T200, p200=p200, robust=robust, call=call, reason_code=code,
            sites_total=int(len(D.chrom)), sites_chr4_used=int((used & is4).sum()), sites_used=int(used.sum()))
json.dump(summ, open(f"{c.RUN}/02_A9_summary.json", "w"), indent=1)
print("A9 DECISION:", json.dumps(summ))
