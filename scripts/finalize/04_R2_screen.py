# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Copy-number screen of the 224 survey trees. Not reported in the current version of the
# paper (the other aneuploid survey trees are cited from Nikolaeva et al. 2025, G3); 09_R0_fig6_table1.py reads
# its output for the survey range drawn in Figure 3.
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
"""R2 population screen (SEGI counts), R2-WPR, V5 reconciliation, G3 mapping; plus the R1a positive-control rows
(decision 2.5, 2.2).

Screen: every non-clone sample with > 1 M SEGI reads x chr1-11; z vs PANEL_ALL (sample removed from the panel when a
member); threshold |z| >= z_B, Bonferroni two-sided 0.05 for (n samples x 11) tests.
Classes: whole-copy-like if |c_hat - 5| <= 0.25 or |c_hat - 7| <= 0.25; otherwise fractional.
G3 Table 1 events (Nikolaeva et al. 2025 G3 jkaf063; decision web lookup): gains chr10 x2, chr1, chr3, chr6, chr7, chr8,
chr9; loss chr4 (a southern tree = 33_c6, a clone ramet). Matched to hits by chromosome and direction.
Fig 6 rule (d04): alpha = share(c)/median_PANEL_ALL share(c), chr1-11 denominator, outside [0.95, 1.05].
Outputs: R2_screen_cells.tsv (all cells), R2_hits.tsv, R2_G3_mapping.tsv, R2_V5_reconciliation.tsv, R2_WPR.tsv,
R1a_posctrl.tsv; stdout -> 04_R2_screen.log
"""
import csv
import os
import sys

import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *  # noqa

S = sample_sets()
rc, inv, P = S["rc"], S["inv"], S["PANEL_ALL"]
Pset = set(P)
screen = S["screen"]
n = len(screen)
ntest = n * 11
zB = stats.norm.isf(0.05 / (2 * ntest))
print(f"screen: {n} non-clone samples with > 1 M reads x 11 = {ntest} tests; Bonferroni two-sided 0.05 |z| >= {zB:.3f}")

R = {c: {s: segi_ratio(rc, s, c) for s in S["big"]} for c in CH}
pan = {c: np.array([R[c][x] for x in P]) for c in CH}
pidx = {x: i for i, x in enumerate(P)}


def zc(s, c):
    v = pan[c]
    if s in Pset:
        v = np.delete(v, pidx[s])
    return (R[c][s] - v.mean()) / v.std(ddof=1), 6 * R[c][s] / np.median(v)


share = {s: {c: rc[s][c] / sum(rc[s].values()) for c in CH} for s in S["big"]}
smed = {c: np.median([share[x][c] for x in P]) for c in CH}
alpha = lambda s, c: share[s][c] / smed[c]

v5 = {}
with open(INPUTS["v5_sheet"], encoding="utf-8-sig") as fh:
    for x in csv.DictReader(fh):
        k = [c for c in x if "name" in c.lower()][0]
        v5[short(x[k]).strip()] = x.get("Aneuploid", "").strip().upper()

cells, hits = [], []
for s in S["plot"]:
    for c in CH:
        z, ch = zc(s, c)
        cls = "whole-copy-like" if (abs(ch - 5) <= 0.25 or abs(ch - 7) <= 0.25) else "fractional"
        hit = abs(z) >= zB and s not in CLONE
        cells.append((s, c, R[c][s], z, ch, alpha(s, c), s in CLONE, s in Pset, hit, cls if abs(z) >= zB else ""))
        if hit:
            hits.append((s, c, z, ch, cls))
with open(f"{OUT}/R2_screen_cells.tsv", "w") as f:
    f.write("sample\tchrom\tR\tz_vs_PANEL_ALL\tc_hat\talpha_fig6\tis_clone\tin_PANEL_ALL\tscreen_hit\tclass\n")
    for r in cells:
        f.write(f"{r[0]}\t{r[1]}\t{r[2]:.5f}\t{r[3]:+.3f}\t{r[4]:.3f}\t{r[5]:.4f}\t{r[6]}\t{r[7]}\t{r[8]}\t{r[9]}\n")
with open(f"{OUT}/R2_hits.tsv", "w") as f:
    f.write("sample\tchrom\tz\tc_hat\tclass\tdirection\tin_PANEL_ALL\tV5_flag\tregion_lat\tpcr_pool\tplate_well\tin_wild_mcmc_2024\n")
    print("\n## (i) screen hits (non-clone)")
    for s, c, z, ch, cls in sorted(hits, key=lambda h: (h[0], CH.index(h[1]))):
        d = "gain" if z > 0 else "loss"
        f.write(f"{s}\t{c}\t{z:+.2f}\t{ch:.3f}\t{cls}\t{d}\t{s in Pset}\t{v5.get(s,'-')}\t{inv[s]['region_lat']}\t"
                f"{inv[s]['pcr_pool']}\t{inv[s]['plate_well']}\t{inv[s]['in_wild_mcmc_input_2024']}\n")
        print(f"{s}\t{c}\tz {z:+.2f}\tc_hat {ch:.2f}\t{cls}\t{'panel' if s in Pset else 'excluded sample'}\tV5 {v5.get(s,'-')}")
hs = sorted({h[0] for h in hits})
wc = sorted({h[0] for h in hits if h[4] == "whole-copy-like"})
print(f"hits: {len(hits)} cells in {len(hs)} samples; samples with a whole-copy-like hit: {len(wc)} {wc}; "
      f"samples with only fractional hits: {sorted(set(hs) - set(wc))}")
print(f"(scope: {n} non-clone samples; {sum(s in Pset for s in hs)} of the hit samples are PANEL_ALL genets)")

# clone cells for reference
cz = [zc(s, "chr4") for s in CLONE]
print(f"clone ramets chr4 (not part of the screen): z {min(z for z, _ in cz):+.2f} to {max(z for z, _ in cz):+.2f}; "
      f"c_hat {min(c for _, c in cz):.3f}-{max(c for _, c in cz):.3f}; all |z| >= z_B: {all(abs(z) >= zB for z, _ in cz)}")
other_clone = [(s, c, *zc(s, c)) for s in CLONE for c in CH if c != "chr4" and abs(zc(s, c)[0]) >= zB]
print(f"clone ramets, other chromosomes beyond z_B: {other_clone if other_clone else 'none'}")

# (ii) G3 mapping
G3 = [("chr10", "gain"), ("chr10", "gain"), ("chr1", "gain"), ("chr3", "gain"), ("chr6", "gain"), ("chr7", "gain"),
      ("chr8", "gain"), ("chr9", "gain"), ("chr4", "loss")]
cand = [(s, c, z, ch, cls) for s, c, z, ch, cls in hits]
cand += [("CLONE(9 ramets; G3 ran 33_c6)", "chr4", st.mean(z for z, _ in cz), st.mean(c for _, c in cz), "whole-copy-like")]
used = set()
with open(f"{OUT}/R2_G3_mapping.tsv", "w") as f:
    f.write("G3_event\tG3_chrom\tG3_direction\tmatched_screen_hit\tz\tc_hat\tclass\tin_wild_mcmc_2024\tnote\n")
    print("\n## (ii) G3 Table 1 events vs screen hits (match by chromosome and direction)")
    for i, (c, d) in enumerate(G3, 1):
        m = [h for h in cand if h[1] == c and (h[2] > 0) == (d == "gain") and (h[0], h[1]) not in used]
        m.sort(key=lambda h: -abs(h[2]))
        if m:
            h = m[0]; used.add((h[0], h[1]))
            mc = inv[h[0]]["in_wild_mcmc_input_2024"] if h[0] in inv else "True (33_c6)"
            note = "other candidates: " + ", ".join(f"{x[0]} {x[3]:.2f}" for x in m[1:]) if len(m) > 1 else ""
            f.write(f"G3_{i}\t{c}\t{d}\t{h[0]}\t{h[2]:+.2f}\t{h[3]:.3f}\t{h[4]}\t{mc}\t{note}\n")
            print(f"G3 {c} {d}: {h[0]} (z {h[2]:+.2f}, c_hat {h[3]:.2f}, {h[4]}) {note}")
        else:
            f.write(f"G3_{i}\t{c}\t{d}\tUNMATCHED\t\t\t\t\t\n")
            print(f"G3 {c} {d}: UNMATCHED")
    un = [h for h in cand if (h[0], h[1]) not in used]
    for h in un:
        mc = inv[h[0]]["in_wild_mcmc_input_2024"] if h[0] in inv else ""
        f.write(f"screen_only\t{h[1]}\t{'gain' if h[2] > 0 else 'loss'}\t{h[0]}\t{h[2]:+.2f}\t{h[3]:.3f}\t{h[4]}\t{mc}\t"
                f"screen hit with no G3 event\n")
    print("screen hits with no G3 event:", ", ".join(f"{h[0]} {h[1]} {h[3]:.2f} ({h[4]})" for h in un) or "none")

# (iii) non-clones with chr4 c_hat < 5.75
low = sorted((zc(s, "chr4")[1], s) for s in screen if zc(s, "chr4")[1] < 5.75)
print("\n## (iii) non-clone samples with chr4 c_hat < 5.75:", ", ".join(f"{s} {c:.3f} (z {zc(s,'chr4')[0]:+.2f})" for c, s in low) or "none")
whole4 = [s for c, s in low if abs(c - 5) <= 0.25]
print(f"non-clone whole-copy-like chr4 loss (|c_hat - 5| <= 0.25 and a hit): {whole4 if whole4 else 'none'}")

# (iv) WPR and C1; R2-WPR
with open(f"{OUT}/R2_WPR.tsv", "w") as f:
    f.write("sample\tchrom\tz\tc_hat\talpha_fig6\tbeyond_zB\n")
    for s in [WPR, C1]:
        out = []
        for c in CH:
            z, ch = zc(s, c)
            f.write(f"{s}\t{c}\t{z:+.2f}\t{ch:.3f}\t{alpha(s,c):.4f}\t{abs(z) >= zB}\n")
            out.append(f"{c} z {z:+.1f} a {alpha(s,c):.3f}")
        mz = max(abs(zc(s, c)[0]) for c in CH)
        print(f"\n## (iv) {s}: max |z| {mz:.2f} (z_B {zB:.2f}) -> "
              f"{'deviates' if mz >= zB else 'euploid by depth'}\n  " + "; ".join(out))

# (v) V5 vs screen vs Fig 6 rule
with open(f"{OUT}/R2_V5_reconciliation.tsv", "w") as f:
    f.write("sample\tis_clone\tV5_flag\tscreen_hit_chroms\tfig6_rule_outside\tagreement\n")
    tab = {}
    for s in S["plot"]:
        fl = v5.get(s, "-")
        hc = [c for c in CH if abs(zc(s, c)[0]) >= zB] if s not in CLONE else ["chr4"]
        out = [f"{c}:{alpha(s,c):.3f}" for c in CH if not 0.95 <= alpha(s, c) <= 1.05]
        if fl != "YES" and not hc and not out:
            continue
        key = (fl == "YES", bool(hc), bool(out))
        tab[key] = tab.get(key, 0) + 1
        f.write(f"{s}\t{s in CLONE}\t{fl}\t{','.join(hc)}\t{' '.join(out)}\t"
                f"V5={'Y' if key[0] else 'N'} screen={'Y' if key[1] else 'N'} fig6={'Y' if key[2] else 'N'}\n")
print("\n## (v) V5 flag x screen hit x Fig 6 +-5% rule (counts of samples; clone ramets counted as chr4 hits)")
for k in sorted(tab, reverse=True):
    print(f"  V5 {'YES' if k[0] else 'no '} | screen {'hit' if k[1] else '-  '} | fig6 {'out' if k[2] else 'in '}: {tab[k]}")

# R1a positive-control rows
with open(f"{OUT}/R1a_posctrl.tsv", "w") as f:
    f.write("sample\tchrom\trole\tz\tc_hat\n")
    seen = set()
    for t, c in GAINCTRL.items():
        z, ch = zc(t, c); f.write(f"{t}\t{c}\tGAINCTRL\t{z:+.2f}\t{ch:.3f}\n"); seen.add((t, c))
    for s, c, z, ch, cls in hits:
        if z > 0 and (s, c) not in seen:
            f.write(f"{s}\t{c}\tR2_gain_hit({cls})\t{z:+.2f}\t{ch:.3f}\n")
print("\nR1a positive controls:", open(f"{OUT}/R1a_posctrl.tsv").read().replace("\n", " | "))
