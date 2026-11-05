#!/usr/bin/env python
# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Flow-cytometry figure: DNA-content histograms of the clone stem and two reference trees,
# each normalized to its own giant sequoia standard peak. The ratios and genome sizes in its legend are typed in
# from a table made by a script that is not on the server (see README).
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
"""
Generate the SEGI-anchored clone-vs-reference flow-cytometry figure for the clone paper.
Run with:  ~/miniconda3/envs/flowcyto/bin/python make_clone_vs_reference_figure.py
Reads 3 tubes from Experiment_002 and writes:
    Experiment_002/clones_vs_reference_SEGIanchored_no_title.png

Normalization (both axes anchored to each tube's own SEGI G1 peak, so tubes run at
different PMT voltages are comparable):
  x-axis = each nucleus's BUV 496-A signal / that tube's SEGI G1 peak  (SEGI G1 = 1.0)
  y-axis = smoothed histogram / its SEGI G1 peak height               (SEGI G1 peak = 1)
Gate: debris > 5000, relaxed singlet band (A/H within 0.6-1.7x median).
Legend ratios/2C are the established (validated) values for each tube.
"""
import os
import fcsparser, numpy as np
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

DIR = "/Users/sashanikolaeva/Desktop/Flow_data/Experiment_002"
A, H = "BUV 496-A", "BUV 496-H"

# (file key, filename, kind, color, linestyle, established ratio, established 2C Gb)
TUBES = [
    ("SEGI+C2_005",    "Specimen_001_SEGI+C2_005.fcs",     "clone",     "#1b7a2f", "-",  2.99, 48.6),
    ("SEGI+NR1_002",   "Specimen_001_SEGI+NR1_002.fcs",    "reference", "#8c1d2b", "--", 3.06, 49.8),
    ("SEGI+NR2_1_015", "Specimen_001_SEGI+NR2_1_015.fcs",  "reference", "#f0a030", "--", 3.06, 49.8),
]
# C7 (Specimen_001_SEGI+C7_010.fcs) was dropped: redwood peak only ~3-4 nuclei, flagged atypical.

def load(fn):
    meta, df = fcsparser.parse(fn, reformat_meta=True)
    a = df[A].values.astype(float); h = df[H].values.astype(float)
    g = np.isfinite(a) & np.isfinite(h) & (a > 0) & (h > 0); a, h = a[g], h[g]
    sel = a > 5000; a2, h2 = a[sel], h[sel]           # debris gate
    r = a2 / h2; m = np.median(r); a2 = a2[(r > m * 0.6) & (r < m * 1.7)]  # relaxed singlet band
    return a2

def segi_peak(a):
    """Lowest prominent nuclei peak = SEGI standard (within-tube rule)."""
    top = min(a.max(), 262143)
    cnt, edg = np.histogram(a, bins=250, range=(5000, top * 1.02)); ctr = 0.5 * (edg[:-1] + edg[1:])
    sc = gaussian_filter1d(cnt.astype(float), 2)
    p, pr = find_peaks(sc, prominence=max(sc.max() * 0.03, 2), distance=4)
    cand = sorted([(ctr[i], pr['prominences'][j]) for j, i in enumerate(p)])
    return cand[0][0] if cand else np.median(a)

fig, ax = plt.subplots(figsize=(13, 6))
ratios = []
for key, fn, kind, col, ls, er, eg in TUBES:
    a = load(os.path.join(DIR, fn))
    s = segi_peak(a)
    x = a / s                                          # x-normalization -> SEGI G1 = 1.0
    cnt, edg = np.histogram(x, bins=160, range=(0.4, 4.3)); ctr = 0.5 * (edg[:-1] + edg[1:])
    sc = gaussian_filter1d(cnt.astype(float), 2.2)     # smoothing -> single clean peaks
    i1 = np.argmin(np.abs(ctr - 1.0)); peak_segi = sc[max(0, i1 - 4):i1 + 5].max()
    y = sc / peak_segi                                 # y-normalization -> SEGI G1 peak = 1
    ax.plot(ctr, y, ls, color=col, lw=1.9,
            label=f"Specimen_001_{key}   ({kind}, {er:.2f}x, {eg:.1f} Gb)")
    ratios.append(er)

ax.axvline(1.0, color="k", ls="--", lw=1); ax.text(1.02, 1.03, "SEGI = 1.0", fontsize=10)
ax.axvspan(min(ratios), max(ratios), color="0.6", alpha=0.18)
ax.annotate("redwood G1\n(≈ 3× SEGI)", xy=(3.0, 0.62), xytext=(3.45, 0.42),
            fontsize=10, color="0.4", arrowprops=dict(arrowstyle="->", color="0.5"))
ax.set_xlim(0.4, 4.3); ax.set_ylim(0, 1.08)
ax.set_xlabel("DNA content relative to SEGI")
ax.set_ylabel("normalized nuclei count")
ax.legend(fontsize=9, loc="upper right", framealpha=0.9)
plt.tight_layout()
out = os.path.join(DIR, "clones_vs_reference_SEGIanchored_no_title.png")
plt.savefig(out, dpi=140)
print("wrote", out)
