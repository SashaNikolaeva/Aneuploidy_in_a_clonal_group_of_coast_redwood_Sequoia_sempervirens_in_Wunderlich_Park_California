# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Check of how the singlet gate changes the peak positions (strict versus relaxed gate;
# Results: the clone peak lies between 2.99x and 3.13x depending on how the peaks were found).
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
"""Relax the singlet gate on sparse tubes to recover fittable peaks.

STRICT  = singlet band 0.7-1.4x median A/H  (current protocol)
RELAXED = singlet band 0.6-1.7x            (keeps more events; we showed doublets
          shift size only ~1%, so this is safe)
Report, per tube, events kept and whether a 2-peak (SEGI + redwood) fit + size
is obtained each way. Flag tubes that become NEWLY fittable under the relaxed gate.
SEGI-only controls legitimately have ONE peak -> not a failure.
"""
import fcsparser, numpy as np, glob, os, csv
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import curve_fit

DIR = "/Users/sashanikolaeva/Desktop/Flow_data/Experiment_002"
A, H = "BUV 496-A", "BUV 496-H"
SEGI_2C_GBP, CVMAX, DEBRIS = 16.25, 22.0, 5000
SPARSE_N = 800   # tubes with fewer kept nuclei than this are "sparse"

def label(fn): return os.path.basename(fn).replace("Specimen_001_", "").replace(".fcs", "")
def gauss(x, a, mu, sig): return a*np.exp(-0.5*((x-mu)/sig)**2)

def load(fn, lo, hi):
    meta, df = fcsparser.parse(fn, reformat_meta=True)
    a = df[A].values.astype(float); h = df[H].values.astype(float)
    g = np.isfinite(a) & np.isfinite(h) & (a > 0) & (h > 0); a, h = a[g], h[g]
    a, h = a[a > DEBRIS], h[a > DEBRIS]
    if len(a) > 50:
        r = a/h; m = np.median(r); a = a[(r > lo*m) & (r < hi*m)]
    return a

def fit_peak(a, center, half=9000):
    seg = a[(a >= center-half) & (a <= center+half)]
    if len(seg) < 20: return None
    cnt, edg = np.histogram(seg, bins=40, range=(center-half, center+half)); ctr = 0.5*(edg[:-1]+edg[1:])
    try:
        popt, _ = curve_fit(gauss, ctr, cnt, p0=[cnt.max(), center, half/3], maxfev=5000)
        mu, sig = popt[1], abs(popt[2])
        if not (center-half < mu < center+half): mu, sig = np.median(seg), seg.std()
    except Exception:
        mu, sig = np.median(seg), seg.std()
    return mu, sig/mu*100.0

def two_peaks(a):
    if len(a) < 60: return []
    top = min(a.max(), 262143)
    cnt, edg = np.histogram(a, bins=400, range=(10000, top*1.02)); ctr = 0.5*(edg[:-1]+edg[1:])
    sc = gaussian_filter1d(cnt.astype(float), 2)
    pk, props = find_peaks(sc, prominence=max(sc.max()*0.06, 3), distance=6)
    cand = [(ctr[p], props['prominences'][k]) for k, p in enumerate(pk) if ctr[p] < 255000]
    fits = []
    for pos, _ in sorted(cand, key=lambda t: -t[1])[:4]:
        f = fit_peak(a, pos)
        if f and f[1] <= CVMAX: fits.append(f)
    return sorted(fits, key=lambda t: t[0])

def result(a):
    pk = two_peaks(a)
    if len(pk) < 2: return len(a), None, None, (pk[0][1] if pk else None)
    segi, samp = pk[0], pk[-1]
    ratio = samp[0]/segi[0]
    return len(a), round(ratio, 3), round(ratio*SEGI_2C_GBP, 2), round(max(segi[1], samp[1]), 1)

rows = []
for fn in sorted(glob.glob(os.path.join(DIR, "*.fcs"))):
    name = label(fn)
    nS, rS, gbS, cvS = result(load(fn, 0.7, 1.4))   # strict
    nR, rR, gbR, cvR = result(load(fn, 0.6, 1.7))   # relaxed
    sparse = nS < SPARSE_N
    newly = (rS is None) and (rR is not None)
    rows.append(dict(sample=name, sparse="yes" if sparse else "", n_strict=nS, n_relax=nR,
                     ratio_strict=rS if rS else "", ratio_relax=rR if rR else "",
                     gb_strict=gbS if gbS else "", gb_relax=gbR if gbR else "",
                     maxCV_relax=cvR if cvR else "", newly_fittable="YES" if newly else ""))

with open(os.path.join(DIR, "relax_gate_sparse.csv"), "w", newline="") as f:
    cols = ["sample","sparse","n_strict","n_relax","ratio_strict","ratio_relax",
            "gb_strict","gb_relax","maxCV_relax","newly_fittable"]
    w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(rows)

print(f"{'sample':24s}{'sparse':>7s}{'nStrict':>8s}{'nRelax':>7s}{'rStrict':>8s}{'rRelax':>8s}{'gbRelax':>8s}{'newly':>7s}")
print("-"*77)
newcount = 0
for r in rows:
    if r['newly_fittable']: newcount += 1
    print(f"{r['sample']:24s}{r['sparse']:>7s}{str(r['n_strict']):>8s}{str(r['n_relax']):>7s}"
          f"{str(r['ratio_strict']):>8s}{str(r['ratio_relax']):>8s}{str(r['gb_relax']):>8s}{r['newly_fittable']:>7s}")
print("-"*77)
print(f"Newly fittable under relaxed gate: {newcount} tube(s).  (SEGI-only controls legitimately have 1 peak.)")
print("saved relax_gate_sparse.csv to", DIR)
