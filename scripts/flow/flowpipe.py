# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Flow-cytometry genome-size pipeline with quality grades set in advance. Its grades (all
# three tubes failing) are no longer reported in the paper; kept as a record of the flow analysis.
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
"""Flow-cytometry genome-size pipeline (decision 2.9 / proposal 4b): FCS in; per-tube table, QC tiers, summaries,
figures out. Library used by run_flow.py; parameters come from a per-run PREREG_flow_<run>.tsv (frozen) or the
defaults below.

Steps
 1 ingest      header voltages, $TOT, date, cytometer; DNA channel auto-detect (BUV 496-A > Pacific Blue-A > BV421-A,
               or the sheet's dna_channel); share of events at the ceiling on DNA / FSC-A / SSC-A
 2 gates       per stratum (run date x DNA channel x voltage set) from the standard-only tube(s): debris cut at
               debris_frac x standard G1; singlets = DNA-A/DNA-H within median +- singlet_mad MAD.
               A stratum without a standard-only tube uses a voltage-scaled reference (gain ~ V^gain_exp) and a
               per-tube singlet window, flagged gate_fallback (DEVIATIONS DV08)
 3 peaks       log10 histogram maxima; standard G1 = most prominent maximum in std_win x reference; redwood G1 =
               most prominent maximum in red_win x standard G1; Gaussian + exponential background fitted by least
               squares in two windows (+-fit_w1, +-fit_w2); G2 searched at 2 x G1 +- g2_win
 4 ratio       redwood G1 / in-tube standard G1 (primary = window 1, fallback window 2); suppression = in-tube /
               standard-only G1
 5 QC tiers    PASS-A (proposal criteria), PASS-B (declared, looser), FAIL
 6 summaries   per stratum; per tree x stratum (PASS-B and above; unit = independent chop); clone vs euploid
               delta only within a stratum, only when every tube has a known tree and the run's PREREG is frozen
 7 figures     per-tube histograms with gates and fits; ratio by stratum / tree; suppression vs ratio
"""
import math
import os

import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

import fcsio

DEFAULTS = dict(
    dna_channels="BUV 496-A,Pacific Blue-A,BV421-A",
    debris_frac=0.5, singlet_mad=4.0, gain_exp=7.65,
    std_win_lo=0.70, std_win_hi=1.30, std_below_lo=0.45,
    red_win_lo=2.6, red_win_hi=3.4, g2_win=0.06,
    fit_w1=0.15, fit_w2=0.25, fit_sensitive_pct=1.0,
    A_red_n=1000, A_std_n=500, A_cv=5.0, A_g2g1_lo=1.95, A_g2g1_hi=2.05, A_suppr=0.85, A_ceiling_pct=1.0,
    B_n=100, B_cv=8.0,
    std_2C_pg_assembly=16.62, std_2C_pg_literature="NA",
    boot_reps=2000, seed=20261026,
    out_consistent_lo=-2.3, out_consistent_hi=-0.8, out_null_halfwidth=0.8, out_inconsistent_below=-3.0,
    expected_delta_pct=-1.54,
)


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return v


# ---------------------------------------------------------------- 1 ingest
def ingest(path, dna_override=""):
    names, X, kv = fcsio.read(path)
    pm = fcsio.param_meta(kv)
    prefs = [dna_override] if dna_override else DEFAULTS["dna_channels"].split(",")
    dna = next((c for c in prefs if c in names), None)
    if dna is None:
        raise ValueError(f"{path}: no DNA channel among {prefs}; channels {names}")
    hch = dna[:-2] + "-H" if dna.endswith("-A") and (dna[:-2] + "-H") in names else None
    col = {n: X[:, i] for i, n in enumerate(names)}
    rng = pm[dna]["R"]
    ceil = rng - 144 if rng > 1000 else rng
    def share(c):
        return 100 * float((col[c] >= ceil).mean()) if c in col else float("nan")
    vs = f"FSC{pm.get('FSC-A', {}).get('V', '')}/SSC{pm.get('SSC-A', {}).get('V', '')}/DNA{pm[dna]['V']}"
    meta = dict(events=int(kv["$TOT"].strip()), date=kv.get("$DATE", "").strip(), btim=kv.get("$BTIM", "").strip(),
                cytometer=kv.get("$CYT", "").strip(), tube_name=kv.get("TUBE NAME", "").strip(), dna_channel=dna,
                dna_h_channel=hch or "", dna_V=pm[dna]["V"], fsc_V=pm.get("FSC-A", {}).get("V", ""),
                ssc_V=pm.get("SSC-A", {}).get("V", ""), voltage_set=vs, range=rng,
                ceil_dna_pct=share(dna), ceil_fsc_pct=share("FSC-A"), ceil_ssc_pct=share("SSC-A"),
                all_voltages=";".join(f"{n}={pm[n]['V']}" for n in names if pm[n]["V"]))
    return meta, col[dna], (col[hch] if hch else None), ceil


# ---------------------------------------------------------------- 3 peaks
def maxima(x):
    if len(x) < 50:
        return []
    lx = np.log10(x)
    bins = np.arange(lx.min(), lx.max() + 0.003, 0.003)
    if len(bins) < 5:
        return []
    c, e = np.histogram(lx, bins)
    sm = gaussian_filter1d(c.astype(float), 2)
    ctr = 10 ** (0.5 * (e[1:] + e[:-1]))
    pk, pr = find_peaks(sm, prominence=max(3, 0.02 * sm.max()))
    return [(float(ctr[p]), float(pr["prominences"][k])) for k, p in enumerate(pk)]


def pick(cands, lo, hi):
    c = [t for t in cands if lo <= t[0] <= hi]
    return (max(c, key=lambda t: t[1])[0], len(c)) if c else (None, 0)


def _model(x, A, mu, s, B, k, lo, span):
    return A * np.exp(-0.5 * ((x - mu) / s) ** 2) + B * np.exp(-k * (x - lo) / span)


def fit_peak(x, c0, W):
    """Gaussian + exponential background on the linear histogram within c0 * (1 +- W)."""
    out = dict(mean=np.nan, cv=np.nan, n=0, chi2=np.nan, converged=False)
    if c0 is None:
        return out
    lo, hi = c0 * (1 - W), c0 * (1 + W)
    seg = x[(x > lo) & (x < hi)]
    if len(seg) < 30:
        return out
    c, e = np.histogram(seg, 60, (lo, hi))
    ctr = 0.5 * (e[1:] + e[:-1])
    f = lambda xx, A, mu, s, B, k: _model(xx, A, mu, s, B, k, lo, hi - lo)
    try:
        p, _ = curve_fit(f, ctr, c, p0=[max(c.max(), 1), c0, 0.04 * c0, max(np.median(c[:5]), 0.1), 1.0],
                         bounds=([0, lo, 1e-3 * c0, 0, -20], [np.inf, hi, W * c0, np.inf, 20]), maxfev=20000)
    except Exception:
        return out
    mu, s = p[1], abs(p[2])
    resid = c - f(ctr, *p)
    chi2 = float(np.sum(resid ** 2 / np.maximum(f(ctr, *p), 1)) / max(len(c) - 5, 1))
    conv = bool(lo < mu < hi and s < 0.5 * W * c0 and p[0] > 0)
    out.update(mean=float(mu), cv=float(100 * s / mu), n=int(((x > mu - 2 * s) & (x < mu + 2 * s)).sum()),
               chi2=chi2, converged=conv)
    return out


def std_reference(a, h, ceil, P):
    """Standard-only tube: G1 and the A/H singlet window. Most prominent maximum above 1% of range."""
    m0 = (a > 0.01 * ceil) & (a < ceil)
    mx = maxima(a[m0])
    if not mx:
        return None
    c0 = max(mx, key=lambda t: t[1])[0]
    f = fit_peak(a[m0], c0, P["fit_w1"])
    if not f["converged"]:
        return None
    sel = (a > P["debris_frac"] * f["mean"]) & (a < ceil)
    win = None
    if h is not None:
        ok = sel & (h > 0)
        r = a[ok] / h[ok]
        med = float(np.median(r)); mad = float(1.4826 * np.median(np.abs(r - med)))
        win = (med - P["singlet_mad"] * mad, med + P["singlet_mad"] * mad)
    return dict(g1=f["mean"], cv=f["cv"], n=f["n"], ah_window=win)


def analyse_tube(a, h, ceil, ref_g1, ah_window, P):
    """Steps 2-4 for one mixed tube. Returns a dict of per-tube results and flags."""
    r = dict(gated_events=0, singlets=0, std_g1_w1=np.nan, std_cv_w1=np.nan, std_n_w1=0, std_g1_w2=np.nan,
             std_cv_w2=np.nan, std_n_w2=0, red_g1_w1=np.nan, red_cv_w1=np.nan, red_n_w1=0, red_g1_w2=np.nan,
             red_cv_w2=np.nan, red_n_w2=0, std_g2g1=np.nan, red_g2g1=np.nan, ratio_w1=np.nan, ratio_w2=np.nan,
             ratio=np.nan, ratio_window_used="", ratio_diff_pct=np.nan, suppression=np.nan, ceiling_gated_pct=np.nan,
             converged=False, flags=[], x=None, fits={})
    m = (a > P["debris_frac"] * ref_g1) & (a < ceil)
    atc = (a >= ceil).sum()
    if h is not None:
        m &= h > 0
        rr = np.where(m, a / np.where(h > 0, h, 1), np.nan)
        if ah_window is None:   # per-tube window (fallback)
            v = rr[m]
            med = np.median(v); mad = 1.4826 * np.median(np.abs(v - med))
            ah_window = (med - P["singlet_mad"] * mad, med + P["singlet_mad"] * mad)
        s = m & (rr >= ah_window[0]) & (rr <= ah_window[1])
    else:
        s = m
        r["flags"].append("no_height_channel")
    x = a[s]
    r.update(gated_events=int(m.sum()), singlets=int(s.sum()), x=x)
    r["ceiling_gated_pct"] = 100 * atc / max(m.sum() + atc, 1)
    mx = maxima(x)
    s0, ns = pick(mx, P["std_win_lo"] * ref_g1, P["std_win_hi"] * ref_g1)
    below, _ = pick(mx, P["std_below_lo"] * ref_g1, P["std_win_lo"] * ref_g1)
    if below is not None and (s0 is None or max(t[1] for t in mx if t[0] == below) >
                              max(t[1] for t in mx if t[0] == s0)):
        r["flags"].append("std_peak_below_window")
    if ns > 1:
        r["flags"].append("multi_max_std_window")
    if s0 is None:
        r["flags"].append("no_std_peak")
        return r
    r0, nr = pick(mx, P["red_win_lo"] * s0, P["red_win_hi"] * s0)
    if nr > 1:
        r["flags"].append("multi_max_red_window")
    for tag, W in (("w1", P["fit_w1"]), ("w2", P["fit_w2"])):
        fs = fit_peak(x, s0, W)
        fr = fit_peak(x, r0, W) if r0 else fit_peak(x, None, W)
        r.update({f"std_g1_{tag}": fs["mean"], f"std_cv_{tag}": fs["cv"], f"std_n_{tag}": fs["n"] if fs["converged"] else 0,
                  f"red_g1_{tag}": fr["mean"], f"red_cv_{tag}": fr["cv"], f"red_n_{tag}": fr["n"] if fr["converged"] else 0})
        r["fits"][tag] = (fs, fr)
        if fs["converged"] and fr["converged"]:
            r[f"ratio_{tag}"] = fr["mean"] / fs["mean"]
    if r0 is None:
        r["flags"].append("no_redwood_peak")
    if not math.isnan(r["ratio_w1"]):
        r["ratio"], r["ratio_window_used"] = r["ratio_w1"], "w1"
        r["converged"] = True
    elif not math.isnan(r["ratio_w2"]):
        r["ratio"], r["ratio_window_used"] = r["ratio_w2"], "w2"
        r["converged"] = True
        r["flags"].append("w1_fit_failed")
    if not (math.isnan(r["ratio_w1"]) or math.isnan(r["ratio_w2"])):
        r["ratio_diff_pct"] = 100 * abs(r["ratio_w1"] / r["ratio_w2"] - 1)
        if r["ratio_diff_pct"] > P["fit_sensitive_pct"]:
            r["flags"].append("fit_sensitive")
    # G2 / G1
    tw = r["ratio_window_used"] or "w1"
    sg = r[f"std_g1_{tw}"]
    if not math.isnan(sg):
        g2, _ = pick(mx, 2 * sg * (1 - P["g2_win"]), 2 * sg * (1 + P["g2_win"]))
        f2 = fit_peak(x, g2, 0.10) if g2 else None
        if f2 and f2["converged"]:
            r["std_g2g1"] = f2["mean"] / sg
    rg = r[f"red_g1_{tw}"]
    if not math.isnan(rg):
        g2, _ = pick(mx, 2 * rg * (1 - P["g2_win"]), 2 * rg * (1 + P["g2_win"]))
        f2 = fit_peak(x, g2, 0.10) if g2 else None
        if f2 and f2["converged"]:
            r["red_g2g1"] = f2["mean"] / rg
    return r


# ---------------------------------------------------------------- 5 tiers
def tier(r, P, suppression_known=True):
    w = r["ratio_window_used"] or "w1"
    sn, rn = r[f"std_n_{w}"], r[f"red_n_{w}"]
    scv, rcv = r[f"std_cv_{w}"], r[f"red_cv_{w}"]
    if not r["converged"]:
        return "FAIL", "no converged ratio"
    reasons = []
    a_ok = (rn >= P["A_red_n"] and sn >= P["A_std_n"] and scv <= P["A_cv"] and rcv <= P["A_cv"]
            and P["A_g2g1_lo"] <= (r["std_g2g1"] if not math.isnan(r["std_g2g1"]) else -1) <= P["A_g2g1_hi"]
            and suppression_known and not math.isnan(r["suppression"]) and r["suppression"] >= P["A_suppr"]
            and r["ceiling_gated_pct"] < P["A_ceiling_pct"])
    if a_ok:
        return "PASS-A", ""
    if rn < P["A_red_n"]: reasons.append(f"red n {rn}<{P['A_red_n']}")
    if sn < P["A_std_n"]: reasons.append(f"std n {sn}<{P['A_std_n']}")
    if scv > P["A_cv"] or rcv > P["A_cv"]: reasons.append(f"CV {scv:.1f}/{rcv:.1f}>{P['A_cv']}")
    if math.isnan(r["std_g2g1"]) or not P["A_g2g1_lo"] <= r["std_g2g1"] <= P["A_g2g1_hi"]: reasons.append("std G2/G1")
    if not suppression_known or math.isnan(r["suppression"]) or r["suppression"] < P["A_suppr"]: reasons.append("suppression")
    b_ok = (rn >= P["B_n"] and sn >= P["B_n"] and scv <= P["B_cv"] and rcv <= P["B_cv"]
            and "fit_sensitive" not in r["flags"] and r["ratio_window_used"] == "w1")
    if b_ok:
        return "PASS-B", "not A: " + "; ".join(reasons)
    rb = []
    if rn < P["B_n"] or sn < P["B_n"]: rb.append(f"n {sn}/{rn}<{P['B_n']}")
    if scv > P["B_cv"] or rcv > P["B_cv"]: rb.append(f"CV {scv:.1f}/{rcv:.1f}>{P['B_cv']}")
    if "fit_sensitive" in r["flags"]: rb.append("fit-sensitive")
    if r["ratio_window_used"] != "w1": rb.append("w1 fit failed")
    return "FAIL", "; ".join(rb)


# ---------------------------------------------------------------- 6 contrast (designed runs only)
def contrast(rows, P, rng):
    """Clone vs euploid delta % within one stratum. rows: dicts with tree, tree_class, chop_id, ratio.
    Unit = independent chop (tubes of one chop averaged); two-level bootstrap (trees, then chops)."""
    chops = {}
    for r in rows:
        chops.setdefault((r["tree_class"], r["tree"], r["chop_id"]), []).append(r["ratio"])
    tree = {}
    for (cl, t, c), v in chops.items():
        tree.setdefault((cl, t), []).append(float(np.mean(v)))
    clone = [v for (cl, t), v in tree.items() if cl == "clone"]
    eup = [v for (cl, t), v in tree.items() if cl == "euploid"]
    if len(clone) != 1 or len(eup) < 2:
        return None
    c = np.array(clone[0])
    E = [np.array(v) for v in eup]
    est = lambda cc, ee: 100 * (cc.mean() - np.mean([e.mean() for e in ee])) / np.mean([e.mean() for e in ee])
    d = est(c, E)
    bs = []
    for _ in range(int(P["boot_reps"])):
        cc = rng.choice(c, len(c))
        ee = [E[i] for i in rng.integers(0, len(E), len(E))]
        ee = [rng.choice(e, len(e)) for e in ee]
        bs.append(est(cc, ee))
    lo, hi = np.percentile(bs, [2.5, 97.5])
    if lo > 0 or hi < 0:
        if P["out_consistent_lo"] <= d <= P["out_consistent_hi"] and hi < 0:
            v = "consistent with one-copy loss"
        elif d < P["out_inconsistent_below"] and not lo <= P["expected_delta_pct"] <= hi:
            v = "inconsistent with a single chr4 copy (re-examine)"
        else:
            v = "inconclusive"
    else:
        v = "does not support" if (hi - lo) / 2 <= P["out_null_halfwidth"] else "inconclusive"
    return dict(delta_pct=d, ci_lo=lo, ci_hi=hi, n_clone_chops=len(c), n_euploid_trees=len(E), verdict=v)


def power_table():
    """Genet-level power (review flow_04): one clone vs k euploid trees, one-sided alpha 0.05, 80% power."""
    from scipy import stats, optimize
    rows = []
    for k in (6, 8, 10, None):
        if k is None:
            need = stats.norm.isf(0.05) + stats.norm.isf(0.20)
        else:
            df = k - 1
            tc = stats.t.isf(0.05, df)
            f = lambda nc: stats.nct.sf(tc, df, nc) - 0.80
            nc = optimize.brentq(f, 0.1, 20)
            need = nc * math.sqrt(1 + 1 / k)   # delta / sigma_T for prediction-interval t
        rows.append((k if k else "inf", need, 1.54 / need))
    return rows
