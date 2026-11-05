#!/usr/bin/env python3
# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Batch-aware clonality test: removes PCR-group offsets from all pairwise distances and
# compares the stems with the bound from same-tree pairs (Methods and Results: clonality).
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
"""04: B1 (batch-corrected identity-side rule) and B4 (same on the uncorrected matrix), plus B3 inputs (PREREG_v2).

Correction: OLS over all pairs of distinct trees in the popgen 224 of d ~ piecewise-linear spline in log1p(km)
(knots at the 20/40/60/80% quantiles) + region-pair factor + count indicators for PCR groups 2-5 (group 1 = pilot capture
baseline). d_c = d - beta_gi - beta_gj. Identity bound U = mean + t(0.995, 6) SD sqrt(1 + 1/7) over the 7 identity pairs.
Outputs: RUN/04_B1_pairs_{corrected,uncorrected}.tsv, RUN/04_B1_summary.json; stdout -> RUN/04_B1_clonality.log
"""
import itertools
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dc_common as c  # noqa: E402

inv = c.load_inventory()
P224 = list(inv.index[inv.in_canonical_224 == True])  # noqa: E712
SOUTH44 = [s for s in P224 if inv.loc[s, "lat"] < c.SOUTH and s not in (c.WPR, c.C1)]
print(f"popgen 224: {len(P224)}; southern anchors (lat < {c.SOUTH}, excl. WPR, C1): {len(SOUTH44)}")


def readmat(p):
    m = pd.read_csv(p, index_col=0)
    m.index = [s.replace("RedwoodCap_", "") for s in m.index]
    m.columns = list(m.index)
    return m


def run(label, path):
    M = readmat(path)
    g = {s: c.pcr_group(s) for s in M.index}
    lat, lon, reg = inv.lat, inv.lon, inv.region_lat
    rows = []
    for a, b in itertools.combinations(P224, 2):
        rows.append((a, b, M.loc[a, b]))
    F = pd.DataFrame(rows, columns=["a", "b", "d"])
    F["km"] = c.hav_km(lat[F.a].values, lon[F.a].values, lat[F.b].values, lon[F.b].values)
    F["rp"] = ["-".join(sorted((reg[a], reg[b]))) for a, b in zip(F.a, F.b)]
    F = F.dropna(subset=["km", "d"])
    x = np.log1p(F.km.values)
    knots = np.quantile(x, [0.2, 0.4, 0.6, 0.8])
    X = [np.ones(len(F)), x] + [np.maximum(x - k, 0) for k in knots]
    rps = sorted(F.rp.unique())
    X += [(F.rp.values == r).astype(float) for r in rps[1:]]
    ga, gb = np.array([g[s] for s in F.a]), np.array([g[s] for s in F.b])
    X += [((ga == gg).astype(float) + (gb == gg)) for gg in (2, 3, 4, 5)]
    X = np.column_stack(X)
    coef, *_ = np.linalg.lstsq(X, F.d.values, rcond=None)
    beta = {1: 0.0, **{gg: float(coef[-4 + i]) for i, gg in enumerate((2, 3, 4, 5))}}
    resid_sd = float(np.std(F.d.values - X @ coef))
    print(f"[{label}] OLS over {len(F):,} pairs; region pairs {rps}; batch offsets beta (group 1 = 0): "
          + ", ".join(f"g{k} {v:+.5f}" for k, v in beta.items()) + f"; residual SD {resid_sd:.5f}")

    def dc(a, b):
        return M.loc[a, b] - beta[g[a]] - beta[g[b]]

    out = []
    for a, b in c.IDENTITY:
        out.append(("identity", a, b, M.loc[a, b], dc(a, b)))
    for a, b in itertools.combinations(c.CLONE, 2):
        out.append(("clone", a, b, M.loc[a, b], dc(a, b)))
    for s in c.CLONE:
        out.append(("anchor_WPR", s, c.WPR, M.loc[s, c.WPR], dc(s, c.WPR)))
        out.append(("anchor_C1", s, c.C1, M.loc[s, c.C1], dc(s, c.C1)))
        for t in SOUTH44:
            out.append(("anchor_south", s, t, M.loc[s, t], dc(s, t)))
    O = pd.DataFrame(out, columns=["cls", "a", "b", "d_raw", "d_c"])
    O["ga"], O["gb"] = [g[a] for a in O.a], [g[b] for b in O.b]
    idc = O[O.cls == "identity"]
    U = float(idc.d_c.mean() + stats.t.ppf(0.995, 6) * idc.d_c.std(ddof=1) * np.sqrt(1 + 1 / 7))
    sd_raw, sd_c = float(idc.d_raw.std(ddof=1)), float(idc.d_c.std(ddof=1))
    cl = O[O.cls == "clone"].copy()
    anc = O[O.cls.str.startswith("anchor")]
    min_anchor = float(anc.d_c.min())
    cl["out"] = cl.d_c > U
    one = bool((~cl.out).all() and cl.d_c.max() < min_anchor)
    # B1 part of the CONTRADICTED rule (B2 agreement checked in 06)
    flag = {}
    for s in c.CLONE:
        sub = cl[(cl.a == s) | (cl.b == s)]
        partners = [(b if a == s else a) for a, b in zip(sub.a, sub.b)]
        outp = [p for p, o in zip(partners, sub.out) if o]
        flag[s] = sum(c.pcr_group(p) != c.pcr_group(s) for p in outp) >= 2
    hiblk = {s for s in c.CLONE if c.pcr_group(s) in (1, 3)}
    outset = cl[cl.out]
    batch_aligned = bool(len(outset) > 0 and all(a in hiblk and b in hiblk for a, b in zip(outset.a, outset.b)))
    b1_contra_candidates = [s for s, f in flag.items() if f] if not batch_aligned else []
    call = "ONE GENET" if one else "AMBIGUOUS"
    O.to_csv(f"{c.RUN}/04_B1_pairs_{label}.tsv", sep="\t", index=False, float_format="%.6f")
    blk = {"low (2/4/5)": cl[cl.ga.isin([2, 4, 5]) & cl.gb.isin([2, 4, 5])].d_c,
           "high (1/3)": cl[cl.ga.isin([1, 3]) & cl.gb.isin([1, 3])].d_c}
    blk["mixed"] = cl.d_c.drop(blk["low (2/4/5)"].index).drop(blk["high (1/3)"].index)
    stem = {s: float(cl[(cl.a == s) | (cl.b == s)].d_c.median()) for s in c.CLONE}
    res = dict(label=label, beta=beta, U=U, identity_sd_raw=sd_raw, identity_sd_corrected=sd_c,
               correction_reduces_identity_spread=sd_c < sd_raw, identity_dc=idc.d_c.round(5).tolist(),
               clone_dc_min=float(cl.d_c.min()), clone_dc_max=float(cl.d_c.max()), n_clone_out=int(cl.out.sum()),
               min_anchor_dc=min_anchor, min_anchor_by=dict(anc.groupby("cls").d_c.min().round(5)),
               call_B1_part=call, batch_aligned=batch_aligned, contradicted_candidates_B1=b1_contra_candidates,
               B3_block_medians={k: float(v.median()) for k, v in blk.items()}, B5_stem_median_dc=stem)
    print(json.dumps(res, indent=1, default=float))
    return res


res = {"corrected": run("corrected", c.MAT_C), "uncorrected": run("uncorrected", c.MAT_U)}
json.dump(res, open(f"{c.RUN}/04_B1_summary.json", "w"), indent=1, default=float)
