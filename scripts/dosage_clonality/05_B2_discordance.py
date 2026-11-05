#!/usr/bin/env python3
# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Allele-presence discordance between libraries after thinning to 30 reads, a clonality
# measure that does not use the distance matrix (Methods and Results: clonality).
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
"""05: B2 matrix-free allele-presence discordance with depth thinning (PREREG_v2 B2).

Sites: panel depth ratio (median over the 224 of DP / tree median DP) in [0.7, 2.5] and RO+AO >= 30 in both libraries.
Each library thinned to 30 reads (hypergeometric; seeds 20261027 + k, k = 0..19). Eligible: at least one library has
min(ref, alt) >= 3 after thinning. Discordant: one library has >= 3 reads (fraction >= 0.10) of an allele with 0 reads in
the other. Rate = discordant / eligible, averaged over the 20 seeds.
Outputs: RUN/05_B2_pairs.tsv, RUN/05_B2_summary.json; stdout -> RUN/05_B2_discordance.log
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

arr = c.load_array()
inv = c.load_inventory()
ix, P224, R = c.sets(arr["samples"], inv)
RO, AO, DP = arr["RO"], arr["AO"], arr["DP"]
pr = c.panel_ratio(DP, [ix[s] for s in P224])
strat = (pr >= 0.7) & (pr <= 2.5)
RO_s, AO_s = RO[strat].astype(np.int64), AO[strat].astype(np.int64)
print(f"sites in the panel-ratio stratum [0.7, 2.5]: {strat.sum():,} of {len(pr):,}")
SOUTH44 = [s for s in P224 if inv.loc[s, "lat"] < c.SOUTH and s not in (c.WPR, c.C1)]

pairs = [("identity", a, b) for a, b in c.IDENTITY]
pairs += [("clone", a, b) for a, b in itertools.combinations(c.CLONE, 2)]
pairs += [("anchor_WPR", s, c.WPR) for s in c.CLONE] + [("anchor_C1", s, c.C1) for s in c.CLONE]
pairs += [("anchor_south", s, t) for s in c.CLONE for t in SOUTH44]


def rate(a, b):
    ia, ib = ix[a], ix[b]
    ra, aa, rb, ab = RO_s[:, ia], AO_s[:, ia], RO_s[:, ib], AO_s[:, ib]
    m = (ra + aa >= 30) & (rb + ab >= 30)
    ra, aa, rb, ab = ra[m], aa[m], rb[m], ab[m]
    rates, elig = [], []
    for k in range(20):
        rng = np.random.default_rng(c.SEED + k)
        xa = rng.hypergeometric(aa, ra, 30)
        xb = rng.hypergeometric(ab, rb, 30)
        ya, yb = 30 - xa, 30 - xb
        el = (np.minimum(xa, ya) >= 3) | (np.minimum(xb, yb) >= 3)
        dis = ((xa >= 3) & (xb == 0)) | ((xb >= 3) & (xa == 0)) | ((ya >= 3) & (yb == 0)) | ((yb >= 3) & (ya == 0))
        dis &= el
        rates.append(dis.sum() / max(el.sum(), 1))
        elig.append(el.sum())
    return float(np.mean(rates)), int(np.mean(elig)), int(m.sum())


rows = []
for cls, a, b in pairs:
    r_, e_, n_ = rate(a, b)
    rows.append(dict(cls=cls, a=a, b=b, rate=r_, eligible_mean=e_, sites_dp30_both=n_,
                     ga=c.pcr_group(a), gb=c.pcr_group(b)))
O = pd.DataFrame(rows)
O.to_csv(f"{c.RUN}/05_B2_pairs.tsv", sep="\t", index=False, float_format="%.6f")
idr = O[O.cls == "identity"]
U = float(idr.rate.mean() + stats.t.ppf(0.995, 6) * idr.rate.std(ddof=1) * np.sqrt(1 + 1 / 7))
cl = O[O.cls == "clone"].copy()
anc = O[O.cls.str.startswith("anchor")]
cl["out"] = cl.rate > U
one = bool((~cl.out).all() and cl.rate.max() < anc.rate.min())
flag = {}
for s in c.CLONE:
    sub = cl[(cl.a == s) | (cl.b == s)]
    partners = [(b if a == s else a) for a, b in zip(sub.a, sub.b)]
    outp = [p for p, o in zip(partners, sub.out) if o]
    flag[s] = sum(c.pcr_group(p) != c.pcr_group(s) for p in outp) >= 2
res = dict(U=U, identity_rates=idr.rate.round(5).tolist(), identity_eligible=idr.eligible_mean.tolist(),
           clone_rate_min=float(cl.rate.min()), clone_rate_max=float(cl.rate.max()), n_clone_out=int(cl.out.sum()),
           clone_eligible_range=[int(cl.eligible_mean.min()), int(cl.eligible_mean.max())],
           min_anchor_rate=float(anc.rate.min()), min_anchor_by=dict(anc.groupby("cls").rate.min().round(5)),
           call_B2="ONE GENET" if one else "AMBIGUOUS", B2_flagged_stems=[s for s, f in flag.items() if f],
           B5_stem_median_rate={s: float(cl[(cl.a == s) | (cl.b == s)].rate.median()) for s in c.CLONE})
json.dump(res, open(f"{c.RUN}/05_B2_summary.json", "w"), indent=1, default=float)
print(json.dumps(res, indent=1, default=float))
print(O.groupby("cls")[["rate", "eligible_mean"]].describe().round(4).to_string())
