#!/usr/bin/env python3
# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Secondary lattice mixture model of allele fractions with its identifiability gate (Results:
# not identifiable, not applied to chromosome 4).
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
"""03: SECONDARY lattice mixture with staged unblinding (PREREG_v2 L0-L3).

Stage 1 fits only the ten non-chr4 chromosomes and applies the identifiability gate (bootstrap SD of W5 <= 0.03 on every
non-chr4 chromosome AND no non-chr4 W5 2.5th percentile > 0.03). chr4 is fitted (stage 2) only if the gate passes.
Model: beta-binomial mixture, 11 lattice peaks + uniform background; r grid 0.80-1.25 (0.05), rho grid {0.002,0.005,0.01,
0.02,0.04} by profile likelihood; EM from 20 Dirichlet starts (seed 20261027) + 1 equal-weight start, <= 2,000 iterations,
tol 1e-8; 200 block-bootstrap replicates over 10-Mb SEGI windows with weights refitted at fixed (r, rho).
Outputs: RUN/03_lattice_stage1.tsv, RUN/03_lattice_gate.json (+ stage-2 files if PASS); stdout -> RUN/03_lattice.log
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.special import betaln, gammaln

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dc_common as c  # noqa: E402

PEAKS = np.array([1/6, 1/5, 1/4, 1/3, 2/5, 1/2, 3/5, 2/3, 3/4, 4/5, 5/6])
FIVE = np.array([1, 4, 6, 9])
RGRID = np.round(np.arange(0.80, 1.2501, 0.05), 2)
RHOGRID = [0.002, 0.005, 0.01, 0.02, 0.04]
arr = c.load_array()
inv = c.load_inventory()
ix, P224, R = c.sets(arr["samples"], inv)
S = np.load(f"{c.RUN}/sites.npz")
m0 = S["L0"]
ci = [ix[s] for s in c.CLONE]
k_all = arr["AO"][m0][:, ci].sum(1).astype(float)
n_all = (arr["AO"][m0][:, ci] + arr["RO"][m0][:, ci]).sum(1).astype(float)
chrom = arr["chrom"][m0]
win = (arr["pos"][m0] // 10_000_000).astype(int)


def loglik_matrix(k, n, r, rho):
    mu = PEAKS / (PEAKS + r * (1 - PEAKS))
    a, b = mu * (1 / rho - 1), (1 - mu) * (1 / rho - 1)
    lc = gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)
    L = lc[:, None] + betaln(k[:, None] + a, (n - k)[:, None] + b) - betaln(a, b)
    bg = -np.log(n + 1)
    return np.column_stack([L, bg])


def em(L, W0, iters=2000, tol=1e-8, fixed_zero=None):
    """L: sites x K log-lik; W0: starts x K weights. Returns best (loglik, weights)."""
    mx = L.max(1, keepdims=True)
    E = np.exp(L - mx)
    W = W0.copy()
    if fixed_zero is not None:
        W[:, fixed_zero] = 0
        W /= W.sum(1, keepdims=True)
    prev = np.full(W.shape[0], -np.inf)
    for _ in range(iters):
        M = E @ W.T                       # sites x starts
        ll = np.log(M).sum(0) + mx.sum()
        R_ = E[:, None, :] * W[None, :, :] / M[:, :, None]
        W = R_.mean(0)
        if fixed_zero is not None:
            W[:, fixed_zero] = 0
            W /= W.sum(1, keepdims=True)
        if np.all(np.abs(ll - prev) < tol):
            break
        prev = ll
    M = E @ W.T
    ll = np.log(M).sum(0) + mx.sum()
    j = int(np.argmax(ll))
    return float(ll[j]), W[j]


def starts(K, rng):
    return np.vstack([np.full(K, 1 / K), rng.dirichlet(np.ones(K), 20)])


def fit_chrom(ch, rng):
    m = chrom == ch
    k, n = k_all[m], n_all[m]
    best = None
    for r in RGRID:
        for rho in RHOGRID:
            L = loglik_matrix(k, n, r, rho)
            ll, w = em(L, starts(12, rng))
            if best is None or ll > best[0]:
                best = (ll, w, r, rho)
    ll, w, r, rho = best
    L = loglik_matrix(k, n, r, rho)
    lln, wn = em(L, w[None, :], fixed_zero=FIVE)
    lln = min(lln, ll)
    # block bootstrap over 10-Mb windows, weights refitted at fixed (r, rho) from the best full-fit weights
    wins = win[m]
    uw = np.unique(wins)
    W5b = []
    for _ in range(200):
        pick = rng.choice(uw, len(uw), replace=True)
        idx = np.concatenate([np.where(wins == u)[0] for u in pick])
        _, wb = em(L[idx], w[None, :])
        W5b.append(wb[FIVE].sum())
    W5b = np.array(W5b)
    return dict(chrom=ch, sites=int(m.sum()), r=r, rho=rho, loglik_full=ll, loglik_null=lln, LLR=2 * (ll - lln),
                W5=float(w[FIVE].sum()), W_background=float(w[11]), W6u=float(w[[0, 10]].sum()),
                W5_boot_sd=float(W5b.std(ddof=1)), W5_lo=float(np.percentile(W5b, 2.5)), W5_hi=float(np.percentile(W5b, 97.5)))


rng = np.random.default_rng(c.SEED)
rows = [fit_chrom(ch, rng) for ch in c.CH if ch != "chr4"]
T1 = pd.DataFrame(rows)
T1.to_csv(f"{c.RUN}/03_lattice_stage1.tsv", sep="\t", index=False, float_format="%.5f")
print(T1.to_string(index=False))
gate = bool((T1.W5_boot_sd <= 0.03).all() and not (T1.W5_lo > 0.03).any())
res = dict(gate="PASS" if gate else "FAIL", max_boot_sd=float(T1.W5_boot_sd.max()),
           n_flagged_lo_gt_0p03=int((T1.W5_lo > 0.03).sum()))
print("IDENTIFIABILITY GATE:", res)
if gate:
    r4 = fit_chrom("chr4", rng)
    dW5 = r4["W5"] - float(T1.W5.median())
    res.update(stage2=r4, dW5=dW5)
    print("stage 2 chr4:", r4, "dW5", dW5)
else:
    res["stage2"] = "not run: lattice not identifiable in these data (chr4 not fitted)"
json.dump(res, open(f"{c.RUN}/03_lattice_gate.json", "w"), indent=1, default=float)
