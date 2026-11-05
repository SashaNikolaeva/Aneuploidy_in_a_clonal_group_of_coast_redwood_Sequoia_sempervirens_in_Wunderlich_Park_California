#!/usr/bin/env python3
# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Input checksums and the SNP sets of the allele-dosage test (1,584 SNPs, 180 on chromosome
# 4) and of the lattice model.
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
"""01: input hashes (V0) and site sets for A9 (primary) and the lattice (L0), per PREREG_v2. No allele fraction of the
stems is used here except clone heterozygosity (allele presence) for the secondary lattice set L0.
Outputs: RUN/01_input_hashes.tsv, RUN/01_site_counts.tsv, RUN/sites.npz; stdout -> RUN/01_inputs_sites.log
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dc_common as c  # noqa: E402

arr = c.load_array()
inv = c.load_inventory()
ix, P224, R = c.sets(arr["samples"], inv)
chrom, DP, RO, AO = arr["chrom"], arr["DP"], arr["RO"], arr["AO"]
with open(f"{c.RUN}/01_input_hashes.tsv", "w") as f:
    for lab, p in [("npz", c.NPZ), ("extract", c.EXTRACT), ("inventory", c.INV), ("ibd_corrected", c.MAT_C),
                   ("ibd_uncorrected", c.MAT_U), ("prereg_v2", f"{c.REPO}/analysis/dosage_clonality_run/PREREG_v2.tsv")]:
        f.write(f"{lab}\t{p}\t{c.sha256(p)}\n")
print(open(f"{c.RUN}/01_input_hashes.tsv").read())
for s in c.CLONE + R:
    g = c.pcr_group(s)
    assert g == int(inv.loc[s, "pcr_group"]), (s, g, inv.loc[s, "pcr_group"])
print("PCR group rule agrees with the inventory for the stems and R")

meta = c.site_meta(arr)
pur = c.purity_ok(meta, chrom)
ri = np.array([ix[s] for s in R])
cov = DP[:, ri] >= 15
carr = (cov & (RO[:, ri] >= 2) & (AO[:, ri] >= 2)).sum(1)
ncov = cov.sum(1)
dmask = np.where(cov, DP[:, ri].astype(float), np.nan)
with np.errstate(all="ignore"):
    pmed = np.nanmedian(dmask, axis=1)
pmed9 = np.nan_to_num(pmed) * 9
fixed = (ncov >= 50) & (carr >= 0.9 * ncov) & (carr >= 10) & pur
cap = np.percentile(pmed9[fixed & (chrom != "chr4") & (pmed9 >= 100)], 99)
A9 = fixed & (pmed9 >= 100) & (pmed9 <= cap)
A9_200 = A9 & (pmed9 >= 200)
print(f"A9: cap (99th pct of panel median DP x9 over qualifying non-chr4 sites, after floor) = {cap:.0f}")

# lattice L0: clone-het (allele presence), same panel depth floor/cap, >= 10 panel carriers among the 224
ci = [ix[s] for s in c.CLONE]
ro, ao = RO[:, ci].sum(1), AO[:, ci].sum(1)
het = (np.minimum(ro, ao) >= 5) & ((AO[:, ci] >= 1).sum(1) >= 2) & ((RO[:, ci] >= 1).sum(1) >= 2)
pi = np.array([ix[s] for s in P224])
cov224 = DP[:, pi] >= 15
carr224 = (cov224 & (RO[:, pi] >= 2) & (AO[:, pi] >= 2)).sum(1)
L0 = het & (pmed9 >= 100) & (pmed9 <= cap) & (carr224 >= 10) & pur

rows = []
for ch in c.CH:
    m = chrom == ch
    rows.append(dict(chrom=ch, A9=int((A9 & m).sum()), A9_single_locus=int((A9 & m & (meta["nloci"] == 1)).sum()),
                     A9_depth200=int((A9_200 & m).sum()), L0=int((L0 & m).sum()),
                     L0_single_locus=int((L0 & m & (meta["nloci"] == 1)).sum())))
T = pd.DataFrame(rows)
T.to_csv(f"{c.RUN}/01_site_counts.tsv", sep="\t", index=False)
print(T.to_string(index=False))
print("totals:", T.drop(columns="chrom").sum().to_dict())
np.savez_compressed(f"{c.RUN}/sites.npz", A9=A9, A9_200=A9_200, L0=L0, pmed9=pmed9, cap=cap)
