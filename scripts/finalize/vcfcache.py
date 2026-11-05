# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Converts the VCF extract of 03_vcf_pass.sh into numpy arrays with giant sequoia chromosome
# and position per SNP.
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
"""Parse the step-4 VCF extract once into numpy arrays (scratch cache), with SEGI chromosome and position per SNP.

SEGI position = best PAF hit of the 601-bp flank on its assigned SEGI chromosome (max MAPQ, then matches);
pos = tstart + (300 - qstart) on '+', tend - (300 - qstart) on '-'  (d02 code).
Only SNPs assigned to SEGI chr1-11 are kept. Missing values ('.') -> 0.
load() returns dict(samples, chrom (str array), pos (int), DP, RO, AO (int32, SNP x sample)).
"""
import gzip
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import CH, INPUTS, SCRATCH, short  # noqa

SRC = f"{SCRATCH}/vcf_dp_ad.tsv.gz"
CACHE = f"{SCRATCH}/vcf_dp_ad_segi.npz"


def _positions():
    assigned = {}
    with open(INPUTS["snp_to_segi"]) as fh:
        next(fh)
        for line in fh:
            s, p, c = line.rstrip("\n").split("\t")
            if c in CH:
                assigned[f"{s}_{p}"] = c
    best = {}
    with open(INPUTS["snp_paf"]) as fh:
        for line in fh:
            f = line.split("\t")
            q = f[0]
            c = assigned.get(q)
            if c is None or f[5] != c:
                continue
            key = (int(f[11]), int(f[9]))
            if q not in best or key > best[q][0]:
                qs, ts, te = int(f[2]), int(f[7]), int(f[8])
                pos = ts + (300 - qs) if f[4] == "+" else te - (300 - qs)
                best[q] = (key, c, pos)
    return best


def build():
    best = _positions()
    chrom, pos, DP, RO, AO = [], [], [], [], []
    with gzip.open(SRC, "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ns = (len(hdr) - 2) // 2   # DP fields, then one "RO,AO" field per sample
        samples = [short(h.split("]", 1)[1].split(":")[0]) for h in hdr[2:2 + ns]]
        for line in fh:
            f = line.rstrip("\n").split("\t")
            b = best.get(f"{f[0]}_{f[1]}")
            if b is None:
                continue
            chrom.append(b[1]); pos.append(b[2])
            DP.append([int(x) if x != "." else 0 for x in f[2:2 + ns]])
            ro, ao = [], []
            assert len(f) == 2 + 2 * ns
            for x in f[2 + ns:]:
                y = x.split(",")
                ro.append(int(y[0]) if y[0] != "." else 0)
                ao.append(int(y[1]) if len(y) > 1 and y[1] != "." else 0)
            RO.append(ro); AO.append(ao)
    np.savez(CACHE, samples=np.array(samples), chrom=np.array(chrom), pos=np.array(pos, dtype=np.int64),
             DP=np.array(DP, dtype=np.int32), RO=np.array(RO, dtype=np.int32), AO=np.array(AO, dtype=np.int32))
    print(f"cache built: {len(chrom)} SNPs x {len(samples)} samples -> {CACHE}")


def load():
    if not os.path.exists(CACHE) or os.path.getmtime(CACHE) < os.path.getmtime(SRC):
        build()
    z = np.load(CACHE)
    return {k: z[k] for k in z.files}


if __name__ == "__main__":
    build()
