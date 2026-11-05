# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Shared inputs of the allele-dosage and batch-aware clonality analyses (sample sets, PCR
# groups, SNP site metadata).
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
"""Shared inputs for the dosage + clonality run (PREREG_v2, analysis/dosage_clonality_run/PREREG_v2.tsv).

Loads the SEGI-mapped survey SNP array directly from the npz with a sha256 assertion (V0), the tree inventory, sample
sets (popgen 224, reference R = 224 minus aneuploid trees, the 9 stems, identity pairs), PCR groups, and per-SNP site
metadata (Neale scaffold/position, Neale 1-Mb block purity, single- vs multi-Neale-locus SEGI 1-kb loci, panel depth
ratio), which is built once from the step-03 extract and cached in the run directory.
"""
import gzip
import hashlib
import os
import sys

import numpy as np
import pandas as pd

REPO = "/space/s1/sashanikolaeva/Scripts/GitHub/redwood-clone-paper"
sys.path.insert(0, f"{REPO}/scripts/finalize")
import vcfcache  # noqa: E402  (only path constants and _positions() are used)

S1 = "/space/s1/sashanikolaeva"
RUN = os.environ.get("DC_RUN", f"{REPO}/analysis/dosage_clonality_run")  # DC_RUN: reproducibility re-run directory
FIG = f"{S1}/Scripts/GitHub/clone_paper_figures/dosage_clonality"
NPZ = vcfcache.CACHE
NPZ_SHA = "401752d614304e3e4feb9fedf67a16753bc13324246814ac3912b1de51cb3c4e"
EXTRACT = vcfcache.SRC
INV = f"{S1}/Scripts/GitHub/redwood-reference-genome/results/popmap/redwood_tree_inventory.tsv"
MAT_C = f"{S1}/Results/Population_structure/April2025files/genetic_distance_matrix_IBD_corrected.csv"
MAT_U = f"{S1}/Results/Population_structure/April2025files/genetic_distance_matrix_IBD.csv"
SEED = 20261027
CH = [f"chr{i}" for i in range(1, 12)]
CLONE = ["13_WP1", "6_WP5", "31_c3", "31_WP2_1", "28_WP7", "28_WP6", "14_WP3", "22_WPd2", "33_c6"]
PILOT_STEMS = ["14_WP3", "22_WPd2"]
ANEU = ["13_57_1", "30_311_8", "9_RSP1", "30_NQ3", "14_364"]
IDENTITY = [("2_BH3", "15_BH2"), ("21_17-1", "23_17"), ("22_MRD", "28_MRD"), ("21_279sg", "10_279-1"),
            ("1_256S", "33_256"), ("12_292", "21_292_1"), ("16_416-1", "35_416_1")]
WPR, C1 = "22_WPR", "16_C1"
SOUTH = 37.82
os.makedirs(RUN, exist_ok=True)
os.makedirs(FIG, exist_ok=True)


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def pcr_group(sample):
    """PCR group from the pool number (first field of the sample name), pool -> group rule of the popgen paper."""
    pool = int(sample.split("_")[0])
    for g, pools in {1: [14, 17, 18, 22], 2: [13, 15, 16, 19, 20, 21, 23, 24], 3: [26, 27, 29, 31, 32, 33, 35, 36],
                     4: [1, 2, 3, 4, 5, 6, 7, 8], 5: [9, 10, 11, 12, 25, 28, 30, 34]}.items():
        if pool in pools:
            return g
    return None


def load_array():
    assert sha256(NPZ) == NPZ_SHA, "npz sha256 mismatch"
    z = np.load(NPZ)
    return dict(samples=list(z["samples"]), chrom=z["chrom"].astype(str), pos=z["pos"], DP=z["DP"], RO=z["RO"],
                AO=z["AO"])


def load_inventory():
    return pd.read_csv(INV, sep="\t", low_memory=False).set_index("sample_id")


def sets(samples, inv):
    ix = {s: i for i, s in enumerate(samples)}
    P224 = [s for s in inv.index[inv.in_canonical_224 == True] if s in ix]  # noqa: E712
    assert len(P224) == 224, len(P224)
    R = [s for s in P224 if s not in ANEU]
    for s in CLONE + [WPR, C1] + [x for p in IDENTITY for x in p]:
        assert s in ix, s
    return ix, P224, R


def site_meta(arr):
    """Neale key per cache row (same filter/order as vcfcache.build), Neale 1-Mb block purity, Neale-locus count per
    SEGI 1-kb locus. Cached in RUN/site_meta.npz."""
    fn = f"{RUN}/site_meta.npz"
    if os.path.exists(fn):
        m = np.load(fn)
        return {k: m[k] for k in m.files}
    best = vcfcache._positions()
    scaf, npos = [], []
    with gzip.open(EXTRACT, "rt") as fh:
        fh.readline()
        for line in fh:
            a, b = line.split("\t", 2)[:2]
            if f"{a}_{b}" in best:
                scaf.append(a); npos.append(int(b))
    scaf, npos = np.array(scaf), np.array(npos)
    chrom = arr["chrom"]
    assert len(scaf) == len(chrom)
    blk = np.char.add(np.char.add(scaf, ":"), (npos // 1_000_000).astype(str))
    df = pd.DataFrame(dict(blk=blk, chrom=chrom))
    tab = df.groupby(["blk", "chrom"]).size().rename("n").reset_index()
    tot = tab.groupby("blk").n.transform("sum")
    tab["share"] = tab.n / tot
    maj = tab.sort_values("share", ascending=False).drop_duplicates("blk").set_index("blk")
    purity = maj.share.reindex(blk).values
    majchrom = maj.chrom.reindex(blk).values.astype(str)
    nlocus = np.char.add(np.char.add(scaf, ":"), (npos // 100_000).astype(str))
    sbin = np.char.add(np.char.add(chrom, ":"), (arr["pos"] // 1000).astype(str))
    k = pd.DataFrame(dict(sbin=sbin, nl=nlocus)).groupby("sbin").nl.nunique()
    nloci = k.reindex(sbin).values
    out = dict(scaf=scaf, npos=npos, purity=purity, majchrom=majchrom, nloci=nloci)
    np.savez_compressed(fn, **out)
    return out


def purity_ok(meta, chrom):
    return (meta["purity"] >= 0.90) & (meta["majchrom"] == chrom)


def panel_ratio(DP, cols):
    d = DP[:, cols].astype(float)
    med = np.median(d, axis=0)
    return np.median(d / np.maximum(med, 1), axis=1)


def hav_km(la1, lo1, la2, lo2):
    r = np.pi / 180
    h = np.sin((la2 - la1) * r / 2) ** 2 + np.cos(la1 * r) * np.cos(la2 * r) * np.sin((lo2 - lo1) * r / 2) ** 2
    return 2 * 6371.0 * np.arcsin(np.sqrt(h))
