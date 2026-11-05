# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Shared definitions of the analysis run: input paths, sample sets (the 224 survey trees, the
# nine stems, WPR, C1), the read-count statistic and helpers. Imported by every script in scripts/finalize/ and
# scripts/dosage_clonality/.
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
"""Shared definitions for the clone-paper finalize run (R0-R5), spec analysis/2026-10-26_finalize_analyses_decision.md
(sha256 f24c1d5b...2b34), PREREG frozen in analysis/finalize_run/PREREG.tsv.

Sample sets (decision section 2.0), input paths, the SEGI-count statistic and the prediction-interval t.
Run every script with ~/miniconda3/envs/quast/bin/python3 -I.
"""
import csv
import hashlib
import math
import os
import statistics as st

import numpy as np
from scipy import stats

REPO = "/space/s1/sashanikolaeva/Scripts/GitHub/redwood-clone-paper"
# FINALIZE_SET=popgen224 (Sasha, 2026-10-26; DV14 in analysis/finalize_run_224/DEVIATIONS_224.md): reference panel and
# screen = exactly the popgen paper's 224 trees (inventory in_canonical_224, incl. 7_241-4); tree 286 is no longer a
# replicate pair. Default (unset) = the frozen pre-registered run, outputs in analysis/finalize_run/.
SET = os.environ.get("FINALIZE_SET", "prereg")
assert SET in ("prereg", "popgen224"), SET
OUT = f"{REPO}/analysis/finalize_run" if SET == "prereg" else f"{REPO}/analysis/finalize_run_224"
SCRATCH = "/space/s1/sashanikolaeva/clone_finalize_work/run"
FIGDIR = "/space/s1/sashanikolaeva/Scripts/GitHub/clone_paper_figures"
SEED = 20261026

INPUTS = {
    "segi_counts": "/space/s1/sashanikolaeva/to_narval_survey_calls/wild_readcounts_SEGI_MQ60_nochlor_proppaired_Aug15.csv",
    "vcf": "/space/s1/santang3/prg/results/variants/vcf/allScaffs_allRegions_filtered.vcf.gz",
    "snp_to_segi": "/space/s1/sashanikolaeva/Results/clone/snp_segi_alignment/snp_to_segi_chr.tsv",
    "snp_paf": "/space/s1/sashanikolaeva/Results/clone/snp_segi_alignment/snp_regions_to_segi.paf",
    "inventory": "/space/s1/sashanikolaeva/Scripts/GitHub/redwood-reference-genome/results/popmap/redwood_tree_inventory.tsv",
    "rep_alleles": "/space/s1/sashanikolaeva/Scripts/GitHub/redwood-reference-genome/results/het_cause/05_replicate_alleles.tsv.gz",
    "rep_concord": "/space/s1/sashanikolaeva/Scripts/GitHub/redwood-reference-genome/results/het_cause/05_replicate_concordance.tsv",
    "ibd_corrected": "/space/s1/sashanikolaeva/Results/Population_structure/April2025files/genetic_distance_matrix_IBD_corrected.csv",
    "ibd_uncorrected": "/space/s1/sashanikolaeva/Results/Population_structure/April2025files/genetic_distance_matrix_IBD.csv",
    "ibd_no_chr4": "/space/s1/sashanikolaeva/Results/clone/genetic_distance_matrix_IBD_no_chr4.csv",
    "fig2_pairs": "/space/s1/sashanikolaeva/Scripts/GitHub/clone_paper_figures/clonality/fig2_clonality_pairs.tsv",
    "fig2_summary": "/space/s1/sashanikolaeva/Scripts/GitHub/clone_paper_figures/clonality/fig2_clonality_summary.tsv",
    "fig2_script": "/space/s1/sashanikolaeva/Scripts/GitHub/redwood-clone-paper/scripts/figures/fig2_clonality.py",
    # secondary columns (rev05 definitions) and the V5 sheet (R2 v)
    "vcf_dp_per_chr": "/space/s1/sashanikolaeva/Results/clone/het_density_fullpop/depth_per_chr.tsv",
    "ad_het_per_chr": "/space/s1/sashanikolaeva/Results/clone/het_density_fullpop/het_density_ad_fullpop.tsv",
    "v5_sheet": "/space/s1/sashanikolaeva/to_narval_survey_calls/wild_metadata_Updated_coordinates.V5.csv",
    # library QC (R2-364 duplicate rate): Neale-alignment samtools stats via MultiQC
    "neale_samtools_stats": "/space/s1/santang3/prg/results/qc/multiqc/neale/multiqc_data/multiqc_samtools_stats.txt",
}

CH = [f"chr{i}" for i in range(1, 12)]
CLONE = ["13_WP1", "6_WP5", "31_c3", "31_WP2_1", "28_WP7", "28_WP6", "14_WP3", "22_WPd2", "33_c6"]
# rev01 REMOVE list (the VCF-script exclusions, depth_per_chr_vcf.py SAMPLES_TO_REMOVE)
REMOVE = set("""35_150 36_JSM4 22_MRD 13_286_1 16_416-1 17_111_2 36_M_spl 8_374 14_86 26_46H_sg 24_57_4 12_57 23_57-3
10_NQ2 27_NQH1 35_97_extra 10_368-2 21_63_5 21_45-2 19_382_3 30_157-0X 1_263_1 8_63-8 4_BH4 15_BH2 36_NP7 20_RW3 29_RW1
30_RW5 36_197 33_256 23_FR2 17_380 17_380-2 24_380-3 20_365 16_257-1 13_340-1 30_PG2 23_182 30_500 28_95 29_fl3 23_17
31_J5 32_J4""".split())
DEDUP_DROP = ["7_241-4", "10_279-1", "21_292_1", "28_385", "32_231", "32_63_7"]
TECH = [("22_MRD", "28_MRD"), ("16_416-1", "35_416_1"), ("13_286_1", "3_286-1")]
if SET == "popgen224":
    TECH = TECH[:2]   # popgen keeps only 3_286-1, so 286 is not a replicate pair
GAINCTRL = {"13_57_1": "chr7", "30_311_8": "chr8", "9_RSP1": "chr6"}
WPR, C1, S364 = "22_WPR", "16_C1", "14_364"
SOUTH_LAT = 37.82
LOCAL_KM = 60.0


def short(s):
    return s.replace("RedwoodCap_", "").replace("_SEGI_full", "")


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def load_segi_counts():
    rows = {}
    with open(INPUTS["segi_counts"]) as fh:
        r = csv.reader(fh)
        h = next(r)
        ix = {c: h.index(c) for c in CH}
        for x in r:
            rows[short(x[0])] = {c: float(x[ix[c]]) for c in CH}
    return rows


def load_inventory():
    with open(INPUTS["inventory"]) as fh:
        return {r["sample_id"]: r for r in csv.DictReader(fh, delimiter="\t")}


def hav(a, b):
    R = 6371.0
    la1, lo1, la2, lo2 = map(math.radians, (a[0], a[1], b[0], b[1]))
    return 2 * R * math.asin(math.sqrt(math.sin((la2 - la1) / 2) ** 2
                                       + math.cos(la1) * math.cos(la2) * math.sin((lo2 - lo1) / 2) ** 2))


def sample_sets():
    """Return dict of the fixed sets. Asserts are done by the caller (01_inputs_R0g.py)."""
    rc = load_segi_counts()
    inv = load_inventory()
    vcf238 = set()
    with open(INPUTS["vcf_dp_per_chr"]) as fh:
        next(fh)
        for line in fh:
            vcf238.add(short(line.split("\t")[0]))
    big = {s for s, d in rc.items() if sum(d.values()) > 1e6}
    ref229 = sorted(s for s in big if s not in CLONE and s not in REMOVE and s in vcf238)
    panel_all = [s for s in ref229 if s not in DEDUP_DROP]
    if SET == "popgen224":
        can = sorted(s for s in inv if inv[s].get("in_canonical_224") == "True")
        bad = [s for s in can if s not in big or s not in vcf238 or s not in rc]
        assert not bad, f"popgen224 trees lacking > 1 M SEGI reads or VCF DP: {bad}"
        assert len(can) == 224 and not set(can) & set(CLONE), len(can)
        panel_all = can
    coord = {s: (float(inv[s]["lat"]), float(inv[s]["lon"])) for s in inv if inv[s]["lat"] not in ("", "NA")}
    cen = (st.mean(coord[s][0] for s in CLONE), st.mean(coord[s][1] for s in CLONE))
    local = [s for s in panel_all if s in coord and hav(cen, coord[s]) <= LOCAL_KM]
    south = [s for s in panel_all if s in inv and inv[s]["region_lat"] == "S" and s in coord and coord[s][0] < SOUTH_LAT]
    # screen = samples tested in R2; plot = samples drawn in Fig 6 (prereg: every sample > 1 M reads)
    screen = list(panel_all) if SET == "popgen224" else sorted(s for s in big if s not in CLONE)
    plot = sorted(set(screen) | set(CLONE))
    return dict(rc=rc, inv=inv, big=big, ref229=ref229, PANEL_ALL=panel_all, PANEL_LOCAL=local, PANEL_SOUTH=south,
                centroid=cen, coord=coord, vcf238=vcf238, screen=screen, plot=plot)


def segi_ratio(rc, s, c):
    """R(s,c) = reads(c) / mean(reads of chr1-11 excluding c and chr4)."""
    return rc[s][c] / st.mean(rc[s][k] for k in CH if k not in (c, "chr4"))


def pred_t(x, ref):
    ref = np.asarray(ref, float)
    n = len(ref)
    mu, sd = ref.mean(), ref.std(ddof=1)
    t = (x - mu) / (sd * math.sqrt(1 + 1 / n))
    return t, n - 1, stats.t.cdf(t, n - 1), stats.t.sf(t, n - 1)   # t, df, p_lower, p_upper


def holm(ps):
    m = len(ps)
    order = sorted(range(m), key=lambda i: ps[i])
    adj, run = [0.0] * m, 0.0
    for k, i in enumerate(order):
        run = max(run, min(1.0, (m - k) * ps[i]))
        adj[i] = run
    return adj


def ensure_dirs():
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(SCRATCH, exist_ok=True)
