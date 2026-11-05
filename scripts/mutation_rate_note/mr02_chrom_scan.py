# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Per-chromosome re-run of the stem-private search and of the same filters on the technical-
# replicate pairs, with the depth of every stem (Supplementary Note S1: error calibration).
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
"""mr02 (from reviewer rv_calib_chr1.py, per-chromosome outputs). Reviewer check for Note S2 (rv_ = review). chr1 only, read-only on all inputs.
Re-runs the exact filters of somatic_mutations_mpileup.py (9-stem scan) and fp_calibration_bam.py
(replicate pairs) on chr1 of SEGI v2.0 in one mpileup over 13 BAMs (9 stems + 4 replicate libraries), and records:
 1. 9-stem scan: private calls per stem (must reproduce per_ramet_counts.tsv chr1 n_raw);
    for each call: alt reads of that base in stems that were NOT callable (DP<10, never checked by the scan),
    alt reads in the 4 replicate libraries (other trees), carrier depth.
 2. depth histogram of each stem at callable sites (for detection probability under the real filters).
 3. replicate pairs: pair-private calls split by library (fp_calibration_bam.py counts both libraries together),
    and whether the allele is seen (>=2 reads) in any of the other 11 libraries (population allele -> dropout/paralog,
    not a library-specific error).
Output: rv_calib_chr1.tsv (summary), rv_calib_chr1_private.tsv, rv_calib_chr1_pairs.tsv, rv_calib_chr1_depth.tsv."""
import bisect, csv, re, subprocess, sys
from collections import defaultdict, Counter
OUT = "/space/s1/sashanikolaeva/Scripts/GitHub/redwood-clone-paper/analysis/mutation_rate_note"
B = "/space/s1/sashanikolaeva/Results/clone/nquire/bam"
REF = "/space/s1/sashanikolaeva/Data/SEGI_genome/SEGI.2.0.fa"
NUPT = "/space/s1/sashanikolaeva/Results/clone/nquire/nupt_mask.bed"
STEMS = ["13_WP1", "6_WP5", "31_c3", "31_WP2_1", "28_WP7", "28_WP6", "14_WP3", "22_WPd2", "33_c6"]
REPS = ["22_MRD", "28_MRD", "16_416-1", "35_416_1"]
PAIRS = {"MRD": (9, 10), "416": (11, 12)}  # indices into ALL
ALL = STEMS + REPS
CHROM = sys.argv[1] if len(sys.argv) > 1 else "chr1"
MIN_DP, MIN_ALT, MIN_AF, MAX_SB, MIN_CALLED = 10, 2, 0.10, 0.95, 7

mask = []
for l in open(NUPT):
    p = l.split("\t")
    if p[0] == CHROM: mask.append((int(p[1]), int(p[2])))
mask.sort(); starts = [m[0] for m in mask]
def in_nupt(pos):  # same (1-based pos vs BED inclusive) logic as the original scripts
    i = bisect.bisect_right(starts, pos) - 1
    return i >= 0 and mask[i][0] <= pos <= mask[i][1]

def parse(s):  # identical to parse_pileup_bases in the original scripts
    s = re.sub(r'\^.', '', s).replace('$', '')
    s = re.sub(r'[+-](\d+)', lambda m: ' ' * (int(m.group(1)) + len(m.group(0))), s)
    s = re.sub(r'\s', '', s)
    f = Counter(c for c in s if c in 'ACGT'); r = Counter(c.upper() for c in s if c in 'acgt')
    alt = {b: f[b] + r[b] for b in set(f) | set(r)}
    return alt, sum(f.values()), sum(r.values())

def private_call(dp, alt, nf, nr):
    if dp < MIN_DP or not alt: return None
    b = max(alt, key=alt.get); n = alt[b]
    if n < MIN_ALT or n / dp < MIN_AF: return None
    if nf + nr > 0 and max(nf, nr) / (nf + nr) > MAX_SB: return None
    return b, n

cmd = ["samtools", "mpileup", "-f", REF, "-B", "-q", "20", "-Q", "20", "--no-output-ins", "--no-output-del", "-r", CHROM] + [f"{B}/{s}.no_chlor.sorted.bam" for s in ALL]
proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
n_pos = n_nupt = n_low9 = 0
callable9 = Counter(); depth_hist = defaultdict(Counter)
priv_rows = []; pair_rows = []; pair_callable = Counter()
for line in proc.stdout:
    p = line.rstrip("\n").split("\t"); pos = int(p[1]); n_pos += 1
    if in_nupt(pos): n_nupt += 1; continue
    dps = [int(p[3 + 3 * i]) for i in range(len(ALL))]
    stem_ok = [d >= MIN_DP for d in dps[:9]]
    pair_ok = {k: dps[a] >= MIN_DP and dps[b] >= MIN_DP for k, (a, b) in PAIRS.items()}
    do9 = sum(stem_ok) >= MIN_CALLED
    if not do9: n_low9 += 1
    if not (do9 or any(pair_ok.values())): continue
    parsed = [parse(p[4 + 3 * i]) for i in range(len(ALL))]
    if do9:
        for i in range(9):
            if stem_ok[i]:
                callable9[STEMS[i]] += 1; depth_hist[STEMS[i]][min(dps[i], 200)] += 1
        for i in range(9):
            if not stem_ok[i]: continue
            c = private_call(dps[i], *parsed[i])
            if not c: continue
            b, n = c
            if any(parsed[j][0].get(b, 0) > 0 for j in range(9) if j != i and stem_ok[j]): continue
            nc_alt = sum(parsed[j][0].get(b, 0) for j in range(9) if j != i and not stem_ok[j])
            n_nc = sum(1 for j in range(9) if not stem_ok[j])
            rep_alt = [parsed[j][0].get(b, 0) for j in range(9, 13)]
            priv_rows.append([STEMS[i], CHROM, pos, b, n, dps[i], n_nc, nc_alt, sum(rep_alt), sum(1 for x in rep_alt if x >= 2)])
    for k, (a, bb) in PAIRS.items():
        if not pair_ok[k]: continue
        pair_callable[k] += 1
        for i, j in ((a, bb), (bb, a)):
            c = private_call(dps[i], *parsed[i])
            if not c: continue
            b, n = c
            if parsed[j][0].get(b, 0) > 0: continue
            others = [x for x in range(13) if x not in (a, bb)]
            n_other2 = sum(1 for x in others if parsed[x][0].get(b, 0) >= 2)
            n_other1 = sum(1 for x in others if parsed[x][0].get(b, 0) >= 1)
            pair_rows.append([k, ALL[i], CHROM, pos, b, n, dps[i], dps[j], n_other1, n_other2])
proc.wait()

with open(f"{OUT}/mr02_chrom/{CHROM}_private.tsv", "w") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(["stem", "chrom", "pos", "alt", "alt_reads", "depth", "n_stems_not_callable", "alt_reads_in_not_callable_stems", "alt_reads_in_4_replicate_libs", "n_replicate_libs_alt_ge2"]); w.writerows(priv_rows)
with open(f"{OUT}/mr02_chrom/{CHROM}_pairs.tsv", "w") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(["pair", "library", "chrom", "pos", "alt", "alt_reads", "depth", "depth_other_replicate", "n_other11_alt_ge1", "n_other11_alt_ge2"]); w.writerows(pair_rows)
with open(f"{OUT}/mr02_chrom/{CHROM}_depth.tsv", "w") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n"); w.writerow(["stem", "depth_capped200", "n_sites"])
    for s in STEMS:
        for d in sorted(depth_hist[s]): w.writerow([s, d, depth_hist[s][d]])
S = []
def rec(q, v, scope, unit): S.append((q, v, scope, unit)); print(q, v, scope, unit, sep="\t")
rec("positions", n_pos, f"{CHROM}, positions output by mpileup over 13 BAMs", "positions")
rec("nupt", n_nupt, CHROM, "positions"); rec("low_depth_9stem", n_low9, "fewer than 7 of 9 stems DP>=10 (non-NUPT)", "positions")
for s in STEMS:
    rec(f"callable_{s}", callable9[s], "DP>=10 at sites with >=7/9 stems DP>=10", "sites")
    rec(f"private_{s}", sum(1 for r in priv_rows if r[0] == s), "9-stem scan, exact original filters", "calls")
rec("private_total", len(priv_rows), CHROM, "calls")
rec("private_with_noncallable_stem", sum(1 for r in priv_rows if r[6] > 0), "calls where >=1 other stem had DP<10 (not checked)", "calls")
rec("private_with_alt_in_noncallable_stem", sum(1 for r in priv_rows if r[7] > 0), "calls where an unchecked stem (DP<10) carried >=1 alt read", "calls")
rec("private_with_alt_ge2_in_replicate_lib", sum(1 for r in priv_rows if r[9] > 0), "calls whose allele has >=2 reads in >=1 of the 4 replicate libraries (other trees)", "calls")
for k in PAIRS:
    rec(f"pair_callable_{k}", pair_callable[k], "both libraries DP>=10", "sites")
    for lib in (ALL[PAIRS[k][0]], ALL[PAIRS[k][1]]):
        rr = [r for r in pair_rows if r[0] == k and r[1] == lib]
        rec(f"pair_private_{k}_{lib}", len(rr), "allele in this library, 0 reads in its replicate", "calls")
        rec(f"pair_private_{k}_{lib}_in_other11_ge2", sum(1 for r in rr if r[9] > 0), "of those, allele >=2 reads in >=1 of the other 11 libraries", "calls")
        rec(f"pair_private_{k}_{lib}_absent_other11", sum(1 for r in rr if r[8] == 0), "of those, 0 reads in all other 11 libraries", "calls")
with open(f"{OUT}/mr02_chrom/{CHROM}_summary.tsv", "w") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n"); w.writerow(["quantity", "value", "scope", "unit"]); w.writerows(S)
print("done")
