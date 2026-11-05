# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Combines the outputs above into every number of Supplementary Note S1 and its table of
# expected mutation counts.
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
"""Note S2: combine the 2026-10-27 re-run (current BAMs) over chr1-11 and derive every number used in the note.
Inputs: mr02_chrom/chr*_{summary,private,pairs,depth}.tsv (mr02_chrom_scan.py) and the original-script re-run
Results/clone/somatic_mutations_mpileup_2026-10-27/per_ramet_counts.tsv (mr02_scan_original.py).
Output: mr03_numbers.tsv + stdout (mr03.log)."""
import csv, math, os
from collections import defaultdict
OUT = os.path.dirname(os.path.abspath(__file__)); C = f"{OUT}/mr02_chrom"
NEW = "/space/s1/sashanikolaeva/Results/clone/somatic_mutations_mpileup_2026-10-27"
CHR = [f"chr{i}" for i in range(1, 12)]
G = {**{p: 1 for p in (14, 17, 18, 22)}, **{p: 2 for p in (13, 15, 16, 19, 20, 21, 23, 24)},
     **{p: 3 for p in (26, 27, 29, 31, 32, 33, 35, 36)}, **{p: 4 for p in range(1, 9)}, **{p: 5 for p in (9, 10, 11, 12, 25, 28, 30, 34)}}
rows = []
def rec(q, v, scope, unit, src): rows.append((q, v, scope, unit, "V(run 2026-10-27)", src)); print(q, v, scope, unit, sep="\t")
def tsv(p): return list(csv.DictReader(open(p), delimiter="\t"))
summ = defaultdict(int)
for c in CHR:
    for r in tsv(f"{C}/{c}_summary.tsv"): summ[r["quantity"]] += int(r["value"])
stems = sorted({k[len("private_"):] for k in summ if k.startswith("private_") and k[len("private_"):] in
                {q[len("callable_"):] for q in summ if q.startswith("callable_")}})
rec("positions", summ["positions"], "chr1-11, mpileup over 13 BAMs (positions with reads)", "positions", "mr02_chrom")
rec("nupt_masked", summ["nupt"], "", "positions", "mr02_chrom")
rec("low_depth", summ["low_depth_9stem"], "fewer than 7 of 9 stems with DP>=10", "positions", "mr02_chrom")
rec("low_depth_share", round(summ["low_depth_9stem"] / summ["positions"], 3), "", "fraction", "mr02_chrom")
# Denominators printed in Note S2 come from the original-script run over the 9 stem BAMs (the 13-BAM re-implementation
# above also counts positions covered only by the 4 replicate libraries). Added 2026-10-29.
import re
_m = re.findall(r"Done: ([\d,]+) positions \| NUPT filtered: ([\d,]+) \| Low depth: ([\d,]+) \| Private variants: (\d+)",
                open(f"{OUT}/mr02_scan_original.log").read())[-1]
_pos, _nupt, _low, _priv = (int(v.replace(",", "")) for v in _m)
rec("positions_original_scan", _pos, "chr1-11, original script over the 9 stem BAMs (positions with reads)", "positions", "mr02_scan_original.log")
rec("low_depth_original_scan", _low, "fewer than 7 of 9 stems with DP>=10", "positions", "mr02_scan_original.log")
rec("low_depth_share_original_scan", round(_low / _pos, 3), "", "fraction", "mr02_scan_original.log")
orig = defaultdict(int)
for r in tsv(f"{NEW}/per_ramet_counts.tsv"): orig[r["ramet"]] += int(r["n_raw"])
for s in stems:
    rec(f"private_{s}", orig[s], f"PCR group {G[int(s.split('_')[0])]}; original script (authoritative); 13-library re-implementation gives {summ[f'private_{s}']}", "calls", "per_ramet_counts.tsv (+ mr02_chrom)")
rec("private_total", sum(orig.values()), "9 stems, original script", "calls", "per_ramet_counts.tsv")
rec("private_total_reimpl", sum(summ[f"private_{s}"] for s in stems), "13-library re-implementation (diagnostics below refer to these calls)", "calls", "mr02_chrom")
tot = sum(summ[f"private_{s}"] for s in stems)
for s in stems: summ[f"private_{s}"] = orig[s]  # per-stem counts and rates below use the authoritative counts
cal = {s: summ[f"callable_{s}"] for s in stems}
rec("callable_per_stem_Mb", f"{min(cal.values())/1e6:.1f}-{max(cal.values())/1e6:.1f}", "DP>=10 in the stem at sites with >=7/9 stems DP>=10", "Mb", "")
rate = {s: summ[f"private_{s}"] / cal[s] * 1e6 for s in stems}
g3 = [s for s in stems if G[int(s.split('_')[0])] == 3]; ot = [s for s in stems if s not in g3]
for lab, ss in (("group3", g3), ("other_groups", ot)):
    rec(f"private_range_{lab}", f"{min(summ['private_'+s] for s in ss)}-{max(summ['private_'+s] for s in ss)} (n={len(ss)})", "", "calls", "")
    rec(f"private_per_Mb_{lab}", f"{min(rate[s] for s in ss):.1f}-{max(rate[s] for s in ss):.1f}", "", "per Mb", "")
grp1 = [s for s in stems if G[int(s.split('_')[0])] == 1]
rec("private_group1", ", ".join(f"{s} {summ['private_'+s]}" for s in grp1), "", "calls", "")
rec("private_with_unchecked_stem", f"{summ['private_with_noncallable_stem']} ({summ['private_with_noncallable_stem']/tot:.0%})", ">=1 other stem DP<10, not checked", "calls", "")
rec("private_with_alt_in_unchecked_stem", summ["private_with_alt_in_noncallable_stem"], "unchecked stem carried >=1 alt read", "calls", "")
rec("private_recurring_in_replicate_libs", f"{summ['private_with_alt_ge2_in_replicate_lib']} ({summ['private_with_alt_ge2_in_replicate_lib']/tot:.0%})", "allele >=2 reads in >=1 of 4 libraries of other trees", "calls", "")
# replicate pairs
lib_rates, lib_abs, lib_share = [], [], []
for k in ("MRD", "416"):
    pc = summ[f"pair_callable_{k}"]; libs = [q[len(f"pair_private_{k}_"):] for q in summ if q.startswith(f"pair_private_{k}_") and not q.endswith(("_ge2", "_other11"))]
    pair_tot = sum(summ[f"pair_private_{k}_{l}"] for l in libs)
    rec(f"pair_{k}_callable_Mb", round(pc / 1e6, 2), "both libraries DP>=10", "Mb", "")
    rec(f"pair_{k}_private_both_libs", f"{pair_tot} ({pair_tot/pc*1e6:.0f} per Mb per pair)", "allele in one library, 0 reads in the other", "calls", "")
    for l in libs:
        n = summ[f"pair_private_{k}_{l}"]; lib_rates.append(n / pc * 1e6)
        lib_share.append(summ[f"pair_private_{k}_{l}_in_other11_ge2"] / n); lib_abs.append(summ[f"pair_private_{k}_{l}_absent_other11"] / pc * 1e6)
        rec(f"pair_{k}_{l}", f"{n}; {n/pc*1e6:.0f}/Mb; in other trees {lib_share[-1]:.0%}; absent from all 11 others {lib_abs[-1]:.0f}/Mb", "", "", "")
mr = sum(lib_rates) / len(lib_rates)
rec("pair_rate_per_library", f"{min(lib_rates):.0f}-{max(lib_rates):.0f} (mean {mr:.0f})", "", "per Mb", "")
rec("pair_share_in_other_trees", f"{min(lib_share):.0%}-{max(lib_share):.0%}", "pair-private allele with >=2 reads in >=1 of other 11 libraries", "fraction", "")
rec("pair_absent_other11_per_Mb", f"{min(lib_abs):.0f}-{max(lib_abs):.0f}", "", "per Mb", "")
rec("predicted_per_stem_from_pair_rate", f"{mr*min(cal.values())/1e6:.0f}-{mr*max(cal.values())/1e6:.0f}", "mean per-library rate x callable", "calls", "")
rec("stem_private_per_Mb_all", f"{min(rate.values()):.1f}-{max(rate.values()):.1f}", "", "per Mb", "")
# depth and detection
def binom(n, k, p): return math.comb(n, k) * p**k * (1 - p)**(n - k)
def pdet(dp, f):
    t = 0.0
    for k in range(2, dp + 1):
        if k / dp < 0.10: continue
        t += binom(dp, k, f) * sum(binom(k, j, 0.5) for j in range(k + 1) if max(j, k - j) / k <= 0.95)
    return t
cache = {}
hist = defaultdict(lambda: defaultdict(int))
for c in CHR:
    for r in tsv(f"{C}/{c}_depth.tsv"): hist[r["stem"]][int(r["depth_capped200"])] += int(r["n_sites"])
med = {}; P6 = {}; P12 = {}
for s in stems:
    h = hist[s]; n = sum(h.values()); cum = 0
    for d in sorted(h):
        cum += h[d]
        if cum >= n / 2: med[s] = d; break
    for f, D in ((1/6, P6), (1/12, P12)):
        D[s] = sum(h[d] * cache.setdefault((d, f), pdet(min(d, 200), f)) for d in h) / n
rec("median_depth_per_stem", f"{min(med.values())}-{max(med.values())}", "callable sites, chr1-11", "reads", "")
pool = lambda D: sum(D[s] * cal[s] for s in stems) / sum(cal.values())
rec("detect_f1_6", f"{pool(P6):.2f} ({min(P6.values()):.2f}-{max(P6.values()):.2f} per stem)", "all carrier filters, real depths", "probability", "")
rec("detect_f1_12", f"{pool(P12):.2f} ({min(P12.values()):.2f}-{max(P12.values()):.2f} per stem)", "", "probability", "")
rec("detect_DP10", f"{pdet(10,1/6):.2f} / {pdet(10,1/12):.2f}", "f=1/6 / 1/12 at DP 10", "probability", "")
L = sum(cal.values()) / len(cal); rec("mean_callable_Mb", round(L / 1e6, 1), "", "Mb", "")
for mu in (7.4e-11, 1e-10, 1e-9, 6.9e-9, 3.4e-8):
    rec(f"expected_mu{mu:.1e}", " / ".join(f"{mu*t*6*L:.1f}" for t in (78, 177)), "78 / 177 yr, x6 copies, all detected", "mutations", "assumed rates")
with open(f"{OUT}/mr03_numbers.tsv", "w") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n"); w.writerow(["quantity", "value", "scope", "unit", "evidence", "source"]); w.writerows(rows)
