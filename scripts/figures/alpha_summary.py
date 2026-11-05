# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Prints the alpha numbers of the Abstract, Results (chromosome-4 and other-chromosome alpha
# of the stems, survey minimum, WPR and C1) and Supplementary Table S1 from the saved output of
# scripts/finalize/09_R0_fig6_table1.py.
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
"""Numbers of the copy-number Results, Abstract and Supplementary Table S1 (alpha_ij, 2026-11-05 version of the paper),
read from the saved outputs of finalize/09_R0_fig6_table1.py; nothing is recomputed.

alpha_ij = share of sample i's reads on chromosome j / median share of the 224 survey trees (fig6_source.tsv);
expected one-copy-loss level from fig6_expected_levels.tsv; PCR pools and groups of the stems from R0_table1_additions.tsv.
fig6_source.tsv keeps 4 decimals, so values ending in 5 were rounded from a full-precision recompute of the same counts
(2026-11-05): 13_WP1 0.84849 -> 0.848, 22_WPd2 0.84652 -> 0.847, WPR 0.98446 -> 0.984, stem mean 0.85117 -> 0.851.
Output: analysis/finalize_run_224/alpha_summary.tsv (quantity, value) and stdout.
Run: ~/miniconda3/envs/quast/bin/python3 -I scripts/figures/alpha_summary.py
"""
import pandas as pd

RUN = "/space/s1/sashanikolaeva/Scripts/GitHub/redwood-clone-paper/analysis/finalize_run_224"
WPR, C1 = "22_WPR", "16_C1"
CH = [f"chr{i}" for i in range(1, 12)]

src = pd.read_csv(f"{RUN}/fig6_source.tsv", sep="\t", index_col=0)
lv = pd.read_csv(f"{RUN}/fig6_expected_levels.tsv", sep="\t", index_col=0)
t1 = pd.read_csv(f"{RUN}/R0_table1_additions.tsv", sep="\t", index_col=0)
a = src[[f"alpha_{c}" for c in CH]]
a.columns = CH
stems = src.index[src.group == "clone_ramet"]
survey = src.index[src.group == "PANEL_ALL"]
other = [c for c in CH if c != "chr4"]

rows = [("n stems", len(stems)), ("n survey trees", len(survey)),
        ("stems chr4 alpha min", a.loc[stems, "chr4"].min()), ("stems chr4 alpha max", a.loc[stems, "chr4"].max()),
        ("stems chr4 alpha mean", a.loc[stems, "chr4"].mean()),
        ("expected chr4 alpha, one copy lost", lv.loc["chr4", "alpha5_one_copy_loss"]),
        ("survey chr4 alpha min", a.loc[survey, "chr4"].min()),
        ("stems other chromosomes alpha min", a.loc[stems, other].min().min()),
        ("stems other chromosomes alpha max", a.loc[stems, other].max().max()),
        ("WPR chr4 alpha", a.loc[WPR, "chr4"]), ("C1 chr4 alpha", a.loc[C1, "chr4"]),
        ("stems PCR pools", t1.loc[stems, "pcr_pool"].nunique()), ("stems PCR groups", t1.loc[stems, "pcr_group"].nunique())]
rows += [(f"{s} chr4 alpha", a.loc[s, "chr4"]) for s in stems]
out = pd.DataFrame(rows, columns=["quantity", "value"])
out.to_csv(f"{RUN}/alpha_summary.tsv", sep="\t", index=False, float_format="%.4f")
print(out.to_string(index=False))
