#!/usr/bin/env bash
# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Driver of the analysis run whose numbers are printed in the paper (the 224 trees of the
# population-structure analysis); gives the order of the steps.
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
# Driver for the run whose numbers are printed in the paper: FINALIZE_SET=popgen224 (the popgen paper's 224 trees;
# Sasha 2026-10-26, DV14 in analysis/finalize_run_224/DEVIATIONS_224.md). The 2026-10-26 run was started step by step;
# this driver records the order (written 2026-10-29). run_finalize.sh is the frozen pre-registered run (analysis/finalize_run/).
# Outputs: analysis/finalize_run_224/ and clone_paper_figures/{clonality,aneuploidy}/; scratch: clone_finalize_work/run/.
# Step 00_prereg.py was run once before any data were opened and is not part of this driver.
set -euo pipefail
cd /space/s1/sashanikolaeva/Scripts/GitHub/redwood-clone-paper
PY="$HOME/miniconda3/envs/quast/bin/python3 -I -B"
export FINALIZE_SET=popgen224
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLCONFIGDIR=/space/s1/sashanikolaeva/clone_finalize_work/run/mpl
export SOURCE_DATE_EPOCH=${SOURCE_DATE_EPOCH:-1793232000}   # fixed PDF metadata date, so re-runs give identical PDFs
O=analysis/finalize_run_224
$PY scripts/figures/fig2_clonality.py       > $O/fig2_clonality.log 2>&1     # Fig 2 and its pair tables (hashed by step 01, read by R5)
$PY scripts/finalize/01_inputs_R0g.py      > $O/01_inputs_R0g.log 2>&1      # input hashes, sample sets, R0g gate (stops on failure)
bash scripts/finalize/03_vcf_pass.sh        > $O/03_vcf_pass.log 2>&1        # one VCF pass (DP, RO/AO), all samples
$PY scripts/finalize/vcfcache.py            > $O/03b_vcfcache.log 2>&1
$PY scripts/finalize/02_R1ab.py             > $O/02_R1ab.log 2>&1            # R1a, R1b (Table 2)
$PY scripts/finalize/04_R2_screen.py        > $O/04_R2_screen.log 2>&1       # R2 screen, gain controls
$PY scripts/finalize/05_R1c_positional.py   > $O/05_R1c_positional.log 2>&1  # R1c positional scan, R2-364
$PY scripts/finalize/06_R3_allelic.py       > $O/06_R3_allelic.log 2>&1      # R3 allele balance
$PY scripts/finalize/07_R4_G000.py          > $O/07_R4_G000.log 2>&1         # R4 stem-private alleles
$PY scripts/finalize/08_R5_clonality.py     > $O/08_R5_clonality.log 2>&1    # R5 (gate R5-g; exits 3 on mismatch)
$PY scripts/finalize/09_R0_fig6_table1.py   > $O/09_R0_fig6_table1.log 2>&1  # Fig 3, Fig S1, Table 1 additions
