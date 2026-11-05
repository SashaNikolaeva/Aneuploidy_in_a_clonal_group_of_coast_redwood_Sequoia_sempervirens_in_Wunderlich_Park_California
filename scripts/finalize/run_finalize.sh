#!/usr/bin/env bash
# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Driver of the first, pre-registered run (223 genets); kept as a record. The paper uses
# run_finalize_224.sh.
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
# Clone-paper finalize run R0-R5 (decision analysis/2026-10-26_finalize_analyses_decision.md, sha256 f24c1d5b...2b34).
# Step 0 (00_prereg.py) was run once before any data were opened and refuses to re-run.
# Outputs: analysis/finalize_run/; scratch: /space/s1/sashanikolaeva/clone_finalize_work/run/. Total < 3 CPU-min.
set -euo pipefail
cd /space/s1/sashanikolaeva/Scripts/GitHub/redwood-clone-paper
PY="$HOME/miniconda3/envs/quast/bin/python3 -I -B"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLCONFIGDIR=/space/s1/sashanikolaeva/clone_finalize_work/run/mpl
O=analysis/finalize_run
$PY scripts/finalize/01_inputs_R0g.py      > $O/01_inputs_R0g.log 2>&1      # input hashes, sample sets, R0g gate (stops on failure)
bash scripts/finalize/03_vcf_pass.sh        > $O/03_vcf_pass.log 2>&1        # one VCF pass (DP, RO/AO), all samples
$PY scripts/finalize/vcfcache.py            > $O/03b_vcfcache.log 2>&1
$PY scripts/finalize/02_R1ab.py             > $O/02_R1ab.log 2>&1            # R1a, R1b
$PY scripts/finalize/04_R2_screen.py        > $O/04_R2_screen.log 2>&1       # R2, R2-WPR, V5, G3 mapping, R1a positive controls
$PY scripts/finalize/05_R1c_positional.py   > $O/05_R1c_positional.log 2>&1  # R1c-g, R1c, R2-364
$PY scripts/finalize/06_R3_allelic.py       > $O/06_R3_allelic.log 2>&1      # R3
$PY scripts/finalize/07_R4_G000.py          > $O/07_R4_G000.log 2>&1         # R4
$PY scripts/finalize/08_R5_clonality.py     > $O/08_R5_clonality.log 2>&1    # R5 (gate R5-g; exits 3 on mismatch)
$PY scripts/finalize/09_R0_fig6_table1.py   > $O/09_R0_fig6_table1.log 2>&1  # R0 Fig 6, Table 1 additions, profile figure
