#!/usr/bin/env bash
# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Driver of the two flow-cytometry runs analysed with the pipeline.
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
# The two flow runs analysed on 2026-10-26 (R6), as worked examples of scripts/flow/README.md.
# Sample sheets and PREREGs already exist (make_sample_sheet.py / freeze_prereg.py refuse to overwrite).
set -euo pipefail
cd /space/s1/sashanikolaeva/Scripts/GitHub/redwood-clone-paper
PY="$HOME/miniconda3/envs/quast/bin/python3 -I -B"
export MPLCONFIGDIR=/space/s1/sashanikolaeva/clone_finalize_work/run/mpl OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
D=/space/s1/sashanikolaeva/flow_cytometry_data_clone
OUT=analysis/finalize_run/flow
FIG=/space/s1/sashanikolaeva/Scripts/GitHub/clone_paper_figures/flow
$PY scripts/flow/run_flow.py --data $D --sheet scripts/flow/sample_sheets/2026-10_october.tsv --run 2026-10_october \
    --out $OUT --figdir $FIG --scan "Clone, Oct14/**/*.fcs" --scan "Experiment_002/**/*.fcs"
$PY scripts/flow/run_flow.py --data $D --sheet scripts/flow/sample_sheets/2025-04_berkeley.tsv --run 2025-04_berkeley \
    --out $OUT --figdir $FIG --scan "*.fcs" --scan "SashaNikolaeva/**/*.fcs"
# legacy 2025 figures (voltage labels now read from each tube's header); they import fcsparser -> kmer_analysis env
for s in flowcyto_clone_saturation flowcyto_clone_scatter; do
  "$HOME/miniconda3/envs/kmer_analysis/bin/python" -I -B scripts/$s.py
done
