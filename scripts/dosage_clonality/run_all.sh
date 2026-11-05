#!/usr/bin/env bash
# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Driver of the allele-dosage and clonality analyses.
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
# Dosage + clonality run (PREREG_v2). Usage: bash run_all.sh [output_dir]  (default analysis/dosage_clonality_run)
set -euo pipefail
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLCONFIGDIR=/space/s1/sashanikolaeva/clone_finalize_work/run/mpl
O=${1:-../../analysis/dosage_clonality_run}
mkdir -p "$O"; export DC_RUN=$(cd "$O" && pwd)
PY="$HOME/miniconda3/envs/quast/bin/python3 -I -B"
for s in 01_inputs_sites 02_A9_dosage 04_B1_clonality 05_B2_discordance 03_lattice 06_summary; do
  [ -f "$s.py" ] || continue
  /usr/bin/time -f "$s %U s user %e s wall" $PY $s.py > "$DC_RUN/$s.log" 2>&1
  tail -1 "$DC_RUN/$s.log"
done
