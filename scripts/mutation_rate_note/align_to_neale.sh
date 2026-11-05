#!/bin/bash
# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Aligns the 13 libraries to the coast redwood assembly for the run above.
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
# align_to_neale.sh
# Align 13 samples (9 clone ramets + 4 technical replicates) to the Neale SESE 2.1
# reference genome — the same reference used for the population structure VCF.
# This is required for the somatic mutation calibration (private variant calling)
# to be on the same reference as the VCF-based shared variant scan.
#
# Mirrors the santang3 pipeline exactly:
#   bwa mem | samtools view | collate | fixmate | sort | markdup
# Run in tmux: tmux new-session -s clone_neale_align

set -euo pipefail

# ── Paths ─────────────────────────────────────────────────────────────────────
BWA="/home/sashanikolaeva/miniconda3/envs/ntsynt-310/bin/bwa"
REF="/space/s1/santang3/prg/results/ref/neale/neale.fa"
BWA_PREFIX="/space/s1/sashanikolaeva/Results/clone/somatic_mutations_neale/ref/neale"
FASTQ_DIR="/space/s1/sashanikolaeva/Results/clone/nquire/bam/SEGI"
BAM_DIR="/space/s1/sashanikolaeva/Results/clone/somatic_mutations_neale/bams"
TMP_DIR="${BAM_DIR}/tmp"
LOG="${BAM_DIR}/align_to_neale.log"
THREADS=20

# ── Samples ───────────────────────────────────────────────────────────────────
SAMPLES=(
    13_WP1
    6_WP5
    31_c3
    31_WP2_1
    28_WP7
    28_WP6
    14_WP3
    22_WPd2
    33_c6
    22_MRD
    28_MRD
    16_416-1
    35_416_1
    22_WPR
    16_C1
)

mkdir -p "${BAM_DIR}" "${TMP_DIR}" "$(dirname ${BWA_PREFIX})"
echo "=== Neale alignment pipeline — $(date) ===" | tee -a "${LOG}"

# ── Step 1: BWA index (skip if already done) ──────────────────────────────────
if [[ ! -f "${BWA_PREFIX}.bwt" ]]; then
    echo "[$(date +%H:%M:%S)] Building BWA index for Neale reference (~6-8 hrs, ~45 GB)" | tee -a "${LOG}"
    echo "[$(date +%H:%M:%S)] Index prefix: ${BWA_PREFIX}" | tee -a "${LOG}"
    "${BWA}" index -p "${BWA_PREFIX}" "${REF}" 2>> "${LOG}"
    echo "[$(date +%H:%M:%S)] BWA index complete" | tee -a "${LOG}"
else
    echo "[$(date +%H:%M:%S)] BWA index already exists, skipping" | tee -a "${LOG}"
fi

# ── Step 2: Align each sample ─────────────────────────────────────────────────
for SAMPLE in "${SAMPLES[@]}"; do
    OUT="${BAM_DIR}/RedwoodCap_${SAMPLE}_neale_sorted_dupsMarked.bam"

    if [[ -f "${OUT}" && -f "${OUT}.csi" ]]; then
        echo "[SKIP] ${SAMPLE} — BAM already exists" | tee -a "${LOG}"
        continue
    fi

    R1="${FASTQ_DIR}/RedwoodCap_${SAMPLE}_R1_concatenated.fastq.gz"
    R2="${FASTQ_DIR}/RedwoodCap_${SAMPLE}_R2_concatenated.fastq.gz"

    if [[ ! -f "${R1}" || ! -f "${R2}" ]]; then
        echo "[WARN] ${SAMPLE} — FASTQ not found, skipping" | tee -a "${LOG}"
        continue
    fi

    RG="@RG\tID:RedwoodCap_${SAMPLE}\tCN:MiSeq\tPL:ILLUMINA\tPM:MiSeq\tSM:RedwoodCap_${SAMPLE}"

    echo "" | tee -a "${LOG}"
    echo "[$(date +%H:%M:%S)] ${SAMPLE} — aligning..." | tee -a "${LOG}"

    "${BWA}" mem -t "${THREADS}" -R "${RG}" "${BWA_PREFIX}" "${R1}" "${R2}" 2>>"${LOG}" \
        | samtools view -hb -o - - \
        | samtools collate --threads "${THREADS}" -o - - "${TMP_DIR}/${SAMPLE}_collate" \
        | samtools fixmate --threads "${THREADS}" -m - - \
        | samtools sort --threads "${THREADS}" -T "${TMP_DIR}/${SAMPLE}_sort" -o - \
        | samtools markdup --threads "${THREADS}" -T "${TMP_DIR}" \
            -f "${BAM_DIR}/RedwoodCap_${SAMPLE}_dupStats.txt" - "${OUT}"

    samtools index -c "${OUT}"
    echo "[$(date +%H:%M:%S)] ${SAMPLE} — done" | tee -a "${LOG}"
done

echo "" | tee -a "${LOG}"
echo "=== All samples done — $(date) ===" | tee -a "${LOG}"
echo "NEALE_ALIGN_EXIT:0" | tee -a "${LOG}"
