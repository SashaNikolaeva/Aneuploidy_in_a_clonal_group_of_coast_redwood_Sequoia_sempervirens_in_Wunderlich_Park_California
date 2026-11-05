#!/bin/bash
# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Aligns the nine stems and the two pairs of technical replicates (13 libraries) to giant
# sequoia v2.0 for the stem-private allele search (Supplementary Note S1).
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
# align_somatic_samples.sh
# Align 13 samples (9 clone ramets + 4 technical replicates) to SEGI v2.0,
# mark duplicates, filter to chr1-chr11, index.
# Output: {sample}.no_chlor.sorted.bam in BAM_DIR

set -euo pipefail

# ── Paths ─────────────────────────────────────────────────────────────────────
BWA="/home/sashanikolaeva/miniconda3/pkgs/bwa-0.7.19-h577a1d6_1/bin/bwa"
SAMTOOLS="/usr/local/bin/samtools"
REF="/space/s1/sashanikolaeva/Data/SEGI_genome/SEGI.2.0.fa"
FASTQ_DIR="/space/s1/sashanikolaeva/Results/clone/nquire/bam/SEGI"
BAM_DIR="/space/s1/sashanikolaeva/Results/clone/nquire/bam"
TMP_DIR="${BAM_DIR}/tmp_align"
LOG="${BAM_DIR}/align_somatic_samples.log"
THREADS=20

# Nuclear regions to keep (exclude chloroplast + unplaced scaffolds)
KEEP_REGIONS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11"

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
)

mkdir -p "${TMP_DIR}"
echo "=== Alignment pipeline — $(date) ===" | tee -a "${LOG}"

for SAMPLE in "${SAMPLES[@]}"; do
    OUT="${BAM_DIR}/${SAMPLE}.no_chlor.sorted.bam"

    if [[ -f "${OUT}" ]]; then
        echo "[SKIP] ${SAMPLE} — already exists" | tee -a "${LOG}"
        continue
    fi

    R1="${FASTQ_DIR}/RedwoodCap_${SAMPLE}_R1_concatenated.fastq.gz"
    R2="${FASTQ_DIR}/RedwoodCap_${SAMPLE}_R2_concatenated.fastq.gz"

    if [[ ! -f "${R1}" || ! -f "${R2}" ]]; then
        echo "[WARN] ${SAMPLE} — FASTQ not found, skipping" | tee -a "${LOG}"
        continue
    fi

    echo "" | tee -a "${LOG}"
    echo "[$(date +%H:%M:%S)] ${SAMPLE} — starting..." | tee -a "${LOG}"

    TMP_NAMESORTED="${TMP_DIR}/${SAMPLE}.namesorted.bam"
    TMP_FIXMATE="${TMP_DIR}/${SAMPLE}.fixmate.bam"
    TMP_SORTED="${TMP_DIR}/${SAMPLE}.sorted.bam"
    TMP_DEDUP="${TMP_DIR}/${SAMPLE}.dedup.bam"
    METRICS="${TMP_DIR}/${SAMPLE}.dup_metrics.txt"

    # Step 1: Align + name sort (required for fixmate)
    ${BWA} mem -t ${THREADS} \
        -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" \
        "${REF}" "${R1}" "${R2}" 2>>"${LOG}" | \
        ${SAMTOOLS} sort -n -@ 4 -T "${TMP_DIR}/${SAMPLE}_namesort_tmp" \
            -o "${TMP_NAMESORTED}" -
    echo "  [$(date +%H:%M:%S)] Aligned and name-sorted" | tee -a "${LOG}"

    # Step 2: Fixmate (adds ms mate-score tag required by markdup)
    ${SAMTOOLS} fixmate -m -@ ${THREADS} "${TMP_NAMESORTED}" "${TMP_FIXMATE}"
    rm -f "${TMP_NAMESORTED}"
    echo "  [$(date +%H:%M:%S)] Fixmate done" | tee -a "${LOG}"

    # Step 3: Coordinate sort
    ${SAMTOOLS} sort -@ 4 -T "${TMP_DIR}/${SAMPLE}_sort_tmp" \
        -o "${TMP_SORTED}" "${TMP_FIXMATE}"
    ${SAMTOOLS} index -c "${TMP_SORTED}"
    rm -f "${TMP_FIXMATE}"
    echo "  [$(date +%H:%M:%S)] Coordinate sorted" | tee -a "${LOG}"

    # Step 4: Mark duplicates
    ${SAMTOOLS} markdup -@ ${THREADS} -T "${TMP_DIR}/${SAMPLE}_md_tmp" \
        -f "${METRICS}" "${TMP_SORTED}" "${TMP_DEDUP}"
    ${SAMTOOLS} index -c "${TMP_DEDUP}"
    echo "  [$(date +%H:%M:%S)] Duplicates marked" | tee -a "${LOG}"

    # Step 5: Filter to chr1-chr11 only
    ${SAMTOOLS} view -h -@ ${THREADS} "${TMP_DEDUP}" ${KEEP_REGIONS} | \
        ${SAMTOOLS} sort -@ 4 -T "${TMP_DIR}/${SAMPLE}_filter_tmp" \
            -o "${OUT}" -
    ${SAMTOOLS} index -c "${OUT}"
    echo "  [$(date +%H:%M:%S)] Filtered to chr1-chr11" | tee -a "${LOG}"

    # Stats
    TOTAL=$(${SAMTOOLS} view -c "${TMP_DEDUP}")
    KEPT=$(${SAMTOOLS} view -c "${OUT}")
    PCT_KEPT=$(awk "BEGIN {printf \"%.1f\", 100.0 * ${KEPT} / ${TOTAL}}")
    echo "  Reads: ${TOTAL} total → ${KEPT} nuclear (${PCT_KEPT}%)" | tee -a "${LOG}"

    # Cleanup temp files
    rm -f "${TMP_SORTED}" "${TMP_SORTED}.csi" \
          "${TMP_DEDUP}" "${TMP_DEDUP}.csi" "${METRICS}"
    echo "  [$(date +%H:%M:%S)] Done: ${OUT}" | tee -a "${LOG}"
done

echo "" | tee -a "${LOG}"
echo "=== All done — $(date) ===" | tee -a "${LOG}"
