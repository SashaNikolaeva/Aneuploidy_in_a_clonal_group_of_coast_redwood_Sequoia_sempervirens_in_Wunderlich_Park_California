#!/bin/bash
# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Shows why the re-analysis on the coast redwood assembly returned no sites: the scaffold
# names in the region file do not match those in the alignments (Supplementary Note S1).
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
# Reviewer check: full BED vs BAM contig-name overlap for the Neale re-run (read-only).
BED=/space/s1/santang3/prg/resources/data/Annotation/SESE.2.1.CDS.merge.bed
BAM=/space/s1/sashanikolaeva/Results/clone/somatic_mutations_neale/bams/RedwoodCap_13_WP1_neale_sorted_dupsMarked.bam
echo "BED lines: $(wc -l < $BED)"
echo "BED total bp: $(awk '{s+=$3-$2} END{print s}' $BED)"
cut -f1 $BED | sort -u > ${TMPDIR:-/tmp}/bed_names.txt
samtools view -H $BAM | awk -F'\t' '$1=="@SQ"{sub("SN:","",$2); print $2}' | sort -u > ${TMPDIR:-/tmp}/bam_names.txt
echo "BED contigs: $(wc -l < ${TMPDIR:-/tmp}/bed_names.txt)"
echo "BAM contigs: $(wc -l < ${TMPDIR:-/tmp}/bam_names.txt)"
echo "exact overlap: $(comm -12 ${TMPDIR:-/tmp}/bed_names.txt ${TMPDIR:-/tmp}/bam_names.txt | wc -l)"
echo "overlap after stripping ;HRSCAF=: $(sed 's/;HRSCAF=.*//' ${TMPDIR:-/tmp}/bam_names.txt | sort -u | comm -12 ${TMPDIR:-/tmp}/bed_names.txt - | wc -l)"
echo "BED line 1: $(head -1 $BED)"
echo "mpileup cmd in Neale script:"; grep -n -A8 '"samtools", "mpileup"' /space/s1/sashanikolaeva/Scripts/GitHub/redwood-clone-paper/scripts/somatic_mutations_neale.py
echo "BAM @PG (last 3):"; samtools view -H $BAM | grep '^@PG' | tail -3 | cut -c1-300
