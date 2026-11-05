#!/bin/bash
# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Checks that the stem search and the replicate calibration used the same BAM files.
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
# Reviewer check: were the 9-stem scan and the replicate calibration run on the same BAM build? (read-only)
R=/space/s1/sashanikolaeva/Results/clone/somatic_mutations_mpileup
B=/space/s1/sashanikolaeva/Results/clone/nquire/bam
echo "== output timestamps"; ls -l --time-style=full-iso $R/run.log $R/per_ramet_counts.tsv $R/private_variants.tsv $R/fp_calibration_bam.tsv | awk '{print $6, $7, $9}'
echo "== BAM timestamps (13 libraries)"; ls -l --time-style=full-iso $B/*.no_chlor.sorted.bam | awk '{print $6, $7, $9}'
echo "== alignment log start"; grep -m1 "Alignment pipeline" $B/align_somatic_samples.log
echo "== current BAM @PG includes markdup?"; for s in 13_WP1 31_c3 22_MRD; do echo "$s: $(samtools view -H $B/$s.no_chlor.sorted.bam | grep -c 'samtools markdup')"; done
echo "== chr1 private calls: March scan (per_ramet_counts.tsv) vs exact re-run on current BAMs (rv_calib_chr1.tsv)"
paste <(awk -F'\t' '$2=="chr1"{print $1"\t"$3"\t"$4}' $R/per_ramet_counts.tsv | sort) \
      <(join <(awk -F'\t' '$1~/^private_/ && $1!="private_total" && $1!~/with/{sub("private_","",$1); print $1"\t"$2}' rv_calib_chr1.tsv | sort) \
             <(awk -F'\t' '$1~/^callable_/{sub("callable_","",$1); print $1"\t"$2}' rv_calib_chr1.tsv | sort) -t $'\t') | \
  awk 'BEGIN{print "stem\tmarch_calls\tmarch_callable\tnow_stem\tnow_calls\tnow_callable"} {print}'
echo "== calibration reproduces on current BAMs? fp_calibration_bam.log MRD progress near chr1 callable size:"; grep -m1 "1,945,187 callable" $R/fp_calibration_bam.log
echo "   rv re-run chr1 MRD: $(awk -F'\t' '$1=="pair_private_MRD_22_MRD"||$1=="pair_private_MRD_28_MRD"{s+=$2} END{print s}' rv_calib_chr1.tsv) calls over $(awk -F'\t' '$1=="pair_callable_MRD"{print $2}' rv_calib_chr1.tsv) callable sites"
