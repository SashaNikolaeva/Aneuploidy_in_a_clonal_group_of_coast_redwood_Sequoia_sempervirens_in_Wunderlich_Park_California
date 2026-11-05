#!/usr/bin/env bash
# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: One pass over the survey VCF: per-SNP read depth and reference/alternative read counts at
# SNPs placed on giant sequoia chromosomes. Input of the allele-balance, allele-dosage and allele-presence
# analyses.
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
# Step 4 (decision 2.4): ONE VCF pass, all 284 samples, per-SNP FORMAT/DP and RO,AO at the SEGI-mapped survey SNPs
# (d01/d03 commands without -S). Reused by R1a/R1b (untruncated VCF-DP column), R1c, R2-364 and R3.
# Output (scratch): /space/s1/sashanikolaeva/clone_finalize_work/run/vcf_dp_ad.tsv.gz
set -euo pipefail
V=/space/s1/santang3/prg/results/variants/vcf/allScaffs_allRegions_filtered.vcf.gz
M=/space/s1/sashanikolaeva/Results/clone/snp_segi_alignment/snp_to_segi_chr.tsv
W=/space/s1/sashanikolaeva/clone_finalize_work/run
BCF=/usr/local/bin/bcftools
tail -n +2 $M | cut -f1,2 | sort -k1,1 -k2,2n > $W/targets.tsv
echo "targets: $(wc -l < $W/targets.tsv); samples in VCF: $($BCF query -l $V | wc -l)"
$BCF query -T $W/targets.tsv -H -f '%CHROM\t%POS[\t%DP][\t%RO,%AO]\n' $V | gzip -1 > $W/vcf_dp_ad.tsv.gz
echo "rows (incl. header): $(zcat $W/vcf_dp_ad.tsv.gz | wc -l)"
sha256sum $W/vcf_dp_ad.tsv.gz
echo DONE
