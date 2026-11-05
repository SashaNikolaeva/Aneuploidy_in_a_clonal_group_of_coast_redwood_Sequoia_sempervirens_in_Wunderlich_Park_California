#!/usr/bin/env python3
# ----------------------------------------------------------------------------------------------------
# Paper: Aneuploidy in a clonal group of coast redwood (Sequoia sempervirens) in Wunderlich Park, California.
# Role in the paper: Places the survey SNPs on giant sequoia chromosomes by aligning each SNP's flanking
# sequence (minimap2). Used by the allele-balance, allele-dosage and allele-presence analyses.
# Annotation added 2026-11-05 for this repository; the code below is unchanged from the analysis run.
# Absolute paths refer to the analysis server; see README.md.
# ----------------------------------------------------------------------------------------------------
"""
Extract 600bp regions around all SNPs and align to SEGI genome.

Steps:
1. Parse VCF to get all SNP positions
2. Extract 600bp flanking regions from redwood reference using samtools faidx
3. Align to SEGI genome using minimap2
4. Parse alignments to get SNP -> SEGI chromosome mapping
5. Output list of SNPs mapping to each chromosome (especially chr4)
"""
import subprocess
import sys
import gzip
import os
from collections import defaultdict

# File paths
VCF_FILE = '/space/s1/santang3/prg/results/variants/vcf/allScaffs_allRegions_filtered.vcf.gz'
REDWOOD_REF = '/space/s1/santang3/prg/results/ref/neale/neale.fa'
SEGI_REF = '/space/s1/sashanikolaeva/Data/SEGI_genome/SEGI.2.0.fa'
OUTPUT_DIR = '/space/s1/sashanikolaeva/Results/clone/snp_segi_alignment'
FLANK_SIZE = 300  # 300bp on each side = 600bp total

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    snp_regions_fa = os.path.join(OUTPUT_DIR, 'snp_regions.fa')
    alignment_paf = os.path.join(OUTPUT_DIR, 'snp_regions_to_segi.paf')
    output_mapping = os.path.join(OUTPUT_DIR, 'snp_to_segi_chr.tsv')
    chr4_snps = os.path.join(OUTPUT_DIR, 'chr4_snps.txt')

    # Step 1: Extract SNP positions and create BED-like regions
    print("Step 1: Parsing VCF to get SNP positions...", file=sys.stderr)
    snp_positions = []

    with gzip.open(VCF_FILE, 'rt') as vcf:
        for line in vcf:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            scaffold = fields[0]
            position = int(fields[1])
            snp_positions.append((scaffold, position))

    print(f"Found {len(snp_positions)} SNPs", file=sys.stderr)

    # Step 2: Create regions file and extract with samtools faidx
    print("Step 2: Extracting 600bp flanking regions...", file=sys.stderr)

    # Write regions to a temp file for samtools
    regions_file = os.path.join(OUTPUT_DIR, 'regions.txt')
    with open(regions_file, 'w') as f:
        for scaffold, pos in snp_positions:
            start = max(1, pos - FLANK_SIZE)
            end = pos + FLANK_SIZE
            # Region format: scaffold:start-end
            # Name will encode original SNP info: scaffold_position
            region = f"{scaffold}:{start}-{end}"
            f.write(f"{region}\n")

    # Extract regions using samtools faidx
    print("Running samtools faidx (this may take a while)...", file=sys.stderr)
    cmd = f"samtools faidx {REDWOOD_REF} -r {regions_file} > {snp_regions_fa}"
    subprocess.run(cmd, shell=True, check=True)

    # Rename FASTA headers to include SNP position info
    print("Renaming FASTA headers to include SNP info...", file=sys.stderr)
    renamed_fa = os.path.join(OUTPUT_DIR, 'snp_regions_renamed.fa')
    with open(snp_regions_fa, 'r') as fin, open(renamed_fa, 'w') as fout:
        snp_idx = 0
        for line in fin:
            if line.startswith('>'):
                # Original header: >Scaffold_X;HRSCAF=Y:start-end
                # We need to associate with original SNP position
                scaffold, pos = snp_positions[snp_idx]
                # Create new header with SNP ID
                new_header = f">{scaffold}_{pos}\n"
                fout.write(new_header)
                snp_idx += 1
            else:
                fout.write(line)

    # Step 3: Align to SEGI with minimap2
    print("Step 3: Aligning to SEGI genome with minimap2...", file=sys.stderr)
    cmd = f"minimap2 -x asm20 -t 20 {SEGI_REF} {renamed_fa} > {alignment_paf}"
    subprocess.run(cmd, shell=True, check=True)

    # Step 4: Parse PAF and map SNPs to chromosomes
    print("Step 4: Parsing alignment results...", file=sys.stderr)
    snp_to_chr = {}
    chr_counts = defaultdict(int)

    with open(alignment_paf, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            query_name = fields[0]  # scaffold_position
            target_chr = fields[5]  # SEGI chromosome
            mapq = int(fields[11]) if len(fields) > 11 else 0

            # Only keep primary/best alignments (mapq >= 30)
            if mapq >= 30:
                # Extract original scaffold and position from query name
                # Format: Scaffold_X;HRSCAF=Y_position
                parts = query_name.rsplit('_', 1)
                if len(parts) == 2:
                    scaffold = parts[0]
                    try:
                        position = int(parts[1])
                        snp_id = f"{scaffold}_{position}"
                        if snp_id not in snp_to_chr:
                            snp_to_chr[snp_id] = target_chr
                            chr_counts[target_chr] += 1
                    except ValueError:
                        continue

    # Step 5: Write output
    print("Step 5: Writing results...", file=sys.stderr)

    with open(output_mapping, 'w') as f:
        f.write("scaffold\tposition\tsegi_chr\n")
        for snp_id, segi_chr in sorted(snp_to_chr.items()):
            parts = snp_id.rsplit('_', 1)
            scaffold = parts[0]
            position = parts[1]
            f.write(f"{scaffold}\t{position}\t{segi_chr}\n")

    # Write chr4 SNPs separately
    with open(chr4_snps, 'w') as f:
        f.write("scaffold\tposition\n")
        for snp_id, segi_chr in sorted(snp_to_chr.items()):
            if segi_chr == 'chr4':
                parts = snp_id.rsplit('_', 1)
                scaffold = parts[0]
                position = parts[1]
                f.write(f"{scaffold}\t{position}\n")

    # Summary
    print(f"\n{'='*50}", file=sys.stderr)
    print("SUMMARY", file=sys.stderr)
    print(f"{'='*50}", file=sys.stderr)
    print(f"Total SNPs in VCF: {len(snp_positions)}", file=sys.stderr)
    print(f"SNPs with SEGI alignment (mapq>=30): {len(snp_to_chr)}", file=sys.stderr)
    print(f"\nSNPs per SEGI chromosome:", file=sys.stderr)
    for chr_name in sorted(chr_counts.keys()):
        print(f"  {chr_name}: {chr_counts[chr_name]}", file=sys.stderr)

    chr4_count = chr_counts.get('chr4', 0)
    print(f"\n*** SNPs mapping to chr4: {chr4_count} ***", file=sys.stderr)
    print(f"\nOutput files:", file=sys.stderr)
    print(f"  All mappings: {output_mapping}", file=sys.stderr)
    print(f"  Chr4 SNPs: {chr4_snps}", file=sys.stderr)

if __name__ == '__main__':
    main()
