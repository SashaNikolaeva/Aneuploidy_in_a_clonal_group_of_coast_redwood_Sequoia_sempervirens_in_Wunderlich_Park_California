# Aneuploidy in a clonal group of coast redwood (*Sequoia sempervirens*) in Wunderlich Park, California

Analysis scripts for the paper of the same name (Nikolaeva, Santangelo, Smith, Moore, Stern, Dodd and Nielsen; in preparation).

Nine stems of a clonal group in Wunderlich Park are one genet, and every stem carries five of the six copies of chromosome 4. These scripts produced the numbers, tables and figures of the paper from exome-capture reads of the stems and of 224 coast redwood trees from across the range.

## What is here

- `scripts/`: every script behind a result in the paper, in the folder layout of the analysis (several scripts import each other by path). Each script starts with a header giving its role in the paper; the code below the header is unchanged from the run that produced the results.
- `prereg/`: the tables of tests, thresholds and conclusion wording that were frozen before the analyses were run (Methods), with their SHA-256 checksums: `finalize_PREREG.tsv` (main run), `dosage_clonality_PREREG_v2.tsv` (allele-dosage and batch-aware clonality tests), `flow/` (flow-cytometry runs).

## How the scripts fit together

1. Reads were aligned to the giant sequoia genome v2.0 (`scripts/upstream/alignment_rule4.py`) and counted per chromosome (`scripts/upstream/idxstats_no_chlor_filtered_august.py`).
2. Survey SNPs (from the variant calls of the population-structure analysis) were placed on giant sequoia chromosomes (`scripts/upstream/extract_and_align_snps_to_segi.py`), and identity-by-descent distance matrices were computed (`scripts/upstream/IBD_matrix_args.py`).
3. The main run is `scripts/finalize/run_finalize_224.sh`; the allele-dosage and clonality tests are `scripts/dosage_clonality/run_all.sh`.
4. Figures come from `scripts/figures/`; Figure 3 is plotted from the saved output of `scripts/finalize/09_R0_fig6_table1.py`.

Chromosome copy number is reported as alpha_ij, the share of a sample's reads on chromosome j divided by the median share of the 224 survey trees, following Nikolaeva et al. 2025 (G3 15(5): jkaf063). The frozen main run also computed statistics that the paper no longer reports (a prediction-interval t and bootstrap interval, a positional scan along chromosome 4, and a copy-number screen of the survey); those scripts are included because later steps read their outputs, and their headers say so.

## Scripts

### Upstream: alignments, read counts, SNP placement, distance matrices

| Script | What it does |
|---|---|
| `scripts/upstream/alignment_rule4.py` | Aligns the survey reads to the giant sequoia genome v2.0 (BWA-MEM; samtools markdup -r removes duplicates). Its BAMs are the input of the per-chromosome read counts (Methods: Chromosome copy number). Identified from the output file names; the BAMs have since been deleted. |
| `scripts/upstream/idxstats_no_chlor_filtered_august.py` | Counts primary, properly paired reads with mapping quality >= 60 per giant sequoia chromosome (samtools view -F 256 -f 2 -q 60; samtools idxstats). These counts are n_ij in alpha_ij (Methods: Chromosome copy number; Figure 3; Supplementary Table S1). The input folder in the script was edited after the survey run. |
| `scripts/upstream/extract_and_align_snps_to_segi.py` | Places the survey SNPs on giant sequoia chromosomes by aligning each SNP's flanking sequence (minimap2). Used by the allele-balance, allele-dosage and allele-presence analyses. |
| `scripts/upstream/depth_per_chr_vcf.py` | Read depth per giant sequoia chromosome from the survey VCF. Its output defines the VCF sample set read by scripts/finalize/common.py. |
| `scripts/upstream/het_density_ad_fullpop.py` | Heterozygous-site density per chromosome from VCF allele depths. Read by scripts/finalize/02_R1ab.py as a secondary column; not reported in the paper. |
| `scripts/upstream/IBD_matrix_args.py` | Identity-by-descent distance matrix from the allele-dosage matrix (PCR-corrected and uncorrected). Distances of Figure 2 and of the clonality tests (Methods: Genetic distance matrix). |
| `scripts/upstream/IBD_matrix_exclude_chr4.py` | The same distance matrix without SNPs placed on chromosome 4 (Methods: check that the aneuploid chromosome does not drive the distances). |
| `scripts/upstream/align_somatic_samples.sh` | Aligns the nine stems and the two pairs of technical replicates (13 libraries) to giant sequoia v2.0 for the stem-private allele search (Supplementary Note S1). |
| `scripts/upstream/somatic_mutations_mpileup.py` | Search for alleles carried by one stem only (samtools mpileup on the 13 libraries above; Supplementary Note S1). Re-run on the current BAMs as scripts/mutation_rate_note/mr02_scan_original.py. |

### Main analysis run (copy number, allele balance, stem-private alleles, clonality tables)

| Script | What it does |
|---|---|
| `scripts/finalize/common.py` | Shared definitions of the analysis run: input paths, sample sets (the 224 survey trees, the nine stems, WPR, C1), the read-count statistic and helpers. Imported by every script in scripts/finalize/ and scripts/dosage_clonality/. |
| `scripts/finalize/00_prereg.py` | Freezes the written plan of the analysis run (prereg/finalize_PREREG.tsv) before any data are opened (Methods: the tests and thresholds were written down before the clone statistics were computed). |
| `scripts/finalize/01_inputs_R0g.py` | Records input checksums and the fixed sample sets, and checks that the inputs reproduce the earlier chromosome-4 value before the run continues. |
| `scripts/finalize/02_R1ab.py` | Chromosome-4 copy-number statistic per stem with prediction-interval t, bootstrap CI and the bound on the fraction of chromosome 4 lost. Part of the frozen run; these numbers are NOT reported in the current version of the paper, which reports alpha_ij (scripts/finalize/09_R0_fig6_table1.py). |
| `scripts/finalize/03_vcf_pass.sh` | One pass over the survey VCF: per-SNP read depth and reference/alternative read counts at SNPs placed on giant sequoia chromosomes. Input of the allele-balance, allele-dosage and allele-presence analyses. |
| `scripts/finalize/vcfcache.py` | Converts the VCF extract of 03_vcf_pass.sh into numpy arrays with giant sequoia chromosome and position per SNP. |
| `scripts/finalize/04_R2_screen.py` | Copy-number screen of the 224 survey trees. Not reported in the current version of the paper (the other aneuploid survey trees are cited from Nikolaeva et al. 2025, G3); 09_R0_fig6_table1.py reads its output for the survey range drawn in Figure 3. |
| `scripts/finalize/05_R1c_positional.py` | Read depth along chromosome 4 in 20-Mb bins. Part of the frozen run; NOT reported in the current version of the paper (dropped because it had little power). |
| `scripts/finalize/06_R3_allelic.py` | Allele balance on chromosome 4 relative to the other chromosomes, against depth-matched survey trees (Methods and Results: allele fractions; allele-share ratio in Supplementary Table S1). |
| `scripts/finalize/07_R4_G000.py` | Replication of stem-private alleles (kappa) in the other stems versus same-tree comparison pairs (Methods: Replication of stem-private alleles; Results; Supplementary Note S1). |
| `scripts/finalize/08_R5_clonality.py` | Clonality tables: classes of pairwise distances, within-clone distance versus library depth, and the comparison of distances with and without chromosome 4 (Results: the nine stems are one genet). |
| `scripts/finalize/09_R0_fig6_table1.py` | Computes alpha_ij (read share of a sample on a chromosome divided by the median share of the 224 survey trees) and the expected one-copy levels, saved as fig6_source.tsv and fig6_expected_levels.tsv (Methods and Results: Chromosome copy number; Supplementary Table S1). Also writes the stem locations and PCR pools of Table 1. Its figure was replaced by scripts/figures/fig3_alpha_segi.py. |
| `scripts/finalize/run_finalize_224.sh` | Driver of the analysis run whose numbers are printed in the paper (the 224 trees of the population-structure analysis); gives the order of the steps. |
| `scripts/finalize/run_finalize.sh` | Driver of the first, pre-registered run (223 genets); kept as a record. The paper uses run_finalize_224.sh. |

### Allele-dosage test and batch-aware clonality tests

| Script | What it does |
|---|---|
| `scripts/dosage_clonality/dc_common.py` | Shared inputs of the allele-dosage and batch-aware clonality analyses (sample sets, PCR groups, SNP site metadata). |
| `scripts/dosage_clonality/01_inputs_sites.py` | Input checksums and the SNP sets of the allele-dosage test (1,584 SNPs, 180 on chromosome 4) and of the lattice model. |
| `scripts/dosage_clonality/02_A9_dosage.py` | Allele-dosage test: chromosome 4 versus the other chromosomes, block jackknife, 10,000 PCR-group-matched random sets of nine trees, and its sensitivity (Methods and Results: allele fractions). |
| `scripts/dosage_clonality/03_lattice.py` | Secondary lattice mixture model of allele fractions with its identifiability gate (Results: not identifiable, not applied to chromosome 4). |
| `scripts/dosage_clonality/04_B1_clonality.py` | Batch-aware clonality test: removes PCR-group offsets from all pairwise distances and compares the stems with the bound from same-tree pairs (Methods and Results: clonality). |
| `scripts/dosage_clonality/05_B2_discordance.py` | Allele-presence discordance between libraries after thinning to 30 reads, a clonality measure that does not use the distance matrix (Methods and Results: clonality). |
| `scripts/dosage_clonality/06_summary.py` | Combines the allele-dosage and clonality results into the calls set in advance and draws a summary figure (not in the paper). |
| `scripts/dosage_clonality/run_all.sh` | Driver of the allele-dosage and clonality analyses. |

### Figures

| Script | What it does |
|---|---|
| `scripts/figures/fig1_site.py` | Figure 1 (range map, site map, historical and present imagery, site photos). |
| `scripts/figures/fig2_clonality.py` | Figure 2 (distance heatmap and pairwise distance classes) and the all-pairs comparison with the survey (Results: 24,976 pairs of distinct trees). |
| `scripts/figures/fig3_alpha_segi.py` | Figure 3 (alpha per chromosome for the stems, WPR, C1 and the survey trees), plotted from the saved output of scripts/finalize/09_R0_fig6_table1.py. |

### Stem spacing

| Script | What it does |
|---|---|
| `scripts/exploratory/within_clone_spacing.py` | Within-clone genetic distance versus distance between stems (Spearman rho, 9,999 permutations; Methods and Results). |

### Flow cytometry

| Script | What it does |
|---|---|
| `scripts/flow/fcsio.py` | Minimal FCS file reader used by the flow-cytometry pipeline. |
| `scripts/flow/flowpipe.py` | Flow-cytometry genome-size pipeline with quality grades set in advance (Methods: Flow cytometry; it grades all three tubes as failing). |
| `scripts/flow/freeze_prereg.py` | Freezes the parameters of a flow-cytometry run before the pipeline is run on it. |
| `scripts/flow/make_sample_sheet.py` | Builds the sample sheet of a flow-cytometry run from its FCS files. |
| `scripts/flow/run_flow.py` | Runs the flow-cytometry pipeline on one run folder. |
| `scripts/flow/validate_fcsio.py` | One-off check of fcsio.py against the fcsparser package. |
| `scripts/flow/run_flow_2026-10.sh` | Driver of the two flow-cytometry runs analysed with the pipeline. |
| `scripts/flow/final_figure/make_clone_vs_reference_figure.py` | Flow-cytometry figure: DNA-content histograms of the clone stem and two reference trees, each normalized to its own giant sequoia standard peak. The ratios and genome sizes in its legend are typed in from a table made by a script that is not on the server (see README). |
| `scripts/flow/final_figure/relax_gate_sparse.py` | Check of how the singlet gate changes the peak positions (strict versus relaxed gate; Results: the clone peak lies between 2.99x and 3.13x depending on how the peaks were found). |
| `scripts/flow/final_figure/genome_sizes.py` | Strict-gate genome sizes of the same tubes (gate comparison in the flow figure legend). Not the script behind the plotted values. |

### Supplementary Note S1 (why no somatic mutation rate)

| Script | What it does |
|---|---|
| `scripts/mutation_rate_note/mr02_scan_original.py` | Re-run of scripts/upstream/somatic_mutations_mpileup.py on the current BAMs (only the output folder changed): 607 stem-private alleles, 30-140 per stem (Supplementary Note S1). |
| `scripts/mutation_rate_note/mr02_chrom_scan.py` | Per-chromosome re-run of the stem-private search and of the same filters on the technical-replicate pairs, with the depth of every stem (Supplementary Note S1: error calibration). |
| `scripts/mutation_rate_note/mr03_combine.py` | Combines the outputs above into every number of Supplementary Note S1 and its table of expected mutation counts. |
| `scripts/mutation_rate_note/rv_detect.py` | Probability of detecting a stem-private mutation under the search filters at the stems' depths (Supplementary Note S1: about 85% at fraction 1/6, 28% at 1/12). |
| `scripts/mutation_rate_note/rv_neale_names.sh` | Shows why the re-analysis on the coast redwood assembly returned no sites: the scaffold names in the region file do not match those in the alignments (Supplementary Note S1). |
| `scripts/mutation_rate_note/rv_bam_build_check.sh` | Checks that the stem search and the replicate calibration used the same BAM files. |
| `scripts/mutation_rate_note/somatic_mutations_neale.py` | Stem-private search on the coast redwood assembly (Neale et al. 2022). Returned no sites because of the scaffold-name mismatch (Supplementary Note S1). |
| `scripts/mutation_rate_note/align_to_neale.sh` | Aligns the 13 libraries to the coast redwood assembly for the run above. |
| `scripts/mutation_rate_note/shared_somatic_variants.py` | Alleles carried by all nine stems at low fraction and rare in other trees, candidate mutations of the founding tree (Supplementary Note S1: 14 candidates, 13 inherited). |

`scripts/upstream/SOURCES.tsv` and `scripts/flow/final_figure/SOURCES.tsv` give the original location, checksum (before the header was added) and date of each copied upstream script.

## Code not included

- Read alignment to the coast redwood assembly, variant calling with Freebayes and VCF filtering: the Snakemake pipeline of the population-structure analysis (J. Santangelo).
- Allele-dosage matrix, PCR-pool correction and the original distance-matrix script: the population-structure analysis (Nikolaeva et al., in revision).
- Tree inventory (coordinates, PCR pools and groups) and the replicate/single-library allele tables read by `scripts/finalize/07_R4_G000.py`: a separate repository of the coast redwood genome project.

## Notes and limitations

- Paths are absolute paths on the analysis server; the scripts document what was run and are not a turnkey pipeline. Input data are not included.
- `scripts/upstream/idxstats_no_chlor_filtered_august.py` was edited in place after the survey run (input folder), and the survey alignment BAMs have since been deleted; `scripts/upstream/alignment_rule4.py` is identified as the alignment script from the names of its output files.
- The genome sizes and ratios in the flow-cytometry figure legend come from a table (`all_tubes_genome_sizes.csv`) whose script was run on a laptop and is not on the server; `make_clone_vs_reference_figure.py` contains these values typed in.
- The mask of nuclear sequences of plastid origin used in Supplementary Note S1 was made without a saved script.

## Software

Python 3 with numpy, pandas, scipy and matplotlib (fcsparser for `validate_fcsio.py` and the scripts in `scripts/flow/final_figure/`); BWA-MEM, samtools 1.19.2, bcftools, minimap2 and Freebayes 1.3.6.
