# Flow-cytometry genome-size pipeline

This pipeline takes FCS files and produces a per-tube table, QC tiers, summaries and figures. It was built on
2026-10-26 to the specification in `analysis/2026-10-26_finalize_analyses_decision.md` §2.9, which is based on
proposal §4b. Use it for every new flow run.

**Environment:**
- Run it with `~/miniconda3/envs/quast/bin/python3 -I -B` (numpy, scipy, matplotlib).
- No fcsparser is needed. `fcsio.py` is a small FCS 2.0/3.x reader. It was checked once against fcsparser 0.2.8 by
  `validate_fcsio.py`: events and TEXT keywords were identical on 3 tubes (`analysis/finalize_run/flow/validate_fcsio.log`).
- Write nothing to HOME. Set `MPLCONFIGDIR` to a scratch folder.

## Files

| File | Role |
|---|---|
| `fcsio.py` | Reads the FCS header TEXT and the event matrix |
| `flowpipe.py` | Pipeline steps 1–6, default parameters (`DEFAULTS`) and the genet-level power table |
| `make_sample_sheet.py` | Step 0. Builds a skeleton sample sheet from a folder for the lab to complete |
| `freeze_prereg.py` | Step 8. Freezes `PREREG_flow_<run>.tsv` (parameters, outcomes, power) with its sha256 and chmod 0444 |
| `run_flow.py` | Runs steps 1–7 on one run |
| `run_flow_2026-10.sh` | The two runs analysed on 2026-10-26, as worked examples |
| `validate_fcsio.py` | One-off reader validation (kmer_analysis python) |
| `sample_sheets/` | One sheet per run |

## How to run it on a new folder

Example: the lab uploads `/space/s1/sashanikolaeva/flow_cytometry_data_clone/Experiment_003/`.

```bash
cd /space/s1/sashanikolaeva/Scripts/GitHub/redwood-clone-paper
PY="$HOME/miniconda3/envs/quast/bin/python3 -I -B"
export MPLCONFIGDIR=/space/s1/sashanikolaeva/clone_finalize_work/run/mpl OMP_NUM_THREADS=1
D=/space/s1/sashanikolaeva/flow_cytometry_data_clone
RUN=2026-11_exp003          # one name per run (a run may span several days)
OUT=analysis/finalize_run/flow
FIG=/space/s1/sashanikolaeva/Scripts/GitHub/clone_paper_figures/flow

# 0. skeleton sample sheet (only the files of this run: --scan globs relative to $D)
$PY scripts/flow/make_sample_sheet.py $D $RUN scripts/flow/sample_sheets/$RUN.tsv --standard SEGI --scan "Experiment_003/**/*.fcs"
#    -> the lab / Sasha fill in tree, tree_class, tissue, buffer, chop_id, rerun_of, blinded_code (see below)

# 8. freeze the run's PREREG BEFORE running, and in a blinded design BEFORE the tube key is unblinded
$PY scripts/flow/freeze_prereg.py $RUN $OUT --purpose "designed run: clone vs euploid"   # optional: key=value overrides

# 1-7. run
$PY scripts/flow/run_flow.py --data $D --sheet scripts/flow/sample_sheets/$RUN.tsv --run $RUN --out $OUT --figdir $FIG \
    --scan "Experiment_003/**/*.fcs"
```

**Outputs:**
- In `$OUT/$RUN/`:
  - `ingest.tsv`: header voltages, events, DNA channel, ceiling shares;
  - `standards.tsv`: standard-only G1 and singlet window;
  - `per_tube.tsv`: fits in both windows, ratio, suppression, G2/G1, tier and reasons, flags, provisional 2C;
  - `per_stratum.tsv`, `per_tree.tsv`, `contrast.tsv`;
  - `run.log`.
- In `$FIG/$RUN/`:
  - one histogram grid per stratum;
  - `ratio_by_stratum_and_suppression.png`.

**Safeguards:**
- `run_flow.py` refuses to run without a frozen PREREG.
- It refuses any FCS file in the scanned folders that is not in the sheet. Byte-identical copies of listed files are
  allowed.
- `make_sample_sheet.py` and `freeze_prereg.py` never overwrite. Record changes after freezing in the run's DEVIATIONS.

## Sample sheet columns

| Column | Meaning |
|---|---|
| `file` | Path relative to `--data` |
| `tube` | Tube name (from the header) |
| `run`, `run_date` | Run name; `$DATE` from the header |
| `role` | `standard_only` (standard alone), `sample` (redwood ± standard) |
| `tree` | WP ramet ID (e.g. `28_WP7`), `22_WPR`, `16_C1` or another tree ID; `UNKNOWN` until the lab confirms |
| `tree_class` | `clone`, `euploid`, `standard` or `unknown` |
| `tissue`, `buffer` | As used; `UNKNOWN` if not known |
| `standard` | In-tube standard code (`SEGI`, `Allium`, …), or `none` |
| `chop_id` | One ID per **independent chop**. Re-runs of the same suspension share it. This is the replicate unit |
| `rerun_of` | Tube this one re-runs, if any |
| `voltage_set` | `FSC…/SSC…/DNA…` from the header (checked against the header in `ingest.tsv`) |
| `blinded_code` | Code seen by the operator |
| `dna_channel` | Optional override of the DNA channel auto-detect (BUV 496-A > Pacific Blue-A > BV421-A) |
| `include` | `yes` / `no` |
| `notes` | Free text |

**The clone-vs-euploid contrast** is computed only within a stratum, and only when every PASS-B-or-better tube in it
has a known tree, a `tree_class` of `clone` or `euploid`, and a `chop_id`. It needs one clone and at least 2 euploid
trees; the frozen design asks for 8–10. As long as any tube is `UNKNOWN`, `contrast.tsv` says "not computed".

## What the pipeline does (frozen defaults in `flowpipe.DEFAULTS`)

1. **Ingest.**
   - Reads voltages, `$TOT`, `$DATE` and `$CYT`.
   - Auto-detects the DNA channel.
   - Records the share of events at the ceiling on DNA-A, FSC-A and SSC-A.
   - Stratum = date × DNA channel × voltage set.
2. **Gates** are fixed per stratum from the standard-only tube(s):
   - debris cut at 0.5 × standard G1;
   - singlets = DNA-A/DNA-H within median ± 4 MAD.
   - **Fallback** for a stratum without a standard-only tube:
     - the reference is the standard G1 of another stratum of the run, scaled by (V/V_ref)^7.65;
     - the singlet window is set per tube;
     - the tube is flagged `gate_fallback`, its suppression is NA, and it cannot reach PASS-A.
   - **Avoid the fallback:** run a standard-only tube at every voltage setting, at the start, middle and end of each
     day.
3. **Peaks.**
   - Maxima are found on a log10 histogram.
   - Standard G1 = the most prominent maximum in 0.70–1.30 × the reference.
     - Flag `std_peak_below_window` when a more prominent maximum lies at 0.45–0.70×, i.e. suppression stronger than
       30%.
   - Redwood G1 = the most prominent maximum in 2.6–3.4 × the standard G1.
   - Each peak is fitted with a Gaussian + exponential background in two windows (±15%, ±25%).
   - G2 is searched at 2 × G1 ± 6%.
4. **Ratio and suppression.**
   - Ratio = redwood G1 / in-tube standard G1, from the ±15% window. The ±25% window is the fallback and is flagged.
   - `fit_sensitive` if the two windows differ by > 1%.
   - Suppression = in-tube / standard-only G1.
   - Provisional 2C = ratio × standard 2C. Two standard 2C values are carried:
     - the assembly (16.62 pg, i.e. 2 × 8.125 Gb);
     - a literature flow C-value (`std_2C_pg_literature`), which is `NA` until a citation is supplied.
     Both are **provisional**.
5. **QC tiers.**
   - **PASS-A** (proposal criteria), all of:
     - redwood n ≥ 1,000 and standard n ≥ 500;
     - both CVs ≤ 5%;
     - standard G2/G1 in [1.95, 2.05];
     - suppression ≥ 0.85;
     - < 1% of gated events at the ceiling;
     - converged.
   - **PASS-B**, all of:
     - ≥ 100 nuclei per peak;
     - both CVs ≤ 8%;
     - converged in the ±15% window;
     - not fit-sensitive.
   - **FAIL** otherwise. `NO-STANDARD` for tubes without an in-tube standard, which get ingest and ceiling QC only.
6. **Summaries.**
   - Per stratum, per tree (PASS-B and better, unit = chop), and the contrast.
   - Contrast: Δ% = (clone chop mean − mean of euploid tree means) / mean of euploid tree means, with a two-level
     bootstrap (trees, then chops; 2,000 reps; seed 20261026).
   - The verdict follows the PREREG outcomes:
     - **consistent:** CI excludes 0 and the point estimate is in [−2.3, −0.8]%;
     - **does not support:** CI includes 0 with half-width ≤ 0.8%;
     - **inconsistent:** Δ < −3% with the CI excluding −1.54%;
     - **inconclusive:** anything else.
7. **Figures.**

**Power** (in every frozen PREREG): one clone against k euploid trees, one-sided α 0.05, 80% power.

| k euploid trees | 6 | 8 | 10 | ∞ |
|---|---|---|---|---|
| required Δ/σ_T | 3.14 | 2.93 | 2.83 | 2.49 |
| max σ_T for Δ = 1.54% | 0.49% | 0.52% | 0.54% | 0.62% |

A designed run needs:
- 8–10 euploid trees;
- ≥ 15 independent chops of the clone over ≥ 3 days;
- one fixed voltage set per day;
- a calibrated standard.

## Runs analysed so far (2026-10-26; outputs in `analysis/finalize_run/flow/`)

| Run | Tubes | Result |
|---|---|---|
| `2025-04_berkeley` | 7 (plus 7 identical copies) | No in-tube standard, so ingest and ceiling QC only. 5 tubes at Pacific Blue 450 V (BUV 396-A 400 V); 2 at 300 V |
| `2026-10_october` | 35 (3 standard-only, 32 mixed; 14 and 19 Oct) | Methods development only: tree labels UNKNOWN, so no contrast. See `per_tube.tsv` and the results report |
