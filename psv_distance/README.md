# Performance vs distance to the nearest PSV

Precision and recall of SR and LR DeepVariant calls in segmental duplications, as a
function of distance to the nearest paralogous sequence variant (PSV) on the same
duplication copy.

## Layout

- `scripts/`: PSV database and scope preparation, official hap.py re-quantification,
  per-site PSV-distance annotation, baseline summary.
- `revisions/bins3000_20261014/scripts/`: current presentation. Re-bins the saved
  distances with a 3-kb tail (1001-3000 bp, >3000 bp) and draws the figures. It imports
  from `scripts/`, so the directory structure must be kept.

## External inputs

- Indexed GRCh38 reference FASTA
- DuploMap hg38 v0.9 PSV VCF and extracted duplication database
- GIAB high-confidence, segmental-duplication and low-mappability BEDs
- Full-context hap.py benchmark VCFs for SR and LR DeepVariant
- JSON manifest (`config/inputs.json`) with `reference`, `highconf_bed`, `segdup_bed`,
  `lowmap_bed` and `callsets` (each with `id`, `tech`, `caller`, `benchmark_vcf`;
  ids `SR_DeepVariant` and `LR_DeepVariant`)
- Tools: hap.py 0.3.15 (`qfy.py`, `quantify`), bcftools
- Python >= 3.8 with pysam, numpy, pandas, matplotlib

## Run order

1. `scripts/prepare_psv_database.py`: index PSVs and copy homology
2. `scripts/prepare_scopes.py`: build evaluation scopes and distance-bin areas
3. `scripts/quantify_benchmarks.py`: re-quantify each callset with hap.py `qfy.py`
4. `scripts/annotate_callset.py`: per-site PSV distance, SNP and INDEL for each callset
5. `scripts/summarize_plot.py --skip-plots`: baseline metrics required by the revision
6. `revisions/bins3000_20261014/scripts/rebin_areas.py`: recompute bin areas at 3 kb
7. `revisions/bins3000_20261014/scripts/rebin_and_plot.py`: final tables and figures

Use `--help` on each script for its arguments.

## Notes

- Steps 3 and 4 check for `SLURM_JOB_ID` and are intended to run as Slurm jobs.
- Set the hap.py location with `--happy-bin` (step 3) and the bcftools path in the
  `NORMALIZER` constant of `annotate_callset.py` (step 4).
- Steps 5 and 7 must use the same bootstrap seed and replicate count
  (defaults 20260929 and 1000), otherwise the revision's consistency checks fail.
- The code computes both a `Segdup` scope and a nested `Segdup_lowmap` scope; the
  manuscript reports `Segdup`.
- `scripts/summarize_plot.py` and `revisions/.../summary_core_original.py` are identical
  by design; the revision imports the latter.
