# Benchmark_small_variants_shortreads_vs_longreads
## Background

Short-read variant calling is mature and accurate within high-confidence regions
of the genome, but those regions cover only part of it. In repetitive,
segmentally duplicated and structurally complex sequence, the short read length
restricts accurate alignment, producing ambiguous mapping, reduced recall and
elevated false discovery. Long reads span more of that sequence, and their
per-base accuracy has improved to the point where small-variant calling from
long-read data is practical.

This study benchmarks SNP and indel detection with long-read and short-read
pipelines across diverse genomic contexts, with an emphasis on difficult-to-call
regions, and characterizes the factors behind the performance difference between
the two technologies: read length, alignment accuracy, and local sequence
complexity.

## Table of Contents

1. [Benchmark Setup](#benchmark-setup)
2. [Analyses](#analyses)
3. [Methods](#methods)
4. [Usage](#usage)
5. [Data Availability](#data-availability)

## Benchmark Setup

| | Long read | Long read | Short read |
|---|---|---|---|
| Platform | PacBio HiFi | ONT | Illumina |
| Aligner | minimap2 | minimap2 | BWA-MEM |
| Callers | DeepVariant, Clair3, Longshot, NanoCaller | Clair3, Longshot, NanoCaller | DeepVariant, GATK, Strelka2 |

- **Reference:** GRCh38 no-alt analysis set
- **Samples:** NA24385 (HG002) against the GIAB v4.2.1 truth set; CEPH 1463
  pedigree (NA12877, NA12878, NA12879, NA12881, NA12882) against the
  Platinum-Pedigree v1.2 truth set
- **Region strata:** GIAB genome stratification v3.6 — `AllConf`, `AllDiff`,
  `LowMap`, `SegDup`, `MHC`, `TRandHomo`
- **Evaluation:** `hap.py`, one run per caller per stratum, `PASS` calls
  for the overview summaries. Record-level filtering for the SOM/geometry module
  is specified separately in [Methods](METHODS.md#benchmark-inputs-and-analysis-scope).

## Analyses

| Directory | Question |
|---|---|
| [`radar_plots/`](radar_plots/) | How does caller accuracy differ between platforms within each region stratum? |
| [`mapq0_fraction/`](mapq0_fraction/) | How ambiguously do reads align at the loci that were called, and does that track the accuracy gap? |
| [`som/`](som/README.md) | How do reference-sequence composition and full-tract geometry relate to SR/HiFi calling performance? |
| Publication figures | Generated from analysis tables; figure assets are not bundled in this source repository. |



### Accuracy across region strata

`hap.py` is run once per caller per stratum against the truth set restricted to
that stratum's BED. [`radar_plots/radar_plot.py`](radar_plots/radar_plot.py)
reads the resulting output tree and writes one radar plot and one CSV per
library × variant type × metric, with one axis per stratum and one line per
caller.

### Multi-mapping at called sites

MAPQ 0 means the aligner found more than one equally good placement for a read.
Measuring its frequency per locus quantifies alignment ambiguity directly rather
than inferring it from a region label, and depth at the same loci serves as a
coverage control.

| Script | Description |
|---|---|
| [`compute_mapq0_fraction.py`](mapq0_fraction/compute_mapq0_fraction.py) | VCF + BAM → per-site table of depth, MAPQ=0 count and MAPQ=0 fraction |
| [`plot_mapq0_by_region.py`](mapq0_fraction/plot_mapq0_by_region.py) | Intersects per-site tables with stratification BEDs, keeps `PASS` calls, compares LR vs SR per stratum |
| [`plot_mapq0_LRonly.py`](mapq0_fraction/plot_mapq0_LRonly.py) | Restricts to loci called by long reads only, re-measures MAPQ=0 in the short-read BAM |

### Sequence composition and tract geometry

The completed HG002 analysis uses separate low-mappability and
segmental-duplication SOMs constructed from reference 4-mer frequencies.
Existing variant-bearing windows undergo short-interval merge/expand cleanup
before feature extraction. SR and HiFi outcomes are compared in a common
regional map coordinate system, and exact paired truth records support
conditional recovery comparisons and explanatory logistic models.

A complementary analysis uses full GIAB tract boundaries to evaluate precision,
recall and F1 by tract length, distance to the nearest edge, and relative inward
position. It covers low-mappability regions, segmental duplications, all difficult
regions, tandem repeats/homopolymers and the MHC, separately for SNPs and INDELs.
The historical SOM recall convention and the final truth-side strata recall
are explicitly distinguished in the Methods.

The [SOM analysis directory](som/README.md) lists the actual scripts and input
interfaces. Methods and plots preserve the distinction between feature windows,
full tracts and benchmark-confident intervals.

## Methods

The Methods draft for the **completed SOM, paired-comparison and
tract-geometry module** is in [METHODS.md](METHODS.md), with formulas, model
definitions, parameters and bin boundaries. This contribution does not rewrite
the radar/MAPQ analyses and does not include GC heatmaps.

- [Methods implementation evidence](docs/methods_evidence.md)
- [Author checks before submission](docs/author_checks.md)
- [Exact historical neuron-support masks](docs/neuron_support_masks.tsv)

The narrative describes the executed methods. The author checks identify
training-schedule and historical-metric issues that require scientific review
before treating all legacy outputs as final manuscript results.

## Usage

For the existing MAPQ plotting scripts, edit their configuration blocks as needed.
The SOM scripts accept command-line inputs; see [som/README.md](som/README.md)
for feature extraction, exact truth pairing, tract metrics and redraw commands.

```bash
# Per-site MAPQ=0 tables, one per library. Cluster job: ~6 h / 100 GB for 8.1M
# HiFi loci on 8 cores. Add --max-variants for a test run.
python mapq0_fraction/compute_mapq0_fraction.py --library LR \
    --vcf Hifi_L1/DeepVariant/output.vcf.gz --bam Hifi_L1.bam --out LR.tsv
python mapq0_fraction/plot_mapq0_by_region.py

# Redraw completed strata tables without rereading VCFs.
python som/11_stratified_metrics.py \
  --tables-prefix "$TABLE_PREFIX" --vtype SNP --min-n 50 \
  --region-label "$REGION_LABEL" --out-prefix "$OUT_PREFIX"
```

Requires Python 3 with `pysam`, `pandas`, `numpy`, `matplotlib`, `joblib`,
`minisom`, `laytr`, `scipy`, `statsmodels` and `Pillow`, plus `bedtools` on `PATH`. `hap.py` runs in a
separate conda environment.

## Data Availability

Truth sets and region stratifications are from
[GIAB](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/). 



### Citation


