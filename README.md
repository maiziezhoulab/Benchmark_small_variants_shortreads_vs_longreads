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
3. [Usage](#usage)
4. [Data Availability](#data-availability)

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

## Analyses

| Directory | Question |
|---|---|
| [`radar_plots/`](radar_plots/) | How does caller accuracy differ between platforms within each region stratum? |
| [`mapq0_fraction/`](mapq0_fraction/) | How ambiguously do reads align at the loci that were called, and does that track the accuracy gap? |
| [`som/`](som/) | Do low-mappability regions of similar sequence composition share the same caller error profile? |
| [`figures/`](figures/) | Publication figures assembled from the tables the analyses above produce. |

Further components of the study — caller concordance, coverage subsampling and
B-allele frequency distributions — are not yet in this repository.

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

### Sequence composition of difficult regions

Low-mappability intervals are tiled into 500 bp windows and embedded on a
self-organizing map by 4-mer composition, so windows land next to each other
because their sequence is alike rather than because they carry the same
annotation. `hap.py`-annotated call sets are then projected onto the trained map,
giving per-neuron precision and recall alongside per-neuron GC content and 4-mer
entropy.

| Script | Description |
|---|---|
| [`split_lowmap_to_500.py`](som/split_lowmap_to_500.py) | Tiles the low-mappability BED into 500 bp windows |
| [`run_kfeat.sh`](som/run_kfeat.sh) | `laytr kfeat`, k = 4, giving a 256-dimensional composition vector per window |
| [`build_lowmap_som_15x15.py`](som/build_lowmap_som_15x15.py) | Trains a 15 × 15 hexagonal MiniSom over those vectors |
| [`map_vcf_to_som_15x15.py`](som/map_vcf_to_som_15x15.py) | Projects a `hap.py` VCF onto the map → per-neuron SNP and INDEL precision and recall |
| [`compute_GC.py`](som/compute_GC.py), [`compute_kmer_entropy.py`](som/compute_kmer_entropy.py) | Per-neuron GC content and 4-mer Shannon entropy |
| [`snp_indel_som/`](som/snp_indel_som/) | The same pipeline on windows centered on variant sites instead of a full tiling |

[`00_config.sh`](som/00_config.sh) holds the shared paths and output prefixes;
[`run_som_4mer_500bp.sh`](som/run_som_4mer_500bp.sh) drives training and mapping
end to end, and [`run_map_only.sh`](som/run_map_only.sh) projects an additional
call set onto an already-trained map.

### Figures

`collect_*.py` build tidy tables under [`figures/data/`](figures/data/), `fig_*.py`
render them to `figures/output/`, and `compose_panels.py` assembles multi-panel
layouts. The tables are committed, so the figures can be regenerated without
access to the alignments or call sets.

## Usage

Paths to alignments, call sets and region BEDs are set in a `Config` block at the
top of each Python script, or in [`som/00_config.sh`](som/00_config.sh) for the
SOM pipeline, and need to be edited before running.

```bash
# Per-site MAPQ=0 tables, one per library. Cluster job: ~6 h / 100 GB for 8.1M
# HiFi loci on 8 cores. Add --max-variants for a test run.
python mapq0_fraction/compute_mapq0_fraction.py --library LR \
    --vcf Hifi_L1/DeepVariant/output.vcf.gz --bam Hifi_L1.bam --out LR.tsv
python mapq0_fraction/plot_mapq0_by_region.py

# SOM: tile, featurize, train, then project a call set onto the trained map
python som/split_lowmap_to_500.py $LOWMAP_BED lowmap_500_len4plus.bed
bash   som/run_kfeat.sh
sbatch som/run_som_4mer_500bp.sh
sbatch som/run_map_only.sh happy_results/LowMap/DeepVariant.vcf.gz lowmapDV_SR

# Collect into tidy tables, then render the publication figures
python figures/collect_mapq0_NA24385.py && python figures/fig_mapq0_NA24385.py
python figures/collect_radar_NA24385.py && python figures/fig_radar_NA24385.py
```

Requires Python 3 with `pysam`, `pandas`, `numpy`, `matplotlib`, `joblib`,
`minisom`, `laytr` and `Pillow`, plus `bedtools` on `PATH`. `hap.py` runs in a
separate conda environment.

## Data Availability

Truth sets and region stratifications are from
[GIAB](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/). Alignments and
caller output are held on lab storage and are not distributed here.

Two caveats apply to the committed figures: the ONT panels use library `Nano_L1`,
which is degraded (13.7 % error rate, 61 % of reads unmapped) and is not
representative of ONT performance; and ONT DeepVariant was not run, so the ONT
panels show three callers where HiFi shows four.

### Citation

Manuscript in preparation.
