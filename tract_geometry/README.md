# Performance vs tract geometry

Precision, recall and F1 for SR and LR DeepVariant calls, stratified by the geometry of
the difficult-region tract each variant falls in: tract length, distance to the nearest
tract edge, and relative position within the tract (0 = edge, 1 = midpoint).

## Files

- `11_stratified_metrics.py`: driver. Reads hap.py output and a tract BED, writes
  stratified tables and calls the plotting module.
- `strata_plotting.py`: plotting module imported by the driver. Keep it in the same directory.

## Inputs

- SR and LR hap.py-annotated VCFs with `TRUTH`/`QUERY` samples and `BD` fields
- Full-tract GIAB stratification BED for the region (e.g. `GRCh38_lowmappabilityall.bed`,
  `GRCh38_segdups.bed`): the original tracts, not windowed or high-confidence-clipped intervals

## Usage

Run once per region and variant type (the manuscript uses LowMap and SegDup, SNP and INDEL):

    python 11_stratified_metrics.py \
        --sr-vcf SR.happy.vcf.gz --lr-vcf LR.happy.vcf.gz \
        --tract-bed GRCh38_lowmappabilityall.bed \
        --vtype SNP --min-n 50 --region-label LowMap \
        --out-prefix results/lowmap

`--tables-prefix` redraws figures from existing tables; `--no-plot` writes tables only.
See `python 11_stratified_metrics.py --help` for all options.

## Outputs

Three table families (`*_by_length.tsv`, `*_by_length_x_distance.tsv`,
`*_by_length_x_relpos.tsv`) and PNG figures: a length summary, and precision/recall/F1
versus edge distance and versus relative position.

## Notes

- Counts are tallied per `BD`-labelled VCF record. The script applies no FILTER/PASS
  selection, so supply VCFs filtered as intended.
- INDEL positions use the VCF `POS` anchor.
- Dependencies: numpy, pandas, pysam, matplotlib.
