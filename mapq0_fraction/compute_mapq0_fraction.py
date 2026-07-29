#!/usr/bin/env python3

# MAPQ=0 fraction at variant loci, one library per run.
# For every SNP/INDEL in the VCF, fetch the reads covering that position in the
# BAM and count how many have MAPQ 0. Output feeds all the downstream plots.
#
#   python compute_mapq0_fraction.py --library LR --vcf x.vcf.gz --bam x.bam --out x.tsv
#
# Secondary alignments are counted but kept out of total_mapped_depth: minimap2
# emits them for HiFi and gives them MAPQ 0, so including
# them would inflate mapq0_fraction. BWA-MEM emits none, so SR is unaffected.
#

import argparse
import os
import time
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd
import pysam


# Step 1: Load variants from VCF

# SNP or INDEL, drop symbolic alleles
def get_variant_type(ref, alts):

    if not alts:
        return None
    clean = [a for a in alts if a and not (a.startswith("<") and a.endswith(">"))]
    if not clean:
        return None
    if len(ref) == 1 and all(len(a) == 1 for a in clean):
        return "SNP"
    if any(len(a) != len(ref) for a in clean):
        return "INDEL"
    return None


# handle chr1 vs 1 naming mismatch
def match_chrom(chrom, bam_refs):

    if chrom in bam_refs:
        return chrom
    if chrom.startswith("chr") and chrom[3:] in bam_refs:
        return chrom[3:]
    if not chrom.startswith("chr") and ("chr" + chrom) in bam_refs:
        return "chr" + chrom
    return None


# all SNP/INDEL variants in the VCF, sorted
def load_variants(vcf_path, bam_path, max_variants=None):

    vcf = pysam.VariantFile(vcf_path)
    bam = pysam.AlignmentFile(bam_path, "rb")
    refs = bam.references

    variants = []
    for rec in vcf:
        vtype = get_variant_type(rec.ref, rec.alts)
        if vtype is None:
            continue
        chrom = match_chrom(str(rec.chrom), refs)
        if chrom is None:
            continue
        variants.append((chrom, int(rec.pos), rec.ref,
                         ",".join(rec.alts) if rec.alts else ".", vtype))
        if max_variants and len(variants) >= max_variants:
            break

    vcf.close()
    bam.close()

    # sort by coordinate
    variants.sort(key=lambda x: (x[0], x[1]))
    return variants


# Step 2: Split into batches
def make_batches(variants, batch_size):

    batches = []
    for i in range(0, len(variants), batch_size):
        batches.append((i // batch_size, variants[i:i + batch_size]))
    return batches


# Step 3: Worker function
# each worker opens the BAM itself and walks the batch with fetch():
#   - total_mapped_depth: primary + supplementary, secondary excluded
#   - secondary: is_secondary count, recorded but not in the denominator
#   - supplementary: is_supplementary count
#   - mapq0: how many have mapping_quality == 0
#   - mapq0_fraction: mapq0 / total_mapped_depth
def process_batch(args):

    batch_id, variant_batch, bam_path, meta = args

    bam = pysam.AlignmentFile(bam_path, "rb")
    rows = []

    for chrom, pos, ref, alt, vtype in variant_batch:
        mapped = 0
        secondary = 0
        supplementary = 0
        mapq0 = 0

        # fetch(chrom, start, end) is 0-based half-open, VCF pos is 1-based
        for read in bam.fetch(chrom, pos - 1, pos):
            if read.is_unmapped:
                continue
            if read.is_secondary:
                secondary += 1
                continue
            mapped += 1
            if read.is_supplementary:
                supplementary += 1
            if read.mapping_quality == 0:
                mapq0 += 1

        rows.append({
            **meta,
            "chrom": chrom,
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "variant_type": vtype,
            "total_mapped_depth": mapped,
            "secondary": secondary,
            "supplementary": supplementary,
            "mapq0": mapq0,
            "mapq0_fraction": (mapq0 / mapped) if mapped > 0 else None,
        })

    bam.close()
    return batch_id, rows


# Step 4: Main
def parse_args():
    p = argparse.ArgumentParser(description="Per-site MAPQ=0 fraction at variant loci")
    p.add_argument("--vcf", required=True, help="caller VCF, bgzipped + indexed")
    p.add_argument("--bam", required=True, help="the BAM the VCF was called from")
    p.add_argument("--out", required=True, help="output TSV")
    p.add_argument("--library", required=True, choices=["LR", "SR"],
                   help="label only, does not change the counting")
    p.add_argument("--sample", default="NA24385")
    p.add_argument("--caller", default="DeepVariant")
    p.add_argument("--processes", type=int, default=8)
    p.add_argument("--batch-size", type=int, default=5000)
    p.add_argument("--max-variants", type=int, default=None, help="set for testing")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

    # same on every row, rides along into the workers
    meta = {"sample": args.sample, "caller": args.caller, "library": args.library,
            "vcf": args.vcf, "bam": args.bam}

    # Read variant
    print(f"[{datetime.now()}] Loading variants from VCF ...", flush=True)
    variants = load_variants(args.vcf, args.bam, args.max_variants)
    print(f"[{datetime.now()}] Loaded {len(variants)} variants", flush=True)

    # Split batches
    batches = make_batches(variants, args.batch_size)
    print(f"[{datetime.now()}] Split into {len(batches)} batches "
          f"(batch_size={args.batch_size}, n_processes={args.processes})", flush=True)

    # multiprocessing
    all_rows = []
    done = 0
    t0 = time.time()

    with ProcessPoolExecutor(max_workers=args.processes) as pool:
        futures = [
            pool.submit(process_batch, (bid, batch, args.bam, meta))
            for bid, batch in batches
        ]

        for future in as_completed(futures):
            batch_id, rows = future.result()
            all_rows.extend(rows)
            done += 1

            if done % 50 == 0 or done == len(batches):
                elapsed = time.time() - t0
                pct = 100 * done / len(batches)
                print(f"[{datetime.now()}] Progress: {done}/{len(batches)} "
                      f"({pct:.1f}%), elapsed={elapsed / 60:.1f} min",
                      flush=True)

    # Write out results
    df = pd.DataFrame(all_rows)
    df.to_csv(args.out, sep="\t", index=False)
    print(f"[{datetime.now()}] Done. Wrote {len(df)} rows to {args.out}", flush=True)


if __name__ == "__main__":
    main()
