#!/usr/bin/env python3

# MAPQ=0 fraction at LR-only loci, measured in the short-read BAM.
# Takes the low-mappability loci called by long reads but not by short reads,
# re-queries the short-read BAM at those positions whether or not a variant was
# called there, and overlays the LR and SR MAPQ=0 distributions for SNP and
# INDEL. Swap LR_TSV / SR_TSV to region_plots_PASS/<region>/ for another stratum.

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam


# Config
BASE = "/data/maiziezhou_lab/Xiwen/Projects/SNP_INDEL_LR/Platinum-Pedigree/secondary_alignment/NA24385"
PREFIX = "NA24385_DeepVariant_L1"

LR_TSV = f"{BASE}/lowmap_tsv/{PREFIX}_LR_secalign_datapoint.lowmap.tsv"
SR_TSV = f"{BASE}/lowmap_tsv/{PREFIX}_SR_secalign_datapoint.lowmap.tsv"
SR_BAM = "/lfs/archer.accre.vu/maiziezhou_lab/maiziezhou_lab/Yichen/Projects/SNP_INDEL_LR/ALN_GRCh38_no_alt/SR/SR_L1.bam"

OUT_DIR = f"{BASE}/plots/groups"
OUT_TSV = os.path.join(OUT_DIR, "LRonly_with_SR_mapq_from_bam.tsv")

JOIN_KEYS = ["chrom", "pos", "ref", "alt", "variant_type"]
BINS = np.linspace(0, 1, 51)
C_LR = "#2E86C1"
C_SR = "#E67E22"


# Step 1: Find LR-only lowmap loci
def get_lr_only(lr_tsv, sr_tsv):

    lr = pd.read_csv(lr_tsv, sep="\t")
    sr = pd.read_csv(sr_tsv, sep="\t")

    sr_keys = set(zip(*[sr[k] for k in JOIN_KEYS]))
    mask = [tuple(row) not in sr_keys
            for row in zip(*[lr[k] for k in JOIN_KEYS])]
    lr_only = lr[mask].copy().reset_index(drop=True)

    print(f"LR lowmap loci:      {len(lr):,}")
    print(f"SR lowmap loci:      {len(sr):,}")
    print(f"LR-only loci:        {len(lr_only):,}  "
          f"(SNP={(lr_only['variant_type'] == 'SNP').sum():,}  "
          f"INDEL={(lr_only['variant_type'] == 'INDEL').sum():,})")
    return lr_only


# Step 2: Extract SR mapq0_fraction from BAM for each LR-only locus
# note this filters is_duplicate, which compute_mapq0_fraction.py does not.
# that is how the figure was made, left alone; only affects the orange
# distribution here.
def extract_sr_mapq(lr_only, bam_path):

    # sort first so the BAM is read front to back instead of seeking around
    lr_only = lr_only.sort_values(["chrom", "pos"]).reset_index(drop=True)

    sr_mapq0_fractions = []
    bam = pysam.AlignmentFile(bam_path, "rb")

    total = len(lr_only)
    for i, row in lr_only.iterrows():
        chrom = row["chrom"]
        pos = int(row["pos"])   # TSV pos is 1-based, fetch is 0-based
        pos0 = pos - 1

        total_reads = 0
        mapq0_reads = 0

        for read in bam.fetch(chrom, pos0, pos0 + 1):
            if read.is_unmapped or read.is_duplicate:
                continue
            total_reads += 1
            if read.mapping_quality == 0:
                mapq0_reads += 1

        # no short-read coverage at all -> nan, reported separately below
        frac = mapq0_reads / total_reads if total_reads > 0 else np.nan
        sr_mapq0_fractions.append(frac)

        if (i + 1) % 1000 == 0:
            print(f"  processed {i + 1:,} / {total:,} loci...", flush=True)

    bam.close()

    lr_only = lr_only.copy()
    lr_only["mapq0_fraction_SR_bam"] = sr_mapq0_fractions
    return lr_only


# Step 3: Plot
# blue = LR mapq0_fraction (from TSV), orange = SR mapq0_fraction (from BAM)
def plot(lr_only, out_dir):

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle(
        "LR-only lowmap loci  (in LR TSV, not in SR TSV)\n"
        "LR mapq0_fraction (TSV) vs SR mapq0_fraction (extracted from BAM)\n"
        "y = number of loci, log scale",
        fontsize=12
    )

    for ax, vtype in zip(axes, ["SNP", "INDEL"]):
        sub = lr_only[lr_only["variant_type"] == vtype]
        lr_v = sub["mapq0_fraction"].dropna()
        sr_v = sub["mapq0_fraction_SR_bam"].dropna()
        n = len(sub)

        def draw(vals, color, label, zorder):
            if len(vals) == 0:
                return
            counts, edges = np.histogram(vals, bins=BINS)
            centers = (edges[:-1] + edges[1:]) / 2
            ax.bar(centers, counts, width=edges[1] - edges[0],
                   color=color, alpha=0.6, edgecolor="white", linewidth=0.3,
                   label=f"{label}  (n={len(vals):,})", zorder=zorder)
            med = np.median(vals)
            ax.axvline(med, color=color, lw=1.5, ls="--", alpha=0.9,
                       label=f"median = {med:.3f}")

        # SR first, so the narrow LR spike at 0 stays visible on top
        draw(sr_v, C_SR, "SR mapq0_fraction (BAM)", zorder=2)
        draw(lr_v, C_LR, "LR mapq0_fraction (TSV)", zorder=3)

        ax.set_yscale("log")
        ax.set_xlim(0, 1)
        ax.set_xlabel("mapq0_fraction", fontsize=11)
        ax.set_ylabel("Number of variant loci", fontsize=11)
        ax.set_title(f"{vtype}  (n={n:,})", fontsize=11)
        ax.legend(fontsize=9)

    plt.tight_layout()
    path = os.path.join(out_dir, "LRonly_LR_vs_SR_bam_mapq.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: LRonly_LR_vs_SR_bam_mapq.png")


# Main
def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    print("Step 1: Finding LR-only lowmap loci...")
    lr_only = get_lr_only(LR_TSV, SR_TSV)

    print("\nStep 2: Extracting SR mapq0_fraction from BAM...")
    print("(this may take 10-30 minutes for ~28k loci)")
    lr_only = extract_sr_mapq(lr_only, SR_BAM)

    # save it so the next redraw skips the BAM pass
    lr_only.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"Saved intermediate TSV: {OUT_TSV}")

    print("\nStep 3: Plotting...")
    plot(lr_only, OUT_DIR)

    print("\nSummary:")
    for vtype in ["SNP", "INDEL"]:
        sub = lr_only[lr_only["variant_type"] == vtype]
        lr_v = sub["mapq0_fraction"].dropna()
        sr_v = sub["mapq0_fraction_SR_bam"]
        print(f"\n  {vtype} (n={len(sub):,})")
        print(f"    LR  median={np.median(lr_v):.3f}  "
              f"frac_at_0={(lr_v == 0).mean() * 100:.1f}%")
        print(f"    SR  median={np.median(sr_v.dropna()):.3f}  "
              f"frac_at_0={(sr_v.dropna() == 0).mean() * 100:.1f}%  "
              f"nan(no_coverage)={sr_v.isna().sum():,}")


if __name__ == "__main__":
    main()
