#!/usr/bin/env python3

# MAPQ=0 fraction at SR-only lowmap loci, SR vs LR.
# Mirror of plot_mapq0_LRonly.py: loci short reads called and long reads did not,
# re-measured in the HiFi BAM. Writes a TSV and a histogram per variant type.

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam


# Config
LR_TSV = "/data/maiziezhou_lab/Xiwen/Projects/SNP_INDEL_LR/Platinum-Pedigree/secondary_alignment/NA24385/lowmap_tsv/NA24385_DeepVariant_L1_LR_secalign_datapoint.lowmap.tsv"
SR_TSV = "/data/maiziezhou_lab/Xiwen/Projects/SNP_INDEL_LR/Platinum-Pedigree/secondary_alignment/NA24385/lowmap_tsv/NA24385_DeepVariant_L1_SR_secalign_datapoint.lowmap.tsv"
LR_BAM = "/lfs/archer.accre.vu/maiziezhou_lab/maiziezhou_lab/Yichen/Projects/SNP_INDEL_LR/ALN_GRCh38_no_alt/Hifi/Hifi_L1.bam"
OUT_DIR = "/data/maiziezhou_lab/Xiwen/Projects/SNP_INDEL_LR/Platinum-Pedigree/secondary_alignment/NA24385/plots/groups"
OUT_TSV = os.path.join(OUT_DIR, "SRonly_with_LR_mapq_from_bam.tsv")

JOIN_KEYS = ["chrom", "pos", "ref", "alt", "variant_type"]
BINS = np.linspace(0, 1, 51)
C_LR = "#2E86C1"
C_SR = "#E67E22"


# Step 1: Find SR-only lowmap loci
def get_sr_only(lr_tsv, sr_tsv):
    lr = pd.read_csv(lr_tsv, sep="\t")
    sr = pd.read_csv(sr_tsv, sep="\t")

    lr_keys = set(zip(*[lr[k] for k in JOIN_KEYS]))
    mask = [tuple(row) not in lr_keys
               for row in zip(*[sr[k] for k in JOIN_KEYS])]
    sr_only = sr[mask].copy().reset_index(drop=True)

    print(f"LR lowmap loci:      {len(lr):,}")
    print(f"SR lowmap loci:      {len(sr):,}")
    print(f"SR-only loci:        {len(sr_only):,}  "
          f"(SNP={(sr_only['variant_type']=='SNP').sum():,}  "
          f"INDEL={(sr_only['variant_type']=='INDEL').sum():,})")
    return sr_only


# Step 2: LR mapq0_fraction from the BAM at each SR-only locus
def extract_lr_mapq(sr_only, bam_path):
    # sort to keep BAM seeks short
    sr_only = sr_only.sort_values(["chrom", "pos"]).reset_index(drop=True)

    lr_mapq0_fractions = []
    bam = pysam.AlignmentFile(bam_path, "rb")

    total = len(sr_only)
    for i, row in sr_only.iterrows():
        chrom = row["chrom"]
        pos = int(row["pos"])   # TSV pos is 1-based, fetch is 0-based
        pos0 = pos - 1

        total_reads = 0
        mapq0_reads = 0

        try:
            for read in bam.fetch(chrom, pos0, pos0 + 1):
                if read.is_unmapped or read.is_duplicate:
                    continue
                total_reads += 1
                if read.mapping_quality == 0:
                    mapq0_reads += 1
        except ValueError:          # contig absent from the HiFi BAM header
            total_reads = 0

        frac = mapq0_reads / total_reads if total_reads > 0 else np.nan
        lr_mapq0_fractions.append(frac)

        if (i + 1) % 5000 == 0:
            print(f"  processed {i+1:,} / {total:,} loci...", flush=True)

    bam.close()

    sr_only = sr_only.copy()
    sr_only["mapq0_fraction_LR_bam"] = lr_mapq0_fractions
    return sr_only


# Step 3: Plot

# SR (TSV) vs LR (BAM) histograms, SNP and INDEL panels
def plot(sr_only, out_dir):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle(
        "SR-only lowmap loci  (in SR TSV, not in LR TSV)\n"
        "SR mapq0_fraction (TSV) vs LR mapq0_fraction (extracted from BAM)\n"
        "y = number of loci, log scale",
        fontsize=12
    )

    for ax, vtype in zip(axes, ["SNP", "INDEL"]):
        sub = sr_only[sr_only["variant_type"] == vtype]
        sr_v = sub["mapq0_fraction"].dropna()
        lr_v = sub["mapq0_fraction_LR_bam"].dropna()
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

        draw(lr_v, C_LR, "LR mapq0_fraction (BAM)", zorder=2)
        draw(sr_v, C_SR, "SR mapq0_fraction (TSV)", zorder=3)

        ax.set_yscale("log")
        ax.set_xlim(0, 1)
        ax.set_xlabel("mapq0_fraction", fontsize=11)
        ax.set_ylabel("Number of variant loci", fontsize=11)
        ax.set_title(f"{vtype}  (n={n:,})", fontsize=11)
        ax.legend(fontsize=9)

    plt.tight_layout()
    path = os.path.join(out_dir, "SRonly_SR_vs_LR_bam_mapq.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: SRonly_SR_vs_LR_bam_mapq.png")


# Main
def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    print("Step 1: Finding SR-only lowmap loci...")
    sr_only = get_sr_only(LR_TSV, SR_TSV)

    print("\nStep 2: Extracting LR mapq0_fraction from BAM...")
    sr_only = extract_lr_mapq(sr_only, LR_BAM)

    sr_only.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"Saved intermediate TSV: {OUT_TSV}")

    print("\nStep 3: Plotting...")
    plot(sr_only, OUT_DIR)

    print("\nSummary:")
    for vtype in ["SNP", "INDEL"]:
        sub = sr_only[sr_only["variant_type"] == vtype]
        sr_v = sub["mapq0_fraction"].dropna()
        lr_v = sub["mapq0_fraction_LR_bam"].dropna()
        print(f"\n  {vtype} (n={len(sub):,})")
        print(f"    SR  median={np.median(sr_v):.3f}  frac_at_0={(sr_v==0).mean()*100:.1f}%")
        print(f"    LR  median={np.median(lr_v):.3f}  frac_at_0={(lr_v==0).mean()*100:.1f}%  "
              f"nan(no_coverage)={(sub['mapq0_fraction_LR_bam'].isna()).sum():,}")


if __name__ == "__main__":
    main()
