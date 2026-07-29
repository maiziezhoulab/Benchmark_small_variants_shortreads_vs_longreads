#!/usr/bin/env python3

# MAPQ=0 fraction by difficult region, LR vs SR.
# Starts from the genome-wide per-site tables from compute_mapq0_fraction.py.
# For each of 5 GIAB strata (lowmap, segdup, MHC, tandem repeats + homopolymers,
# all-difficult) we intersect the tables against the region BED and keep only
# PASS calls — a DeepVariant VCF also lists RefCall sites, which are hom-ref,
# not variants. Then LR and SR are plotted against each other three ways:
# (1) a histogram split by SNP/INDEL, 
# (2) one with both pooled
# (3) a bar summary of the mean fraction plus the % of loci 

import os
import subprocess

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam


# Config
SITE_DIR = "/data/maiziezhou_lab/Xiwen/Projects/SNP_INDEL_LR/Platinum-Pedigree/secondary_alignment/NA24385"
PREFIX = "NA24385_DeepVariant_L1"

FULL_TSV = {
    "LR": f"{SITE_DIR}/{PREFIX}_LR_secalign_datapoint.tsv",
    "SR": f"{SITE_DIR}/{PREFIX}_SR_secalign_datapoint.tsv",
}

# original VCFs, used here only to get the PASS site sets
CALL_DIR = "/data/maiziezhou_lab/Yichen/Projects/SNP_INDEL_LR/GRCh38_no_alt_Call"
VCF = {
    "LR": f"{CALL_DIR}/LR/Hifi_L1/DeepVariant/output.vcf.gz",
    "SR": f"{CALL_DIR}/SR/SR_L1/DeepVariant/output.vcf.gz",
}

BED_DIR = "/data/maiziezhou_lab/Yichen/Projects/SNP_INDEL_LR/Benchmark_GRCh38/Genome_stratification_v3.6"

REGIONS = {
    "lowmap":          "GRCh38_lowmappabilityall_highconf.bed",
    "segdup":          "GRCh38_segdups_highconf.bed",
    "mhc":             "GRCh38_MHC_highconf.bed",
    "trandhomoploidy": "GRCh38_AllTandemRepeatsandHomopolymers_slop5_highconf.bed",
    "alldifficult":    "GRCh38_alldifficultregions_highconf.bed",
}

OUT_BASE = "/data/maiziezhou_lab/Xiwen/Projects/SNP_INDEL_LR/Platinum-Pedigree/secondary_alignment/NA24385/region_plots_PASS"

METRICS = {
    "mapq0_fraction": "MAPQ=0 Fraction",
}

BINS = np.linspace(0, 1, 51)
C_LR = "steelblue"
C_SR = "sandybrown"


# Step 1: Filter full TSV by BED, then by PASS sites

# (chrom, pos) of every PASS site in the VCF
#   - empty filter ('.') or PASS only - keep
#   - RefCall / LowQual etc - skip
#   store both 'chr1' and '1' so this matches whichever the TSV uses
def load_pass_sites(vcf_path):

    vcf = pysam.VariantFile(vcf_path)
    pass_sites = set()
    n_total = 0
    n_pass = 0

    for rec in vcf:
        n_total += 1
        filt = set(rec.filter.keys())
        if filt and filt != {"PASS"}:
            continue
        n_pass += 1

        chrom = str(rec.chrom)
        pos = int(rec.pos)
        pass_sites.add((chrom, pos))
        if chrom.startswith("chr"):
            pass_sites.add((chrom[3:], pos))
        else:
            pass_sites.add(("chr" + chrom, pos))

    vcf.close()
    print(f"    PASS sites: {n_pass:,} / {n_total:,} records "
          f"({100 * n_pass / max(n_total, 1):.1f}% PASS)")
    return pass_sites


def tsv_to_bed(tsv_path, bed_path):
    df = pd.read_csv(tsv_path, sep="\t", usecols=["chrom", "pos"])
    bed = pd.DataFrame({
        "chrom": df["chrom"].astype(str),
        "start": df["pos"].astype(int) - 1,
        "end":   df["pos"].astype(int),
    })
    bed.to_csv(bed_path, sep="\t", header=False, index=False)


def filter_tsv_by_bed(tsv_path, region_bed, out_tsv, tmp_dir, pass_sites):
    if os.path.exists(out_tsv):
        print(f"    [skip] already exists: {os.path.basename(out_tsv)}")
        return

    os.makedirs(tmp_dir, exist_ok=True)
    tmp_sites_bed = os.path.join(tmp_dir, os.path.basename(tsv_path) + ".sites.bed")
    tmp_intersect = os.path.join(tmp_dir, os.path.basename(tsv_path) + ".intersect.bed")

    tsv_to_bed(tsv_path, tmp_sites_bed)

    cmd = ["bedtools", "intersect", "-u", "-a", tmp_sites_bed, "-b", region_bed]
    print(f"    bedtools intersect -> {os.path.basename(tmp_intersect)}", flush=True)
    with open(tmp_intersect, "w") as f:
        subprocess.check_call(cmd, stdout=f)

    kept = pd.read_csv(tmp_intersect, sep="\t", header=None, names=["chrom", "start", "end"])
    kept_set = set(zip(kept["chrom"].astype(str), kept["end"].astype(int)))

    df = pd.read_csv(tsv_path, sep="\t")
    df["chrom"] = df["chrom"].astype(str)
    df["_key"] = list(zip(df["chrom"], df["pos"].astype(int)))

    n_before = len(df)
    # region filter
    df = df[df["_key"].isin(kept_set)]
    n_after_region = len(df)
    # PASS filter, drops RefCall
    df = df[df["_key"].isin(pass_sites)]
    n_after_pass = len(df)

    df.drop(columns="_key").to_csv(out_tsv, sep="\t", index=False)
    print(f"    {n_before:,} -> region {n_after_region:,} -> PASS {n_after_pass:,}  "
          f"-> {os.path.basename(out_tsv)}", flush=True)


# Step 2: Plots

# LR and SR in one histogram, log y because the MAPQ=0 tail sits orders of
# magnitude below the peak at 0
def overlay(ax, lr_vals, sr_vals, metric_label):
    ax.hist(lr_vals, bins=BINS, alpha=0.6, label=f"LR (n={len(lr_vals):,})", color=C_LR)
    ax.hist(sr_vals, bins=BINS, alpha=0.6, label=f"SR (n={len(sr_vals):,})", color=C_SR)
    ax.set_yscale("log")
    ax.set_xlabel(metric_label)
    ax.set_ylabel("Number of variant loci")
    ax.legend()


def save(fig, out_dir, fname):
    fig.savefig(os.path.join(out_dir, fname), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"    Saved: {fname}")


def plot_by_variant_type(lr_df, sr_df, metric, metric_label, region, out_dir):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle(f"{region} (PASS only) — LR vs SR — {metric_label}", fontsize=14)

    for ax, vtype in zip(axes, ["SNP", "INDEL"]):
        overlay(ax,
                lr_df.loc[lr_df["variant_type"] == vtype, metric].dropna(),
                sr_df.loc[sr_df["variant_type"] == vtype, metric].dropna(),
                metric_label)
        ax.set_title(vtype)

    plt.tight_layout()
    save(fig, out_dir, f"{region}_{metric}_by_vartype.png")


def plot_combined(lr_df, sr_df, metric, metric_label, region, out_dir):
    fig, ax = plt.subplots(figsize=(8, 5))
    overlay(ax, lr_df[metric].dropna(), sr_df[metric].dropna(), metric_label)
    ax.set_title(f"{region} (PASS only) — LR vs SR — {metric_label}")

    plt.tight_layout()
    save(fig, out_dir, f"{region}_{metric}_combined.png")


# left panel is the mean, right panel is how many loci are affected at all
def plot_summary_bars(lr_df, sr_df, metric, metric_label, region, out_dir):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(f"{region} (PASS only) — LR vs SR Summary — {metric_label}", fontsize=14)

    vtypes = ["SNP", "INDEL"]
    x = np.arange(len(vtypes))
    width = 0.35

    def bars(ax, lr_vals, sr_vals, ylabel, title, fmt, pad):
        ax.bar(x - width / 2, lr_vals, width, label="LR", color=C_LR)
        ax.bar(x + width / 2, sr_vals, width, label="SR", color=C_SR)
        ax.set_xticks(x)
        ax.set_xticklabels(vtypes)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend()
        for i, (lv, sv) in enumerate(zip(lr_vals, sr_vals)):
            ax.text(i - width / 2, lv + pad, fmt.format(lv), ha="center", fontsize=9)
            ax.text(i + width / 2, sv + pad, fmt.format(sv), ha="center", fontsize=9)

    def means(df):
        return [df.loc[df["variant_type"] == v, metric].mean() for v in vtypes]

    def pct_above_zero(df):
        return [100 * (df.loc[df["variant_type"] == v, metric] > 0).sum()
                / max(len(df[df["variant_type"] == v]), 1) for v in vtypes]

    bars(axes[0], means(lr_df), means(sr_df),
         f"Mean {metric_label}", "Mean Fraction", "{:.4f}", 0.005)
    bars(axes[1], pct_above_zero(lr_df), pct_above_zero(sr_df),
         "% of loci", f"% Loci with {metric_label} > 0", "{:.1f}%", 0.5)

    plt.tight_layout()
    save(fig, out_dir, f"{region}_{metric}_summary.png")


# Main
def main():
    os.makedirs(OUT_BASE, exist_ok=True)
    tmp_dir = os.path.join(OUT_BASE, "_tmp")

    # build the PASS site set once per library
    print("Loading PASS sites from VCFs ...")
    pass_sites = {}
    for lib, vcf_path in VCF.items():
        print(f"  {lib}: {vcf_path}")
        pass_sites[lib] = load_pass_sites(vcf_path)

    for region, bed_fname in REGIONS.items():
        print(f"\n{'=' * 50}")
        print(f"  Region: {region}")
        print(f"{'=' * 50}")

        region_bed = os.path.join(BED_DIR, bed_fname)
        if not os.path.exists(region_bed):
            print(f"  [WARNING] BED not found, skipping: {region_bed}")
            continue

        region_dir = os.path.join(OUT_BASE, region)
        os.makedirs(region_dir, exist_ok=True)

        tsv = {}
        for lib in ["LR", "SR"]:
            tsv[lib] = os.path.join(
                region_dir, f"{PREFIX}_{lib}_secalign_datapoint.{region}.PASS.tsv")
            print(f"  Filtering {lib} TSV (region + PASS) ...")
            filter_tsv_by_bed(FULL_TSV[lib], region_bed, tsv[lib], tmp_dir, pass_sites[lib])

        lr_df = pd.read_csv(tsv["LR"], sep="\t")
        sr_df = pd.read_csv(tsv["SR"], sep="\t")
        print(f"  LR: {len(lr_df):,} rows  |  SR: {len(sr_df):,} rows")

        for metric, metric_label in METRICS.items():
            print(f"  Plotting: {metric_label}")
            plot_by_variant_type(lr_df, sr_df, metric, metric_label, region, region_dir)
            plot_combined(lr_df, sr_df, metric, metric_label, region, region_dir)
            plot_summary_bars(lr_df, sr_df, metric, metric_label, region, region_dir)

    print("\nAll done.")


if __name__ == "__main__":
    main()
