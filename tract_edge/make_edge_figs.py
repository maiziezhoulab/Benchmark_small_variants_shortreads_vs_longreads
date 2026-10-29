#!/usr/bin/env python3

# Recall and MAPQ 0 rate vs distance to the nearest segdup (tract) edge, split by tract length.
# Reads master.pkl (truth variants) and sd_edge_*.tsv.gz (reads); writes figE1/figE2 PNGs.
#
#   python make_edge_figs.py [segdup|lowmap]

import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# Config
BASE = "/data/maiziezhou_lab/Xiwen/Projects/SNP_INDEL_LR/Platinum-Pedigree/segdup_map"
MASTER = "/data/maiziezhou_lab/Xiwen/Projects/SNP_INDEL_LR/Platinum-Pedigree/segdup_structure/data/master.pkl"
STRAT = "/data/maiziezhou_lab/Yichen/Projects/SNP_INDEL_LR/Benchmark_GRCh38/Genome_stratification_v3.6"
REGION = sys.argv[1] if len(sys.argv) > 1 else "segdup"
SEG, READS, SUFFIX = {
    "segdup": (f"{STRAT}/GRCh38_segdups.bed", "sd_edge", ""),
    "lowmap": (f"{STRAT}/GRCh38_lowmappabilityall.bed", "lowmap_edge", "_lowmap"),
}[REGION]
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze_sd_edge import load_merged

DIST_EDGES = [0, 50, 150, 250, 500, 1000, 2000, 5000, 1e12]
LEN_EDGES = [1000, 2000, 5000, 10000, 1e12]
LEN_LAB = ["1k-2k", "2k-5k", "5k-10k", ">10k"]
LEN_COL = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]   # ordinal: short -> long tract
MIN_N = 50
INK, MUTED, GRID = "#1f1f1e", "#6b6a66", "#e6e5e1"
plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                     "text.color": INK, "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb"})

iv = load_merged(SEG)


# Step 1: Distance to tract edge

# inside flag, distance to nearer edge and tract length for positions p
def place(chrom, p):
    m = iv[chrom]
    st, en = m[:, 0], m[:, 1]
    j = np.clip(np.searchsorted(st, p, side="right") - 1, 0, len(st) - 1)
    inside = (p >= st[j]) & (p < en[j])
    return inside, np.minimum(p - st[j], en[j] - p), en[j] - st[j]


# keep rows inside a tract and bin them by distance and tract length
def annotate(d, pos_col):

    out = []
    for c, g in d.groupby("chrom", sort=False):
        if c not in iv:
            continue
        ins, dist, L = place(c, g[pos_col].values)
        out.append(g.assign(dist=dist, tract=L)[ins])
    x = pd.concat(out, ignore_index=True)
    x["dbin"] = pd.cut(x.dist, DIST_EDGES, right=False)
    x["lbin"] = pd.cut(x.tract, LEN_EDGES, labels=LEN_LAB, right=False)
    return x


# mean of cols per (tract length, distance) bin
def table(x, cols):

    agg = {"n": (cols[0], "size"), "d_med": ("dist", "median")}
    agg.update({c: (c, "mean") for c in cols})
    return x.groupby(["lbin", "dbin"], observed=True).agg(**agg)


# Step 2: Plot

# two panels, one table per panel
def draw(tabs, cols, names, ylabel, title, path, ylim=None, legend_loc="lower left"):

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    for ax, tab, col, name in zip(axes, tabs, cols, names):
        for lab, color in zip(LEN_LAB, LEN_COL):
            if lab not in tab.index.get_level_values(0):
                continue
            t = tab.loc[lab]
            t = t[t.n >= MIN_N]
            ax.plot(t.d_med, t[col], "-", color=color, lw=2, label=lab, zorder=3)
            ax.plot(t.d_med, t[col], "o", ms=8, mfc=color, mec="#fcfcfb", mew=2, zorder=4)
        ax.set_xscale("log")
        ax.set_title(name)
        ax.set_xlabel("distance to nearest tract edge (bp)")
        ax.grid(color=GRID, lw=0.8, zorder=0)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        if ylim:
            ax.set_ylim(*ylim)
    axes[0].set_ylabel(ylabel)
    axes[1].legend(title="tract length", frameon=False, loc=legend_loc)
    fig.suptitle(title, y=1.0)
    fig.tight_layout()
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)


# Step 3: Recall
v = pd.read_pickle(MASTER)
v["p0"] = v.pos - 1
rv = annotate(v, "p0")
rt = table(rv, ["SR", "LR"])
print(f"truth variants inside tracts: {len(rv):,}")
print(rt.to_string(float_format=lambda f: f"{f:.4f}"))
draw([rt, rt], ["SR", "LR"], ["Illumina (SR_L1, DeepVariant)", "HiFi (Hifi_L1, DeepVariant)"],
     "recall", f"Recall vs distance to nearest {REGION} tract edge",
     f"{BASE}/plots/figE1_recall_vs_tract_edge{SUFFIX}.png")


# Step 4: MAPQ
tabs = {}
for t in ("SR", "LR"):
    d = pd.read_csv(f"{BASE}/data/{READS}_{t}.tsv.gz", sep="\t",
                    usecols=["chrom", "rs", "re", "mapq", "flag", "tlen", "clip_l", "clip_r"])
    d = d.drop_duplicates()
    d["mid"] = (d.rs + d.re) // 2
    x = annotate(d, "mid")
    x["mq0"] = (x.mapq == 0).astype(float)
    tabs[t] = table(x, ["mq0", "mapq"])
    print(f"\n{t}: reads with midpoint inside a tract: {len(x):,}")
    print(tabs[t].to_string(float_format=lambda f: f"{f:.4f}"))
draw([tabs["SR"], tabs["LR"]], ["mq0", "mq0"], ["Illumina (SR_L1)", "HiFi (Hifi_L1)"],
     "fraction of reads with MAPQ 0", f"MAPQ 0 rate vs distance to nearest {REGION} tract edge",
     f"{BASE}/plots/figE2_mapq0_vs_tract_edge{SUFFIX}.png", legend_loc="upper left")
print("figures written")
