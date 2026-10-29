#!/usr/bin/env python3

# Bar plot of segdup overlap per mapq0_fraction band.
# Reads sweep_summary.tsv from mapq0_threshold_sweep.py.

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# Config
SUM = "/data/maiziezhou_lab/Xiwen/Projects/SNP_INDEL_LR/Platinum-Pedigree/secondary_alignment/NA24385/mapq0_partners/threshold_sweep/sweep_summary.tsv"
OUT = "/data/maiziezhou_lab/Xiwen/Projects/SNP_INDEL_LR/Platinum-Pedigree/secondary_alignment/NA24385/mapq0_partners/threshold_sweep/segdup_overlap_by_fraction.png"


# band string to a percent label
def lab(b):
    lo, hi = b.strip("[)").split(",")
    return f"{int(round(float(lo)*100))}-{min(100, int(round(float(hi)*100)))}%"


# Plot
sm = pd.read_csv(SUM, sep="\t")
labels = [lab(b) for b in sm.band]
x = np.arange(len(sm))
w = 0.38

fig, ax = plt.subplots(figsize=(8, 5))
b1 = ax.bar(x - w / 2, sm.dest_in_segdup_pct, w, color="#E67E22", label="XA destinations in segdup")
b2 = ax.bar(x + w / 2, sm.src_in_segdup_pct, w, color="#2E86C1", label="Source loci in segdup")
ax.bar_label(b1, fmt="%.1f", fontsize=8, padding=2)
ax.bar_label(b2, fmt="%.1f", fontsize=8, padding=2)

ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.set_xlabel("mapq=0 read fraction at source locus")
ax.set_ylabel("% overlap with segmental duplications")
ax.set_ylim(0, 105)
ax.set_title("Segdup overlap across mapq0-fraction bands (NA24385 SR)")
ax.legend(loc="lower right", fontsize=9)
ax.spines[["top", "right"]].set_visible(False)
ax.grid(axis="y", color="0.9", lw=0.8)
ax.set_axisbelow(True)
plt.tight_layout()
fig.savefig(OUT, dpi=160)
print(f"-> {OUT}")
