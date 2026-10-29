#!/usr/bin/env python3

# Recall vs identity to the most confusable partner, at 150 bp, 1 kb and 10 kb
# window sizes. Reads var_win_{150,1000,10000}.pkl from analyze_recall.py.

import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# Config
BASE = "/data/maiziezhou_lab/Xiwen/Projects/SNP_INDEL_LR/Platinum-Pedigree/segdup_map"
OUT = f"{BASE}/plots"
os.makedirs(OUT, exist_ok=True)
C = {"LR": "#1f77b4", "SR": "#d62728"}
NAME = {"LR": "HiFi", "SR": "Illumina"}
plt.rcParams.update({"font.size": 9, "axes.titlesize": 10})

IDL = ["<90", "90-95", "95-98", "98-99", "99-99.5", "99.5-99.9", ">=99.9"]
IDB = [0, .90, .95, .98, .99, .995, .999, 1.0001]


# Step 1: Plot recall vs identity, one panel per window size
fig, ax = plt.subplots(1, 3, figsize=(13, 3.9), sharey=True)
for k, S in enumerate((150, 1000, 10000)):
    v = pd.read_pickle(f"{BASE}/data/var_win_{S}.pkl")
    v = v[v.has_window & (v.n_partners > 0)]
    v = v.assign(b=pd.cut(v.max_ident, IDB, labels=IDL, right=False))
    g = v.groupby("b", observed=True).agg(n=("LR", "size"), LR=("LR", "mean"), SR=("SR", "mean"))
    x = np.arange(len(g))
    for t in ("LR", "SR"):
        ax[k].plot(x, g[t], "o-", color=C[t], label=NAME[t], ms=5)
    ax[k].set_xticks(x)
    ax[k].set_xticklabels(g.index, rotation=45, ha="right")
    ax[k].set_title(f"{S} bp window" + ("  (Illumina scale)" if S == 150 else
                    "  (HiFi scale)" if S == 10000 else ""))
    ax[k].set_xlabel("identity to most confusable partner (%)")
    ax[k].grid(alpha=.3)
    ax[k].set_ylim(0, 1.02)
ax[0].set_ylabel("recall")
ax[0].legend()
fig.suptitle("B1  Recall against the aligner's own similarity measurement", y=1.02)
fig.tight_layout()
fig.savefig(f"{OUT}/figB1_recall_vs_identity.png", dpi=150, bbox_inches="tight")
