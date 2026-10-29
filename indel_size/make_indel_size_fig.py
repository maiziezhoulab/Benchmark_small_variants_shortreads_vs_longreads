#!/usr/bin/env python3

# INDEL size effect: F1 by indel length per caller, HiFi vs Illumina.
# Reads the size_*.tsv counts from count_bd_by_size.py. Bars are GIAB truth
# INDELs per length bin, lines are F1 per caller in that bin.
# Writes the figure and a per-bin table.

import glob

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Config
B = "/data/maiziezhou_lab/Xiwen/Projects/SNP_INDEL_LR/Platinum-Pedigree/indel_size"
PANELS = [("Hifi_L1", "HiFi (Hifi_L1)"), ("SR_L1", "Illumina (SR_L1)")]
# F1 axis is zoomed, use (0, 1.03) for the full range
F1_LIM = (0.85, 1.002)
# read-based callers only
ORDER = ["Clair3", "DeepVariant", "NanoCaller", "GATK", "Strelka2"]
COL = {"Clair3": "#2a78d6", "DeepVariant": "#eb6834", "NanoCaller": "#1baf7a",
       "GATK": "#eda100", "Strelka2": "#e87ba4", "Dipcall": "#008300", "PAV": "#4a3aa7"}
ASM = {"Dipcall", "PAV"}
DEL_BAR, INS_BAR = "#d3d2ce", "#8a8984"
INK, MUTED, GRID, SURF = "#1f1f1e", "#6b6a66", "#e6e5e1", "#fcfcfb"
MIN_TRUTH = 30                       # a bin thinner than this is not plotted

# (low, high) inclusive, in bp; negative = deletion
BINS = [(-10**6, -31), (-30, -21), (-20, -16), (-15, -11), (-10, -9), (-8, -7), (-6, -5),
        (-4, -4), (-3, -3), (-2, -2), (-1, -1),
        (1, 1), (2, 2), (3, 3), (4, 4), (5, 6), (7, 8), (9, 10), (11, 15), (16, 20),
        (21, 30), (31, 10**6)]
LAB = ["≤-31", "-30~-21", "-20~-16", "-15~-11", "-10~-9", "-8~-7", "-6~-5",
       "-4", "-3", "-2", "-1",
       "1", "2", "3", "4", "5~6", "7~8", "9~10", "11~15", "16~20", "21~30", "≥31"]
plt.rcParams.update({"font.size": 8.5, "axes.titlesize": 10, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                     "text.color": INK, "figure.facecolor": SURF, "axes.facecolor": SURF})


# Step 1: Recall, precision, F1 per length bin
def binned(df):
    out = []
    for (lo, hi), lab in zip(BINS, LAB):
        m = df[(df["size"] >= lo) & (df["size"] <= hi)]
        out.append({"bin": lab, "truth_TP": m.truth_TP.sum(), "truth_FN": m.truth_FN.sum(),
                    "query_TP": m.query_TP.sum(), "query_FP": m.query_FP.sum()})
    t = pd.DataFrame(out).set_index("bin")
    t["truth"] = t.truth_TP + t.truth_FN
    t["recall"] = np.where(t.truth > 0, t.truth_TP / t.truth.replace(0, np.nan), np.nan)
    q = t.query_TP + t.query_FP
    t["precision"] = np.where(q > 0, t.query_TP / q.replace(0, np.nan), np.nan)
    s = t.recall + t.precision
    t["F1"] = np.where(s > 0, 2 * t.recall * t.precision / s.replace(0, np.nan), np.nan)
    t.loc[t.truth < MIN_TRUTH, ["recall", "precision", "F1"]] = np.nan
    return t


# Step 2: Load counts
raw = pd.concat([pd.read_csv(f, sep="\t") for f in sorted(glob.glob(f"{B}/data/size_*.tsv"))])
raw[["lib", "caller_name"]] = raw.caller.str.split("|", expand=True)
tabs = {(lib, cal): binned(g) for (lib, cal), g in raw.groupby(["lib", "caller_name"])}

# Step 3: Plot
x = np.arange(len(LAB))
nd = sum(1 for lo, hi in BINS if hi < 0)
fig, axes = plt.subplots(1, len(PANELS), figsize=(5.6 * len(PANELS), 4.8), sharex=True)
L = []
for ax, (lib, title) in zip(axes, PANELS):
    cals = [c for c in ORDER if (lib, c) in tabs]
    truth = tabs[(lib, cals[0])].truth.values
    ax.bar(x[:nd], truth[:nd], 0.82, color=DEL_BAR, zorder=1)
    ax.bar(x[nd:], truth[nd:], 0.82, color=INS_BAR, zorder=1)
    ax.set_ylim(0, 1.05 * truth.max())
    ax.axvline(nd - 0.5, color=MUTED, lw=1, ls=":", zorder=2)
    ax2 = ax.twinx()
    for c in cals:
        t = tabs[(lib, c)]
        ax2.plot(x, t.F1.values, "--" if c in ASM else "-", color=COL[c], lw=1.8,
                 marker="o", ms=3.2, label=c, zorder=3)
        L.append(t.assign(lib=lib, caller=c).reset_index())
    ax2.set_ylim(*F1_LIM)
    ax.set_title(f"{title}   {' / '.join(cals)}", fontsize=9)
    ax.set_xticks(x)
    ax.set_xticklabels(LAB, rotation=60, ha="right", fontsize=7.5)
    ax.set_xlabel("INDEL length (bp):  deletion ← | → insertion")
    ax.grid(axis="y", color=GRID, lw=0.7, zorder=0)
    for s in ("top",):
        ax.spines[s].set_visible(False)
        ax2.spines[s].set_visible(False)
    if ax is axes[0]:
        ax.set_ylabel("number of GIAB truth INDELs (bars)")
    if ax is axes[-1]:
        ax2.set_ylabel("F1 (lines)  \u2014 note the axis starts at 0.85")

drawn = [c for c in ORDER if any((lib, c) in tabs for lib, _ in PANELS)]
h = ([plt.Line2D([], [], color=COL[c], ls="--" if c in ASM else "-", lw=1.8) for c in drawn]
     + [plt.Rectangle((0, 0), 1, 1, fc=DEL_BAR), plt.Rectangle((0, 0), 1, 1, fc=INS_BAR)])
fig.legend(h, drawn + ["Deletion (bars)", "Insertion (bars)"], frameon=False,
           loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=8.5)
fig.suptitle("INDEL size effect: HG002, GIAB v4.2.1 AllConf, hap.py F1 by indel length", y=1.02)
fig.tight_layout()
fig.savefig(f"{B}/plots/figD1_indel_size_effect.png", dpi=160, bbox_inches="tight")
plt.close(fig)

# Step 4: Write table
rep = pd.concat(L, ignore_index=True)
rep.to_csv(f"{B}/data/indel_size_table.tsv", sep="\t", index=False, float_format="%.4f")
print(rep.pivot_table(index="bin", columns=["lib", "caller"], values="F1", sort=False)
         .reindex(LAB).to_string(float_format=lambda v: f"{v:.3f}"))
print("\ntruth per bin:")
print(tabs[("Hifi_L1", "Clair3")].truth.reindex(LAB).to_string())
