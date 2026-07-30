#!/usr/bin/env python3
"""
Radar plots of variant calling accuracy across difficult genomic regions.

Reads hap.py benchmarking summaries and produces one plot per sequencing
library, variant type and accuracy metric. Each axis is a region stratum and
each line is a variant caller.
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

DEFAULT_LAYOUT = "{tool}/{seq_lib}/happy_results/{region}/{tool}.summary.csv"

DEFAULT_REGIONS = ["AllConf", "AllDiff", "LowMap", "MHC", "SegDup", "TRandHomo"]


def read_data(csv,
              var_type,   # SNP or INDEL
              metric,     # Recall, Precision or F1_Score
              var_filter="PASS",   # PASS or ALL
              ):
    df = pd.read_csv(csv, na_values=["", "NA"])
    subset = df.query('Type == @var_type and Filter == @var_filter')
    return float(subset[f"METRIC.{metric}"].squeeze())


def collect_data(root_dir,
                 seq_lib,
                 tools,
                 regions,
                 var_type,   # SNP or INDEL
                 metric,     # Recall, Precision or F1_Score
                 var_filter="PASS",   # PASS or ALL
                 layout=DEFAULT_LAYOUT,
                 ):
    data = dict()
    for tool in tools:
        data[tool] = dict()
        for region in regions:
            # NOTE: directory layout is not fixed; it comes from `layout`
            #       so the same script works on differently arranged trees.
            csv = os.path.join(root_dir, layout.format(tool=tool,
                                                       seq_lib=seq_lib,
                                                       region=region))
            # skip-if-missing: tolerate callers/regions that aren't benchmarked yet ---
            if not os.path.isfile(csv):
                print(f"[skip] missing file: {csv}")
                data[tool][region] = np.nan
                continue
            try:
                data[tool][region] = read_data(csv, var_type, metric, var_filter)
            except Exception as e:
                print(f"[skip] failed to read ({var_type}/{metric}) from {csv}: {e}")
                data[tool][region] = np.nan

    df = pd.DataFrame(data)
    # drop any tool column that is entirely missing (e.g. caller not run for this lane)
    df = df.dropna(axis=1, how="all")
    return df


def radar_plot(data, seq_lib, var_type, out_dir, metric):
    # ---- Compute angles ----
    labels = data.index.tolist()
    N = len(labels)
    angles = np.linspace(0, 2 * np.pi, N, endpoint=False)

    # ---- Create plot ----
    fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(polar=True))

    # Grid / ticks
    ax.set_xticks(angles)
    ax.set_xticklabels(labels)
    ax.set_rlabel_position(0)
    # yticks = np.linspace(0, 1.0, 5)[1:]
    # ax.set_yticks(yticks)
    # ax.set_yticklabels([f"{t:.2f}" for t in yticks])

    # y axis limit
    ax.set_ylim(0, 1.0)

    # ---- Plot multiple lines (one per tool) ----
    for tool in data.columns:
        values = data[tool].tolist()
        values_loop = np.r_[values, values[0]]
        angles_loop = np.r_[angles, angles[0]]
        ax.plot(angles_loop, values_loop, linewidth=2, label=tool)

    ax.set_title(f"{seq_lib} {var_type} {metric}", va="bottom")
    ax.legend(loc="upper right", bbox_to_anchor=(1.1, 1.1))

    plt.tight_layout()
    os.makedirs(out_dir, exist_ok=True)
    plt.savefig(os.path.join(out_dir, f"radar_plot_{seq_lib}_{var_type}_{metric}.png"),
                dpi=150, bbox_inches="tight")
    plt.close()


def main(seq_lib,
         tools,
         var_type,   # SNP or INDEL
         root_dir="./",
         out_dir="./",
         regions=DEFAULT_REGIONS,
         metric="F1_Score",   # Recall, Precision or F1_Score
         var_filter="PASS",   # PASS or ALL
         layout=DEFAULT_LAYOUT,
         ):
    data = collect_data(root_dir, seq_lib, tools, regions,
                        var_type, metric, var_filter, layout)
    # if no tool has any data for this combo, skip plotting entirely
    if data.empty or data.shape[1] == 0:
        print(f"[skip plot] no data for {seq_lib} {var_type} {metric}")
        return
    # save data as csv
    os.makedirs(out_dir, exist_ok=True)
    data.to_csv(os.path.join(out_dir, f"radar_plot_data_{seq_lib}_{var_type}_{metric}.csv"),
                index=True)
    radar_plot(data, seq_lib, var_type, out_dir, metric)


#  batch driver 

def parse_args():
    p = argparse.ArgumentParser(
        description="Radar plots of variant calling accuracy across difficult regions.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--root", required=True,
                   help="root of the hap.py output tree")
    p.add_argument("--out", required=True,
                   help="output directory; one subdirectory per library")
    p.add_argument("--libs", nargs="+", required=True,
                   help="sequencing libraries, one set of figures each")
    p.add_argument("--tools", nargs="+", required=True,
                   help="variant callers, one line per caller")
    p.add_argument("--regions", nargs="+", default=DEFAULT_REGIONS,
                   help="region strata, one radar axis each")
    p.add_argument("--metrics", nargs="+", default=["Recall", "Precision", "F1_Score"])
    p.add_argument("--var-types", nargs="+", default=["SNP", "INDEL"])
    p.add_argument("--filter", default="PASS", choices=["PASS", "ALL"],
                   help="hap.py Filter column to read")
    p.add_argument("--layout", default=DEFAULT_LAYOUT,
                   help="summary table path under --root")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()

    for metric in args.metrics:
        for seq_lib in args.libs:
            out_dir = os.path.join(args.out, seq_lib)
            for var_type in args.var_types:
                main(seq_lib,
                     tools=args.tools,
                     var_type=var_type,
                     root_dir=args.root,
                     out_dir=out_dir,
                     regions=args.regions,
                     metric=metric,
                     var_filter=args.filter,
                     layout=args.layout)
