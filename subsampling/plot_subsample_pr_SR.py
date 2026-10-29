#!/usr/bin/env python3

# Precision-recall curves across short-read subsampling depths (10x-50x).
# Reads the hap.py summary.csv of each caller at each depth and draws one line per
# caller, one marker per depth, with F1 contours. Writes full and zoomed plots.

import csv
import os

import matplotlib
matplotlib.use('agg')
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

# Config
colors = ['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',
          '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#800000', '#808000', '#ffd8b1',
          '#000075', '#808080', '#ffffff', '#000000']

tools_sr = ['DeepVariant', 'GATK', 'Strelka2']

color_dict = {tool.lower(): colors[i] for i, tool in enumerate(tools_sr)}
markers = ["o", "^", "*", "s", "H"]


# Step 1: Read hap.py results

# recall and precision of the PASS row for one variant type
def read_happy_summary(summary_csv, svtype='SNP'):
    result = {'recall': np.nan, 'precision': np.nan}
    try:
        with open(summary_csv, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                if row.get('Type') == svtype and row.get('Filter') == 'PASS':
                    result['recall']    = float(row['METRIC.Recall'])
                    result['precision'] = float(row['METRIC.Precision'])
                    break
    except Exception:
        pass
    return result


# Step 2: Plot

# PR plot with F1 contours, saved full and zoomed
def main_plot(recalls, precisions, SVtype, subsample_covs, outdir="./"):
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42

    fig, ax = plt.subplots(figsize=(10, 8), dpi=300)
    ax.set_xlim(-0.02, 1)
    ax.set_ylim(-0.05, 1.02)

    f1x = np.linspace(0.0, 1.0, num=100)
    for i in [x * 0.1 for x in range(1, 10)]:
        f1y = i * f1x / (2 * f1x - i)
        ax.plot(f1x, np.ma.masked_where(f1y <= 0, f1y), color='gray', linestyle='--', zorder=0)
        ax.annotate("F1=" + str(np.round(i, 1)), (1.0, f1y[-1] - 0.01), va='top', ha='right',
                    zorder=len(recalls.keys()) + 1)

    dot_alpha = 0.5
    lw = 1
    for j, tool in enumerate(recalls.keys()):
        ax.plot([i for i in recalls[tool] if not np.isnan(i)],
                [i for i in precisions[tool] if not np.isnan(i)],
                label=tool, zorder=j+1, color=color_dict[tool.lower()], linewidth=lw, alpha=0.9)
        for i in range(len(subsample_covs)):
            ax.scatter(recalls[tool][i], precisions[tool][i],
                       s=50, zorder=j+1, c=color_dict[tool.lower()],
                       marker=markers[i], alpha=dot_alpha)

    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")

    # tool legend, kept as an artist so the second legend does not replace it
    legend_lines = ax.legend(loc='lower left', borderaxespad=0.5, fontsize=9)
    ax.add_artist(legend_lines)

    # coverage-marker legend
    f = lambda m, c, s, l, lw: plt.plot([], [], marker=m, color=c, markersize=s, ls=l, linewidth=lw)[0]
    handles = [f(markers[i], "black", 40**0.5, '-', lw) for i in range(len(subsample_covs))]
    ax.legend(handles, subsample_covs, loc='upper left', borderaxespad=0.5, fontsize=9)

    os.makedirs(outdir, exist_ok=True)
    for ext in ['eps', 'png', 'pdf']:
        plt.savefig(os.path.join(outdir, f'{SVtype}_Recall-Precision-F1.{ext}'), bbox_inches='tight')

    ax.set_xlim(0.75, 1)
    ax.set_ylim(0.75, 1.02)
    # re-apply legends for the zoomed plot
    legend_lines2 = ax.legend(loc='lower left', borderaxespad=0.5, fontsize=9)
    ax.add_artist(legend_lines2)
    ax.legend(handles, subsample_covs, loc='upper left', borderaxespad=0.5, fontsize=9)
    for ext in ['eps', 'png', 'pdf']:
        plt.savefig(os.path.join(outdir, f'{SVtype}_Recall-Precision-F1_zoomin.{ext}'), bbox_inches='tight')

    plt.close(fig)
    print(f"  Saved plots to: {outdir}")


# Main
if __name__ == '__main__':

    tools = ['DeepVariant', 'GATK', 'Strelka2']
    subsampling = ['10x', '20x', '30x', '40x', '50x']
    SVtype = 'SNP'

    BASE_DIR = "/data/maiziezhou_lab/Yichen/Projects/SNP_INDEL_LR/GRCh38_no_alt_Call/SR_subsample"
    OUT_BASE = "/data/maiziezhou_lab/Xiwen/Projects/SNP_INDEL_LR/Platinum-Pedigree/Subsampling/SR/plots"

    regions = ['AllConf', 'AllDiff']   # the other strata were dropped

    for region in regions:
        recall = {tool: [] for tool in tools}
        precision = {tool: [] for tool in tools}

        for tool in tools:
            for sub in subsampling:
                summary_csv = os.path.join(
                    BASE_DIR,
                    f"SR_L1_{sub}",
                    tool,
                    "happy_results",
                    region,
                    f"{tool}.summary.csv"
                )
                metrics = read_happy_summary(summary_csv, svtype=SVtype)
                recall[tool].append(metrics['recall'])
                precision[tool].append(metrics['precision'])

        outdir = os.path.join(OUT_BASE, region)
        os.makedirs(outdir, exist_ok=True)

        print(f"\nPlotting region: {region}")
        main_plot(recall, precision, SVtype, subsampling, outdir=outdir)
