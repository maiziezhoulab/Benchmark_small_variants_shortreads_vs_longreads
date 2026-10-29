"""Categorical-bin plots for the tract-stratified calling metrics.

Main figures show the metric panels only; support and method details stay in TSVs / Methods.
Plotting is separate from VCF parsing so saved count tables can be reused.
No interpolation, imputation, or reweighting of the saved metrics is done.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
matplotlib.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42})
from matplotlib.lines import Line2D
from matplotlib.ticker import PercentFormatter

LENGTH_ORDER = ["<200", "200-500", "500-1k", "1k-2k", "2k-5k", "5k-10k", ">10k"]
DIST_ORDER = ["0-50", "50-125", "125-250", "250-500", "500-1k", "1k-2.5k", ">2.5k"]
REL_ORDER = ["0-0.1", "0.1-0.25", "0.25-0.5", "0.5-0.75", "0.75-1"]
DIST_DISPLAY = ["0–50", "51–125", "126–250", "251–500", "501–1000", "1001–2500", ">2500"]
REL_DISPLAY = ["[0, 0.10]", "(0.10, 0.25]", "(0.25, 0.50]", "(0.50, 0.75]", "(0.75, 1]"]
LENGTH_DISPLAY = dict(zip(LENGTH_ORDER, ["≤200", "(200, 500]", "(500, 1k]", "(1k, 2k]", "(2k, 5k]", "(5k, 10k]", ">10k"]))
LENGTH_DISPLAY.update({"0-100": "≤100", "100bp-1kb": "(100, 1k]", "1-10kb": "(1k, 10k]", ">10kb": ">10k"})
COLORS = ["#0072B2", "#56B4E9", "#009E73", "#7B4F2C", "#CC79A7", "#C78600", "#D55E00"]
MARKERS = ["o", "s", "^", "D", "v", "p", "h"]
TECH_COLORS = {"SR": "#C0392B", "LR": "#2463A5"}


def bin_spec(kind):
    if kind == "distance":
        return "dist_bin", DIST_ORDER, DIST_DISPLAY, "Distance to nearest tract edge (bp; bins)"
    return "rel_bin", REL_ORDER, REL_DISPLAY, "Relative position bin (0 = edge; 1 = midpoint)"


def metric_counts(frame, metric):
    if metric == "precision":
        return frame["n_query"]
    if metric == "recall":
        return frame["n_truth"]
    return frame[["n_truth", "n_query"]].min(axis=1, skipna=False)


def plotted_values(frame, metric, min_n):
    """Keep the complete bin grid: NaNs break lines at suppressed bins."""
    y = pd.to_numeric(frame[metric], errors="coerce").to_numpy(float)
    count = metric_counts(frame, metric).to_numpy(float)
    return np.where((count >= min_n) & np.isfinite(y), y, np.nan)


def cell_status(row, metric, min_n):
    """0 = no evaluated denominator, 1 = below threshold, 2 = shown."""
    if row is None:
        return 0
    counts = [row.get("n_query", np.nan)] if metric == "precision" else [row.get("n_truth", np.nan)]
    if metric == "f1":
        counts = [row.get("n_truth", np.nan), row.get("n_query", np.nan)]
    if not all(np.isfinite(x) and x > 0 for x in counts):
        return 0
    return 2 if min(counts) >= min_n and np.isfinite(row.get(metric, np.nan)) else 1


def _save(fig, filename):
    path = Path(filename)
    path.parent.mkdir(parents=True, exist_ok=True)
    for suffix in (".png",):
        target = path.with_suffix(suffix)
        fig.savefig(target, dpi=300, bbox_inches="tight", facecolor="white")
        print("Saved:", target)
    plt.close(fig)



def _metric_axis(ax, metric, values):
    # Uniform 0–100% ticks for all metrics in the selected lowmap/segdup export.
    # Fail rather than silently clipping data if future inputs exceed the range.
    lower = 0.0
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if len(finite) and (finite.min() < lower - 1e-12 or finite.max() > 1 + 1e-12):
        raise ValueError(f"{metric} values exceed the shared [{lower}, 1] plotting range")
    ax.set_ylim(lower, 1.025)
    ax.set_yticks(np.linspace(lower, 1, 6))
    ax.yaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))
    ax.grid(axis="y", color="#d9dfe5", linewidth=.7, alpha=.8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="y", labelsize=11)
    ax.set_axisbelow(True)


def plot_binned_metric(tab, bands, kind, metric, filename, title, min_n=50, ci=False):
    key, order, display, xlabel = bin_spec(kind)
    fig, axes = plt.subplots(1, 2, figsize=(13.4, 5.5), sharex=True, sharey=True)
    all_values = []
    visible_bands = set()
    for ax, tech in zip(axes, ("SR", "LR")):
        for i, band in enumerate(bands):
            sub = tab[(tab["tech"] == tech) & (tab["len_bin"] == band)].set_index(key).reindex(order)
            y = plotted_values(sub, metric, min_n)
            all_values.extend(y)
            if np.isfinite(y).any():
                visible_bands.add(band)
            color, marker = COLORS[i % len(COLORS)], MARKERS[i % len(MARKERS)]
            x = np.arange(len(order))
            # Keep NaNs in the line data: no interpolation across hidden bins.
            ax.plot(x, y, color=color, marker=marker, lw=1.8, ms=6,
                    markeredgewidth=.7, markeredgecolor=color, clip_on=False)
            lo, hi = f"{metric}_boot_lo", f"{metric}_boot_hi"
            if ci and lo in sub and hi in sub:
                mask = np.isfinite(y) & sub[lo].notna().to_numpy() & sub[hi].notna().to_numpy()
                ax.errorbar(x[mask], y[mask],
                    yerr=np.maximum(0, np.vstack([y[mask] - sub[lo].to_numpy()[mask], sub[hi].to_numpy()[mask] - y[mask]])),
                    fmt="none", color=color, alpha=.55, capsize=2)
        ax.set_xticks(np.arange(len(order)), display, fontsize=10)
        ax.set_xlim(-.22, len(order)-.78)
        ax.set_xlabel(xlabel, fontsize=12, labelpad=11)
        ax.set_title("Short reads (SR)" if tech == "SR" else "HiFi / long reads (LR)", fontsize=13, pad=10)
        if not any(np.isfinite(line.get_ydata()).any() for line in ax.lines):
            ax.text(.5, .5, "No supported estimates", transform=ax.transAxes, ha="center", color="#68717c")
    for ax in axes:
        _metric_axis(ax, metric, all_values)
    axes[0].set_ylabel(metric.upper() if metric == "f1" else metric.capitalize(), fontsize=13, labelpad=9)
    handles = [Line2D([], [], color=COLORS[i % len(COLORS)], marker=MARKERS[i % len(MARKERS)],
                      lw=1.8, markersize=6, label=LENGTH_DISPLAY.get(b, b))
               for i, b in enumerate(bands) if b in visible_bands]
    if handles:
        fig.legend(handles=handles, title="Tract length (bp)", loc="upper center", ncol=len(handles),
                   bbox_to_anchor=(.52, .90), frameon=False, fontsize=10, title_fontsize=10,
                   handlelength=2, columnspacing=1.7)
    fig.suptitle(f"{title} | {metric.upper() if metric == 'f1' else metric.capitalize()}", fontsize=16, y=.99)
    fig.subplots_adjust(top=.69, bottom=.16, left=.075, right=.985, wspace=.12)
    _save(fig, filename)


def plot_length_metrics(tab, bands, filename, title, min_n=50):
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 5.0))
    for ax, metric in zip(axes, ("precision", "recall", "f1")):
        values = []
        for tech in ("SR", "LR"):
            sub = tab[tab["tech"] == tech].set_index("len_bin").reindex(bands)
            y = plotted_values(sub, metric, min_n)
            values.extend(y)
            ax.plot(np.arange(len(bands)), y, color=TECH_COLORS[tech], label=tech,
                    marker="o" if tech == "SR" else "s", ms=6, lw=1.8)
        _metric_axis(ax, metric, values)
        ax.set_title(metric.upper() if metric == "f1" else metric.capitalize(), fontsize=13, pad=10)
        ax.set_xticks(np.arange(len(bands)), [LENGTH_DISPLAY.get(b, b) for b in bands], rotation=35, ha="right", fontsize=10)
        ax.set_xlabel("Tract length (bp; bins)", fontsize=12, labelpad=8)
    fig.legend(*axes[0].get_legend_handles_labels(), loc="upper center", ncol=2,
               bbox_to_anchor=(.51, .905), frameon=False, fontsize=11)
    fig.suptitle(f"{title} | Calling performance by tract length", fontsize=16, y=.995)
    fig.subplots_adjust(top=.73, bottom=.24, left=.055, right=.985, wspace=.23)
    _save(fig, filename)


def render_strata_tables(t1, t2, t3, out_prefix, vtype, title, min_n=50, ci=False, bands=None):
    """Seven main-figure families. Support counts remain in the source TSVs."""
    bands = bands or LENGTH_ORDER
    plot_length_metrics(t1, bands, f"{out_prefix}_{vtype}_metrics_by_length.png", title, min_n)
    for kind, tab in (("distance", t2), ("relpos", t3)):
        for metric in ("precision", "recall", "f1"):
            plot_binned_metric(tab, bands, kind, metric, f"{out_prefix}_{vtype}_{metric}_vs_{kind}.png", title, min_n, ci)
