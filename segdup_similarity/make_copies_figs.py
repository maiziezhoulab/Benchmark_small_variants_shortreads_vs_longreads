#!/usr/bin/env python3

# Count how many other places each 1 kb segdup window aligns to (paf_1000.paf,
# self-hit dropped, hits within 1 kb on the same chromosome merged; capped at 50),
# then plot LR / SR recall and precision against that count and against segdup size.
# HG002, Hifi_L1 vs SR_L1, DeepVariant.

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# Config
BASE = "/data/maiziezhou_lab/Xiwen/Projects/SNP_INDEL_LR/Platinum-Pedigree/segdup_map"
SDS = "/data/maiziezhou_lab/Xiwen/Projects/SNP_INDEL_LR/Platinum-Pedigree/segdup_structure/data"

PL_EDGES = [0, 1, 2, 3, 4, 6, 10, 20, 50, 1e9]
PL_LAB = ["0", "1", "2", "3", "4-5", "6-9", "10-19", "20-49", "50+"]
SZ_EDGES = [1000, 2000, 5000, 10000, 50000, 200000, 1e12]
SZ_LAB = ["1-2k", "2-5k", "5-10k", "10-50k", "50-200k", ">200k"]
C = {"LR": "#2a78d6", "SR": "#eb6834"}
NAME = {"LR": "LR (HiFi)", "SR": "SR (Illumina)"}
INK, MUTED, GRID, SURF = "#1f1f1e", "#6b6a66", "#e6e5e1", "#fcfcfb"
plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                     "text.color": INK, "figure.facecolor": SURF, "axes.facecolor": SURF})


# Step 1: Number of distinct other places per window
def places_per_window():
    w = pd.read_pickle(f"{BASE}/data/windows_1000_parent.pkl")
    p = pd.read_pickle(f"{BASE}/data/partners_1000.pkl").sort_values(["win", "t_chrom", "t_start"])
    prev_end = p.groupby(["win", "t_chrom"]).t_end.shift()
    p["new"] = prev_end.isna() | (p.t_start > prev_end + 1000)
    w["n_places"] = w.win.map(p.groupby("win").new.sum()).fillna(0).astype(int)
    w["pl_bin"] = pd.cut(w.n_places, PL_EDGES, labels=PL_LAB, right=False)
    return w


# Step 2: Attach each call (1-based pos) to the 1 kb window containing it
def to_windows(calls, w):
    out = []
    for c, g in calls.groupby("chrom", sort=False):
        ww = w[w.chrom == c].sort_values("start")
        if ww.empty:
            continue
        p0 = g.pos.values - 1
        i = np.searchsorted(ww.start.values, p0, side="right") - 1
        ok = (i >= 0) & (p0 < ww.end.values[np.clip(i, 0, None)])
        g = g[ok].copy()
        g["pl_bin"] = ww.pl_bin.values[i[ok]]
        out.append(g)
    return pd.concat(out, ignore_index=True)


# Wilson 95% interval for k successes out of n
def wilson(k, n, z=1.96):
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return c - h, c + h


# Step 3: Rate plot, tabs = {platform: DataFrame by pl_bin with n, k}
def rate_plot(tabs, ylabel, title, path, count_note):
    fig, ax = plt.subplots(figsize=(8.6, 4.4))
    x = np.arange(len(PL_LAB))
    off = {"LR": -0.08, "SR": 0.08}
    for t in ("SR", "LR"):
        tb = tabs[t].reindex(PL_LAB)
        r = tb.k / tb.n
        lo, hi = wilson(tb.k.values, tb.n.values)
        xx = x + off[t]
        ax.vlines(xx, lo, hi, color=C[t], lw=1.2, alpha=0.6, zorder=2)
        ax.plot(xx, r, "-", color=C[t], lw=2, zorder=3, label=NAME[t])
        ax.plot(xx, r, "o", ms=8, mfc=C[t], mec=SURF, mew=2, zorder=4)
        ax.text(xx[-1] + 0.25, r.values[-1], NAME[t], va="center", fontsize=8, color=INK)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{lab}\n{count_note(lab)}" for lab in PL_LAB], fontsize=8)
    ax.set_xlim(-0.5, len(PL_LAB) + 0.6)
    ax.set_xlabel("number of other places this piece of segdup aligns to in GRCh38 (1 kb pieces)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(axis="y", color=GRID, lw=0.8, zorder=0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(frameon=False, loc="lower left")
    fig.tight_layout()
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)


# Main
def main():
    w = places_per_window()
    L = []

    # recall
    v = pd.read_pickle(f"{BASE}/data/var_win_1000.pkl")
    v = v[v.has_window].copy()
    v["pl_bin"] = pd.Categorical(v.win.map(w.set_index("win").pl_bin), categories=PL_LAB)
    rec = {t: v.groupby("pl_bin", observed=False).agg(n=(t, "size"), k=(t, "sum")) for t in ("LR", "SR")}
    n_truth = rec["LR"].n
    L.append("## Recall by number of other places (truth variants in 1 kb segdup windows)")
    L.append(pd.DataFrame({"n_truth": n_truth, "recall_LR": rec["LR"].k / n_truth,
                           "recall_SR": rec["SR"].k / n_truth}).to_string(float_format=lambda f: f"{f:.4f}"))
    rate_plot(rec, "recall", "Recall vs how many places the segdup aligns back to",
              f"{BASE}/plots/figC1_recall_vs_n_places.png",
              lambda lab: f"n={int(n_truth[lab]):,}")

    # precision
    prec, nq = {}, {}
    for t, f in (("LR", "query_LR_HifiL1_DV.tsv"), ("SR", "query_SR_L1_DV.tsv")):
        q = to_windows(pd.read_csv(f"{SDS}/{f}", sep="\t"), w)
        q["tp"] = q.bd == "TP"
        q["pl_bin"] = pd.Categorical(q.pl_bin, categories=PL_LAB)
        prec[t] = q.groupby("pl_bin", observed=False).agg(n=("tp", "size"), k=("tp", "sum"))
        prec[t]["fp"] = prec[t].n - prec[t].k
    L.append("")
    L.append("## Precision by number of other places (query calls in 1 kb segdup windows)")
    L.append(pd.DataFrame({"calls_LR": prec["LR"].n, "FP_LR": prec["LR"].fp,
                           "precision_LR": prec["LR"].k / prec["LR"].n,
                           "calls_SR": prec["SR"].n, "FP_SR": prec["SR"].fp,
                           "precision_SR": prec["SR"].k / prec["SR"].n})
             .to_string(float_format=lambda f: f"{f:.4f}"))
    rate_plot(prec, "precision", "Precision vs how many places the segdup aligns back to",
              f"{BASE}/plots/figC2_precision_vs_n_places.png",
              lambda lab: f"FP {int(prec['LR'].fp[lab])}/{int(prec['SR'].fp[lab])}")

    # places vs segdup size
    ww = w[w.parent_len.notna()].copy()
    ww["sz_bin"] = pd.cut(ww.parent_len, SZ_EDGES, labels=SZ_LAB, right=False)
    ct = pd.crosstab(ww.pl_bin, ww.sz_bin, normalize="columns").reindex(index=PL_LAB, columns=SZ_LAB)
    ncol = ww.sz_bin.value_counts().reindex(SZ_LAB)
    med = ww.groupby("sz_bin", observed=True).n_places.median().reindex(SZ_LAB)
    L.append("")
    L.append("## Share of 1 kb windows by number of other places, within each segdup size")
    L.append(ct.to_string(float_format=lambda f: f"{100*f:.1f}%"))
    L.append("windows per size: " + "  ".join(f"{k}={v:,}" for k, v in ncol.items()))
    L.append("median places:    " + "  ".join(f"{k}={v:.0f}" for k, v in med.items()))

    fig, ax = plt.subplots(figsize=(7.6, 5.2))
    im = ax.imshow(ct.values * 100, origin="lower", cmap="Blues", vmin=0, vmax=60, aspect="auto")
    for i in range(ct.shape[0]):
        for j in range(ct.shape[1]):
            val = ct.values[i, j] * 100
            ax.text(j, i, f"{val:.0f}%" if val >= 0.5 else "<1%", ha="center", va="center",
                    fontsize=8, color=SURF if val > 35 else INK)
    ax.set_xticks(range(len(SZ_LAB)))
    ax.set_xticklabels([f"{s}\n{int(ncol[s]):,} kb\nmedian {med[s]:.0f}" for s in SZ_LAB], fontsize=8)
    ax.set_yticks(range(len(PL_LAB)))
    ax.set_yticklabels(PL_LAB)
    ax.set_xlabel("size of the segdup region the 1 kb piece sits in")
    ax.set_ylabel("number of other places it aligns to")
    ax.set_title("How many places a segdup aligns back to, by segdup size\n"
                 "(each column sums to 100%)")
    for s in ax.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label("% of the column's 1 kb pieces")
    cb.outline.set_visible(False)
    fig.tight_layout()
    fig.savefig(f"{BASE}/plots/figC3_n_places_vs_segdup_size.png", dpi=160, bbox_inches="tight")
    plt.close(fig)

    # same thing as line charts
    FE = [1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 1e12]
    FL = ["1-2k", "2-5k", "5-10k", "10-20k", "20-50k", "50-100k", "100-200k", "200-500k", ">500k"]
    ww["fz"] = pd.cut(ww.parent_len, FE, labels=FL, right=False)
    f = ww.groupby("fz", observed=True).agg(
        windows=("win", "size"), sz=("parent_len", "median"),
        mean=("n_places", "mean"), med=("n_places", "median"),
        q25=("n_places", lambda s: s.quantile(.25)), q75=("n_places", lambda s: s.quantile(.75)),
        ge1=("n_places", lambda s: (s >= 1).mean()), ge2=("n_places", lambda s: (s >= 2).mean()),
        ge5=("n_places", lambda s: (s >= 5).mean()), ge10=("n_places", lambda s: (s >= 10).mean()))
    L.append("")
    L.append("## Same, as a curve over segdup size")
    L.append(f.to_string(float_format=lambda x: f"{x:.3f}"))

    # (a) mean number of lookalikes
    fig, a1 = plt.subplots(figsize=(7.8, 4.6))
    a1.plot(f.sz, f["mean"], "-", color="#2a78d6", lw=2.2, zorder=3)
    a1.plot(f.sz, f["mean"], "o", ms=9, mfc="#2a78d6", mec=SURF, mew=2, zorder=4)
    for xx, yy in zip(f.sz, f["mean"]):
        a1.annotate(f"{yy:.1f}", (xx, yy), textcoords="offset points", xytext=(0, 11),
                    ha="center", fontsize=8, color=INK)
    a1.set_ylim(0, 6.4)
    a1.set_ylabel("average number of lookalikes\n(other places a 1 kb piece aligns to)")
    a1.set_title("A piece of a big segdup has more lookalikes elsewhere in the genome\n"
                 "each point: take a random 1 kb piece from a segdup of this size,\n"
                 "how many other places does it align to on average", fontsize=9.5)

    # (b) share with at least N lookalikes
    fig2, a2 = plt.subplots(figsize=(7.8, 4.6))
    ramp = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]
    for col, lab, c in zip(["ge1", "ge2", "ge5", "ge10"],
                           ["1 or more", "2 or more", "5 or more", "10 or more"], ramp):
        a2.plot(f.sz, 100 * f[col], "-", color=c, lw=2.2, zorder=3, label=lab)
        a2.plot(f.sz, 100 * f[col], "o", ms=9, mfc=c, mec=SURF, mew=2, zorder=4)
        a2.text(f.sz.values[-1] * 1.3, 100 * f[col].values[-1], lab, va="center",
                fontsize=8.5, color=INK)
    a2.set_ylim(0, 105)
    a2.set_ylabel("% of 1 kb pieces that have this many lookalikes")
    a2.set_title("Same thing, split by how many lookalikes\n"
                 "read one line at a time: of the 1 kb pieces in a segdup of this size,\n"
                 "what share has at least N lookalikes", fontsize=9.5)
    a2.legend(title="lookalikes", frameon=False, fontsize=8.5, loc="center left")

    for a, path in ((a1, "figC4_n_places_vs_size_lines.png"),
                    (a2, "figC4b_share_by_n_lookalikes.png")):
        a.set_xscale("log")
        a.set_xlim(1100, 3.2e6)
        a.set_xticks([1000, 10000, 100000, 1000000])
        a.set_xticklabels(["1 kb", "10 kb", "100 kb", "1 Mb"])
        a.set_xlabel("size of the whole segdup region the piece comes from")
        a.grid(color=GRID, lw=0.8, zorder=0)
        for sp in ("top", "right"):
            a.spines[sp].set_visible(False)
        a.get_figure().tight_layout()
        a.get_figure().savefig(f"{BASE}/plots/{path}", dpi=160, bbox_inches="tight")
    plt.close(fig)
    plt.close(fig2)

    txt = "\n".join(L)
    open(f"{BASE}/data/report_n_places.txt", "w").write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
