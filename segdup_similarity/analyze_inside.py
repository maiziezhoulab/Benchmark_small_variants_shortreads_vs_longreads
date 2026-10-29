#!/usr/bin/env python3

# Is a segdup interval uniform inside?
# Measures how much partner identity varies between windows of one GIAB interval,
# and whether position (edge vs core) or identity better predicts recall.
#   python analyze_inside.py --windows w.pkl --varwin var_win_1000.pkl --segdup segdups.bed --report r.txt --out out.pkl

import argparse
from bisect import bisect_right

import numpy as np
import pandas as pd


# Config
MAIN = set(f"chr{i}" for i in range(1, 23)) | {"chrX", "chrY"}


# Step 1: Load segdup intervals per chromosome
def load_intervals(path):
    iv = {}
    for line in open(path):
        if line.startswith(("#", "track", "browser")):
            continue
        c, s, e = line.split("\t")[:3]
        if c not in MAIN:
            continue
        iv.setdefault(c, []).append((int(s), int(e)))
    for c in iv:
        iv[c].sort()
    return iv


# Main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--windows", required=True)
    ap.add_argument("--varwin", required=True)
    ap.add_argument("--segdup", required=True)
    ap.add_argument("--report", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    iv = load_intervals(a.segdup)
    w = pd.read_pickle(a.windows)

    pid, plen, prel, pedge = [], [], [], []
    for chrom, g in w.groupby("chrom", sort=False):
        arr = iv.get(chrom, [])
        st = np.array([x[0] for x in arr])
        en = np.array([x[1] for x in arr])
        i = np.searchsorted(st, g.start.values, side="right") - 1
        i = np.clip(i, 0, max(0, len(st) - 1))
        mid = (g.start.values + g.end.values) // 2
        ok = (len(st) > 0) & (mid >= st[i]) & (mid < en[i])
        pid.extend(np.where(ok, [f"{chrom}:{x}" for x in i], None))
        L = np.where(ok, en[i] - st[i], np.nan)
        plen.extend(L)
        prel.extend(np.where(ok, (mid - st[i]) / np.maximum(1, en[i] - st[i]), np.nan))
        pedge.extend(np.where(ok, np.minimum(mid - st[i], en[i] - mid), np.nan))
    w = w.assign(parent=pid, parent_len=plen, rel_in_parent=prel, dist_edge=pedge)
    w.to_pickle(a.out)

    L = []
    def p(s=""):
        L.append(s)

    p("# Is a segdup interval uniform inside?")
    p()
    p(f"{len(w):,} windows mapped to {w.parent.nunique():,} GIAB intervals")
    p()

    p("## Spread of confusability inside one interval")
    p("For every interval with at least 10 windows, the range of max_ident across")
    p("its own windows. If intervals were uniform these would all be near zero.")
    p()
    big = w[w.parent.notna()].groupby("parent").filter(lambda g: len(g) >= 10)
    gg = big.groupby("parent").agg(
        n=("max_ident", "size"),
        lo=("max_ident", "min"), hi=("max_ident", "max"),
        iqr=("max_ident", lambda s: s.quantile(.75) - s.quantile(.25)),
        nopart=("n_partners", lambda s: (s == 0).mean()),
        plen=("parent_len", "first"))
    gg["span"] = gg.hi - gg.lo
    p(f"intervals with >=10 windows : {len(gg):,}")
    p(f"median within-interval SPAN of max_ident : {gg.span.median():.4f}")
    p(f"median within-interval IQR  of max_ident : {gg.iqr.median():.4f}")
    p(f"intervals whose windows span >5 identity points  : "
      f"{(gg.span > 0.05).mean()*100:.1f}%")
    p(f"intervals that contain BOTH a no-partner window AND a >=99.5% window : "
      f"{((gg.nopart > 0) & (gg.hi >= 0.995)).mean()*100:.1f}%")
    p()
    p("by interval size:")
    gg["size_bin"] = pd.cut(gg.plen, [0, 1e4, 5e4, 2e5, 1e9],
                            labels=["<10kb", "10-50kb", "50-200kb", ">200kb"])
    p(gg.groupby("size_bin", observed=True)
        .agg(n=("span", "size"), median_span=("span", "median"),
             median_iqr=("iqr", "median"), frac_nopart=("nopart", "mean"))
        .to_string(float_format=lambda x: f"{x:.4f}"))
    p()

    p("## Confusability vs position inside the interval")
    p("rel = 0 at one edge, 0.5 in the core. Folded, so 0 = nearest edge.")
    p()
    sub = w[w.parent.notna() & (w.parent_len >= 20000)].copy()
    sub["fold"] = np.minimum(sub.rel_in_parent, 1 - sub.rel_in_parent)
    g = sub.groupby(pd.cut(sub.fold, np.linspace(0, .5, 11)), observed=True).agg(
        n=("max_ident", "size"),
        frac_nopartner=("n_partners", lambda s: (s == 0).mean()),
        median_ident=("max_ident", "median"),
        frac_ge995=("max_ident", lambda s: (s >= 0.995).mean()))
    p(g.to_string(float_format=lambda x: f"{x:.4f}"))
    p()
    p("### the same, in bp from the nearest edge")
    g = sub.groupby(pd.cut(sub.dist_edge,
                           [0, 500, 2000, 10000, 50000, 2e5, 1e9]), observed=True).agg(
        n=("max_ident", "size"),
        frac_nopartner=("n_partners", lambda s: (s == 0).mean()),
        median_ident=("max_ident", "median"),
        frac_ge995=("max_ident", lambda s: (s >= 0.995).mean()))
    p(g.to_string(float_format=lambda x: f"{x:.4f}"))
    p()

    p("## Which predicts recall: where you sit, or what you look like?")
    v = pd.read_pickle(a.varwin)
    v = v[v.has_window].merge(
        w[["win", "parent", "parent_len", "rel_in_parent", "dist_edge"]],
        on="win", how="left")
    v = v[v.parent_len >= 20000]
    v["fold"] = np.minimum(v.rel_in_parent, 1 - v.rel_in_parent)
    v["pos_bin"] = pd.cut(v.fold, [0, .1, .25, .5], labels=["outer", "mid", "core"])
    v["id_bin"] = pd.cut(v.max_ident, [0, .98, .995, 1.001],
                         labels=["<98%", "98-99.5%", ">=99.5%"])
    p("Illumina recall, position x aligner-measured identity")
    p(v.pivot_table(index="id_bin", columns="pos_bin", values="SR",
                    aggfunc="mean", observed=True)
       .to_string(float_format=lambda x: f"{x:.3f}"))
    p()
    p("cell counts")
    p(v.pivot_table(index="id_bin", columns="pos_bin", values="SR",
                    aggfunc="size", observed=True).to_string())
    p()
    p("marginals:")
    p("  by position only : " + "  ".join(
        f"{i}={r:.3f}" for i, r in v.groupby("pos_bin", observed=True).SR.mean().items()))
    p("  by identity only : " + "  ".join(
        f"{i}={r:.3f}" for i, r in v.groupby("id_bin", observed=True).SR.mean().items()))
    p()

    txt = "\n".join(L)
    open(a.report, "w").write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
