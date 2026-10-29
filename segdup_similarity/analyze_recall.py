#!/usr/bin/env python3

# Join the window partner table onto the truth variants and compute recall
# by number of partners and by identity to the closest partner, compared with
# the plain yes/no GIAB segdup label. Writes a text report.
#   python analyze_recall.py --master master.pkl --windows w.pkl --size 1000 --out var_win_1000.pkl --report r.txt

import argparse

import numpy as np
import pandas as pd


# Config
IDENT_BINS = [0, 0.90, 0.95, 0.98, 0.99, 0.995, 0.999, 1.0001]
IDENT_LAB = ["<90%", "90-95%", "95-98%", "98-99%", "99-99.5%", "99.5-99.9%", ">=99.9%"]
NPART_BINS = [-0.5, 0.5, 1.5, 2.5, 4.5, 9.5, 49.5, 1e9]
NPART_LAB = ["0", "1", "2", "3-4", "5-9", "10-49", ">=50 (capped)"]


# Step 1: Map each variant to the window containing it
def assign_windows(var, win):
    out = []
    for chrom, g in var.groupby("chrom", sort=False):
        w = win[win.chrom == chrom].sort_values("start")
        if w.empty:
            continue
        i = np.searchsorted(w.start.values, g.pos.values - 1, side="right") - 1
        ok = i >= 0
        idx = np.where(ok, i, 0)
        inside = ok & ((g.pos.values - 1) < w.end.values[idx])
        sub = g.copy()
        for c in ("win", "n_partners", "max_ident", "n_ge99", "n_chrom",
                  "psv_per_kb", "start", "end"):
            v = w[c].values[idx]
            sub[c if c not in ("start", "end") else "win_" + c] = np.where(inside, v, np.nan)
        sub["has_window"] = inside
        out.append(sub)
    return pd.concat(out, ignore_index=True)


# Step 2: LR / SR recall per group
def recall_table(d, by, labels=None):
    r = d.groupby(by, observed=True).agg(
        n=("LR", "size"), LR=("LR", "mean"), SR=("SR", "mean"))
    r["gap"] = r.LR - r.SR
    return r


# Main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--master", required=True)
    ap.add_argument("--windows", required=True)
    ap.add_argument("--size", type=int, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--report", required=True)
    a = ap.parse_args()

    var = pd.read_pickle(a.master)
    win = pd.read_pickle(a.windows)
    d = assign_windows(var, win)

    d["ident_bin"] = pd.cut(d.max_ident, IDENT_BINS, labels=IDENT_LAB, right=False)
    d["npart_bin"] = pd.cut(d.n_partners, NPART_BINS, labels=NPART_LAB)
    d.to_pickle(a.out)

    L = []
    def w(s=""):
        L.append(s)

    w(f"# Confusability map vs recall  (window size {a.size} bp)")
    w()
    w(f"truth variants            : {len(d):,}")
    w(f"inside a segdup window    : {d.has_window.sum():,}")
    w(f"GIAB in_segdup flag       : {d.in_segdup.sum():,}")
    w()

    w("## The baseline: the yes/no label")
    w(recall_table(d, "in_segdup").to_string(float_format=lambda x: f"{x:.4f}"))
    w()

    w("## Windows the aligner finds no partner for")
    w("GIAB calls these segdups. minimap2 says nothing else in the genome looks")
    w("like them at this scale. If the label were the thing that mattered, recall")
    w("here should be low too.")
    sub = d[d.has_window]
    w(recall_table(sub.assign(grp=np.where(sub.n_partners == 0,
                                           "no partner", "has partner")), "grp")
      .to_string(float_format=lambda x: f"{x:.4f}"))
    w()

    w("## Recall vs number of alternative placements")
    w(recall_table(sub, "npart_bin").to_string(float_format=lambda x: f"{x:.4f}"))
    w()

    w("## Recall vs identity of the most confusable partner")
    w(recall_table(sub[sub.n_partners > 0], "ident_bin")
      .to_string(float_format=lambda x: f"{x:.4f}"))
    w()

    w("## Same, SNPs only")
    snp = sub[(sub.n_partners > 0) & (sub.bvt == "SNP")]
    w(recall_table(snp, "ident_bin").to_string(float_format=lambda x: f"{x:.4f}"))
    w()

    w("## Distinguishing positions the read gets to see")
    w("psv_per_kb = (1 - identity) x 1000, i.e. how many letters separate this")
    w("window from its most similar competitor, per kb.")
    x = sub[sub.n_partners > 0].copy()
    x["psv_bin"] = pd.cut(x.psv_per_kb, [-0.001, 1, 2, 5, 10, 20, 50, 1e9],
                          labels=["<=1", "1-2", "2-5", "5-10", "10-20", "20-50", ">50"])
    w(recall_table(x, "psv_bin").to_string(float_format=lambda x: f"{x:.4f}"))
    w()

    w("## Aligner-measured identity vs the UCSC annotation")
    w("sd_id_hard comes from the precomputed genomicSuperDups table; max_ident is")
    w("what minimap2 actually found for this exact 1 kb of sequence.")
    both = sub[(sub.n_partners > 0) & sub.sd_id_hard.notna()]
    w(f"n = {len(both):,}   Pearson r = {both.max_ident.corr(both.sd_id_hard):.3f}")
    w(f"median |difference| = {(both.max_ident - both.sd_id_hard).abs().median():.4f}")
    w()

    w("## Two axes crossed: how many competitors x how similar")
    piv = sub[sub.n_partners > 0].pivot_table(
        index="ident_bin", columns="npart_bin", values="SR",
        aggfunc="mean", observed=True)
    w("Illumina recall")
    w(piv.to_string(float_format=lambda x: f"{x:.3f}"))
    w()
    piv = sub[sub.n_partners > 0].pivot_table(
        index="ident_bin", columns="npart_bin", values="LR",
        aggfunc="mean", observed=True)
    w("HiFi recall")
    w(piv.to_string(float_format=lambda x: f"{x:.3f}"))
    w()

    txt = "\n".join(L)
    open(a.report, "w").write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
