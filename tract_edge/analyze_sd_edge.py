#!/usr/bin/env python3

# Where a read sits relative to the segdup boundary, and what MAPQ it gets.
# Each read gets one position: "in" (inside one segdup, bp to the nearer edge),
# "out-reach" (crosses an edge, bp of unique sequence outside) or "unique" (no segdup).
# Illumina is also scored per fragment (read + mate) since bwa lets a unique mate lift MAPQ.
#
#   python analyze_sd_edge.py --dir data --segdup segdups.bed --report r.txt --out-prefix x

import argparse

import numpy as np
import pandas as pd


# Config
MAIN = set(f"chr{i}" for i in range(1, 23)) | {"chrX", "chrY"}
MIN_LEN = 1000
SIZE_EDGES = [1000, 5000, 10000, 50000, 200000, 1e12]
SIZE_LAB = ["1-5kb", "5-10kb", "10-50kb", "50-200kb", ">200kb"]
GAP_EDGES = [0, 50, 100, 250, 500, 1000, 2500, 5000, 10000, 25000, 50000, 1e12]
GAP_LAB = ["0-50", "50-100", "100-250", "250-500", "500-1k", "1k-2.5k",
           "2.5k-5k", "5k-10k", "10k-25k", "25k-50k", ">50k"]
ANC_EDGES = [0, 50, 100, 250, 500, 1000, 2500, 5000, 10000, 1e12]
ANC_LAB = ["1-50", "50-100", "100-250", "250-500", "500-1k", "1k-2.5k",
           "2.5k-5k", "5k-10k", ">10k"]
OUT_EDGES = [0, 1000, 1e12]
OUT_LAB = ["<1k", ">=1k"]
# left -> right: deep inside ... edge ... reaching out ... fully outside
ORDER = ([f"in {g}" for g in GAP_LAB[::-1]] + [f"out-reach {a}" for a in ANC_LAB]
         + [f"unique {o}" for o in OUT_LAB])


# Step 1: Load and merge segdup intervals

# merged intervals per main chromosome, shorter than MIN_LEN dropped
def load_merged(path):

    iv = {}
    for line in open(path):
        if line.startswith(("#", "track", "browser")):
            continue
        c, s, e = line.split("\t")[:3]
        if c in MAIN:
            iv.setdefault(c, []).append((int(s), int(e)))
    out = {}
    for c, v in iv.items():
        v.sort()
        m = []
        for s, e in v:
            if m and s <= m[-1][1]:
                m[-1][1] = max(m[-1][1], e)
            else:
                m.append([s, e])
        m = np.array([x for x in m if x[1] - x[0] >= MIN_LEN], dtype=np.int64)
        out[c] = m
    return out


# Step 2: Read position relative to segdup edges

# kind, distance and segdup size for half-open spans [a, b)
def geometry(iv, chrom, a, b):

    m = iv[chrom]
    st, en = m[:, 0], m[:, 1]
    ln = en - st
    cum = np.concatenate([[0], np.cumsum(ln)])
    n = len(st)
    j = np.searchsorted(en, a, side="right")          # first segdup ending after a
    i = np.searchsorted(st, b, side="left") - 1       # last segdup starting before b
    nhit = i - j + 1
    jc, ic = np.clip(j, 0, n - 1), np.clip(i, 0, n - 1)

    kind = np.full(len(a), "out", dtype=object)
    val = np.zeros(len(a), dtype=np.int64)
    size = np.zeros(len(a), dtype=np.int64)

    hit = nhit >= 1
    al = np.maximum(0, st[jc] - a)
    ar = np.maximum(0, b - en[ic])
    cov = (cum[ic + 1] - cum[jc]) - np.maximum(0, a - st[jc]) - np.maximum(0, en[ic] - b)
    uniq_mid = (b - a) - cov - al - ar                 # unique sequence between two segdups
    inside = hit & (nhit == 1) & (al == 0) & (ar == 0)
    cross = hit & ~inside

    kind[inside] = "in"
    val[inside] = np.minimum(a - st[jc], en[jc] - b)[inside]
    kind[cross] = "out-reach"
    val[cross] = np.maximum(np.maximum(al, ar), uniq_mid)[cross]
    size[hit & (nhit == 1)] = ln[jc][hit & (nhit == 1)]
    for k in np.where(hit & (nhit > 1))[0]:
        size[k] = ln[j[k]:i[k] + 1].max()

    miss = ~hit                                        # i == j - 1: neighbours on each side
    left_d = np.where(i >= 0, a - en[ic], 10**12)
    right_d = np.where(j < n, st[jc] - b, 10**12)
    use_left = left_d <= right_d
    val[miss] = np.where(use_left, left_d, right_d)[miss]
    size[miss] = np.where(use_left, ln[ic], ln[jc])[miss]
    return kind, val, size


# bin (kind, distance) into the labels in ORDER
def label(kind, val):

    lab = np.empty(len(kind), dtype=object)
    for k, edges, labs, pre in (("in", GAP_EDGES, GAP_LAB, "in "),
                                ("out-reach", ANC_EDGES, ANC_LAB, "out-reach "),
                                ("out", OUT_EDGES, OUT_LAB, "unique ")):
        s = kind == k
        idx = np.searchsorted(edges, val[s], side="right" if k != "out-reach" else "left") - 1
        idx = np.clip(idx, 0, len(labs) - 1)
        lab[s] = [pre + labs[x] for x in idx]
    return lab


# load one read table and add per-read and per-fragment positions
def load(path, iv):

    d = pd.read_csv(path, sep="\t")
    d = d.drop_duplicates(["chrom", "rs", "re", "flag", "mapq", "tlen", "clip_l", "clip_r"])
    parts = []
    for c, g in d.groupby("chrom", sort=False):
        if c not in iv:
            continue
        g = g.copy()
        k, v, s = geometry(iv, c, g.rs.values, g.re.values)
        g["kind"], g["val"], g["sd_size"] = k, v, s
        # fragment span for proper pairs, from TLEN
        prop = ((g.flag.values & 0x2) > 0) & (g.tlen.values != 0) & (np.abs(g.tlen.values) < 2000)
        fs = np.where(g.tlen.values > 0, g.rs.values, g.re.values + g.tlen.values)
        fe = np.where(g.tlen.values > 0, g.rs.values + g.tlen.values, g.re.values)
        fs, fe = np.where(prop, fs, g.rs.values), np.where(prop, fe, g.re.values)
        fk, fv, fsz = geometry(iv, c, fs, fe)
        g["f_kind"], g["f_val"], g["f_size"], g["proper"] = fk, fv, fsz, prop
        parts.append(g)
    d = pd.concat(parts, ignore_index=True)
    d["pos"] = label(d.kind.values, d.val.values)
    d["f_pos"] = label(d.f_kind.values, d.f_val.values)
    d["size_bin"] = pd.cut(d.sd_size, SIZE_EDGES, labels=SIZE_LAB, right=False)
    d["f_size_bin"] = pd.cut(d.f_size, SIZE_EDGES, labels=SIZE_LAB, right=False)
    d["mq0"] = d.mapq == 0
    d["mqlt20"] = d.mapq < 20
    return d


# Step 3: Summarise MAPQ by position

# Wilson confidence interval for a proportion
def wilson(k, n, z=1.96):

    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return c - h, c + h


# MAPQ stats per position bin, for one segdup size bin
def curve(d, pos_col, size_col, size):

    s = d[d[size_col] == size]
    g = s.groupby(pos_col).agg(n=("mapq", "size"), k0=("mq0", "sum"),
                               p0=("mq0", "mean"), p20=("mqlt20", "mean"),
                               mapq=("mapq", "mean"))
    return g.reindex([o for o in ORDER if o in g.index])


# Main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--segdup", required=True)
    ap.add_argument("--report", required=True)
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()

    iv = load_merged(a.segdup)
    D = {}
    for t in ("LR", "SR"):
        D[t] = load(f"{a.dir}/sd_edge_{t}.tsv.gz", iv)
        D[t].to_pickle(f"{a.out_prefix}_{t}.pkl")

    L = []
    def w(s=""):
        L.append(s)

    w("# MAPQ against the segdup boundary")
    w()
    w("in X        = read wholly inside a segdup, X bp from the nearer segdup edge")
    w("out-reach X = read crosses a segdup edge, X bp of unique sequence at its longest outside end")
    w("unique X    = read touches no segdup, X bp from the nearest one")
    w()
    for t in ("LR", "SR"):
        d = D[t]
        w(f"{t}: {len(d):,} primary reads after dedup, median span "
          f"{(d.re - d.rs).median():,.0f} bp")
        w("   " + "  ".join(f"{k}={v:,}" for k, v in d.kind.value_counts().items()))
    w()
    w("Intervals under 1 kb do not exist in the GIAB segdup BED (segdups are >=1 kb by")
    w("definition; the 3 records under 100 bp are 1-3 bp debris and were dropped).")
    w()

    for size in SIZE_LAB:
        w("=" * 78)
        w(f"## segdup size {size}")
        for t, pc, sc, tag in (("LR", "pos", "size_bin", "HiFi, per read"),
                               ("SR", "pos", "size_bin", "Illumina, per read"),
                               ("SR", "f_pos", "f_size_bin", "Illumina, per fragment")):
            dd = D[t] if pc == "pos" else D[t][D[t].proper]
            c = curve(dd, pc, sc, size)
            w(f"### {tag}")
            w(c[["n", "p0", "p20", "mapq"]].to_string(
                float_format=lambda x: f"{x:.4f}" if x < 1 else f"{x:.2f}"))
            w()

    txt = "\n".join(L)
    open(a.report, "w").write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
