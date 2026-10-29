#!/usr/bin/env python3
"""
Step 11 -- precision, recall and F1 by low-mappability tract length and by
distance from the variant to the nearest tract edge, for LR and SR.

What is new here
----------------
Everything up to step 09 was truth-side only: it read the TRUTH sample's BD
field (TP or FN) and could therefore only produce recall. False positives live
on the QUERY side, so precision and F1 need a fresh parse of both hap.py VCFs.
That is what this script does.

Plots show the main panels with shared categorical bins and metric-specific
denominator filters. Missing bins break lines; support and medians remain in TSVs.
Use --tables-prefix to redraw without parsing VCFs. PNG and PDF are saved.

Two distance scales, reported side by side
------------------------------------------
absolute (bp)   how far a read pair anchored in flanking unique sequence would
                have to reach. If read reach is the mechanism, curves from
                different tract-length bands should superimpose over the first
                few hundred bp and flatten near the insert size.
relative        dist / (tract_len / 2), so 0 is the edge and 1 is the midpoint.
                Scale-free; tests depth-within-tract rather than reach.

If the absolute curves superimpose across bands, reach explains it. If longer
tracts stay worse at the same absolute distance, tract length is carrying
something extra -- most plausibly paralog copy number.

Clustering caveat
-----------------
Variants inside one tract are not independent. Each cell reports n_tracts
alongside n_variants; when n_variants is large but n_tracts is small, treat the
confidence interval as decorative. The >10 kb band is the one to watch.

Example
-------
python 11_stratified_metrics.py \
    --lr-vcf "$LR_VCF" --sr-vcf "$SR_VCF" \
    --tract-bed $STRAT/GRCh38_lowmappabilityall.bed \
    --vtype SNP --out-prefix lowmap_strat
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
from strata_plotting import render_strata_tables

FINE_LEN_BINS = [0, 200, 500, 1000, 2000, 5000, 10000, np.inf]
FINE_LEN_LABELS = ["<200", "200-500", "500-1k", "1k-2k", "2k-5k", "5k-10k", ">10k"]
COARSE_LEN_BINS = [0, 100, 1000, 10000, np.inf]
COARSE_LEN_LABELS = ["0-100", "100bp-1kb", "1-10kb", ">10kb"]
DIST_BINS = [-1, 50, 125, 250, 500, 1000, 2500, np.inf]
DIST_LABELS = ["0-50", "50-125", "125-250", "250-500", "500-1k", "1k-2.5k", ">2.5k"]
REL_BINS = [-0.01, 0.1, 0.25, 0.5, 0.75, 1.01]
REL_LABELS = ["0-0.1", "0.1-0.25", "0.25-0.5", "0.5-0.75", "0.75-1"]


def _norm(v):
    if isinstance(v, tuple):
        v = v[0] if len(v) else None
    if isinstance(v, bytes):
        v = v.decode()
    return v


def _fmt(sample, key):
    try:
        return _norm(sample[key])
    except Exception:
        return None


def infer_vtype(ref, alts):
    alts = [a for a in alts if a not in (None, "*", "<NON_REF>")]
    if not alts:
        return None
    return "SNP" if (len(ref) == 1 and all(len(a) == 1 for a in alts)) else "INDEL"


def parse_both_sides(vcf_path, tech, vtype="ALL"):
    """One row per (record, side). side='truth' carries TP/FN for recall,
    side='query' carries TP/FP for precision."""
    rows = []
    vf = pysam.VariantFile(vcf_path)
    for rec in vf:
        alts = [a for a in (rec.alts or []) if a not in (None, "*", "<NON_REF>")]
        if not alts:
            continue
        for side, sample_name in (("truth", "TRUTH"), ("query", "QUERY")):
            if sample_name not in rec.samples:
                continue
            s = rec.samples[sample_name]
            bd = _fmt(s, "BD")
            wanted = ("TP", "FN") if side == "truth" else ("TP", "FP")
            if bd not in wanted:
                continue
            bvt = _fmt(s, "BVT")
            vt = bvt if bvt in ("SNP", "INDEL") else infer_vtype(rec.ref, alts)
            if vt is None or (vtype != "ALL" and vt != vtype):
                continue
            rows.append({"chrom": rec.chrom, "pos0": rec.pos - 1,
                         "vtype": vt, "side": side, "bd": bd, "tech": tech})
    vf.close()
    df = pd.DataFrame(rows, columns=["chrom", "pos0", "vtype", "side", "bd", "tech"])
    print(f"{tech}: truth-side {(df['side'] == 'truth').sum()}, "
          f"query-side {(df['side'] == 'query').sum()}")
    return df


def load_merged_intervals(bed_path):
    df = pd.read_csv(bed_path, sep="\t", header=None, usecols=[0, 1, 2],
                     names=["chrom", "start", "end"], comment="#")
    out = {}
    for chrom, sub in df.groupby("chrom"):
        s = sub.sort_values("start")
        starts, ends, cs, ce = [], [], None, None
        for a, b in zip(s["start"].to_numpy(int), s["end"].to_numpy(int)):
            if cs is None:
                cs, ce = a, b
            elif a <= ce:
                ce = max(ce, b)
            else:
                starts.append(cs); ends.append(ce); cs, ce = a, b
        if cs is not None:
            starts.append(cs); ends.append(ce)
        out[chrom] = (np.asarray(starts), np.asarray(ends))
    return out


def tract_geometry(iv, chrom_arr, pos_arr):
    n = len(pos_arr)
    tl = np.full(n, np.nan)
    de = np.full(n, np.nan)
    ts = np.full(n, np.nan)
    for chrom in pd.unique(chrom_arr):
        if chrom not in iv:
            continue
        m = np.where(chrom_arr == chrom)[0]
        S, E = iv[chrom]
        p = pos_arr[m]
        i = np.searchsorted(S, p, side="right") - 1
        ok = np.where((i >= 0) & (i < len(S)))[0]
        inside = ok[p[ok] < E[i[ok]]]
        if len(inside) == 0:
            continue
        gi = i[inside]
        tl[m[inside]] = E[gi] - S[gi]
        ts[m[inside]] = S[gi]
        de[m[inside]] = np.minimum(p[inside] - S[gi], E[gi] - 1 - p[inside])
    return tl, de, ts


def wilson(k, n, z=1.96):
    if n == 0:
        return np.nan, np.nan
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (c - h) / d, (c + h) / d


def cluster_boot_ci(frame, rng, n_boot):
    """Resample TRACTS with replacement and return (lo, hi) for the TP
    fraction. Variants inside one tract share a sequence context and a read
    pileup, so resampling variants would give an interval far too narrow."""
    if not n_boot or len(frame) == 0:
        return np.nan, np.nan
    g = (frame.assign(is_tp=(frame["bd"] == "TP").astype(int))
         .groupby("tract_id")["is_tp"])
    tp_arr = g.sum().to_numpy(float)
    n_arr = g.size().to_numpy(float)
    if len(tp_arr) < 2:
        return np.nan, np.nan
    idx = rng.integers(0, len(tp_arr), size=(n_boot, len(tp_arr)),
                       dtype=np.int32)
    num = tp_arr[idx].sum(axis=1)
    den = n_arr[idx].sum(axis=1)
    vals = np.divide(num, den, out=np.full(n_boot, np.nan), where=den > 0)
    lo, hi = np.nanpercentile(vals, [2.5, 97.5])
    return float(lo), float(hi)


def metrics_table(df, strata_cols, label, out_prefix, n_boot=0, seed=2023):
    rows = []
    rng = np.random.default_rng(seed)
    for keys, sub in df.groupby(strata_cols, observed=True):
        keys = keys if isinstance(keys, tuple) else (keys,)
        rec = dict(zip(strata_cols, keys))
        t = sub[sub["side"] == "truth"]
        q = sub[sub["side"] == "query"]
        tp_t = int((t["bd"] == "TP").sum())
        fn = int((t["bd"] == "FN").sum())
        tp_q = int((q["bd"] == "TP").sum())
        fp = int((q["bd"] == "FP").sum())
        recall = tp_t / (tp_t + fn) if tp_t + fn else np.nan
        precision = tp_q / (tp_q + fp) if tp_q + fp else np.nan
        f1 = (2 * precision * recall / (precision + recall)
              if np.isfinite(precision) and np.isfinite(recall)
              and (precision + recall) > 0 else np.nan)
        r_lo, r_hi = wilson(tp_t, tp_t + fn)
        p_lo, p_hi = wilson(tp_q, tp_q + fp)
        b_lo, b_hi = cluster_boot_ci(t, rng, n_boot)
        pb_lo, pb_hi = cluster_boot_ci(q, rng, n_boot)
        rec.update({
            "n_truth": tp_t + fn, "n_query": tp_q + fp,
            "n_tracts": sub["tract_id"].nunique(),
            "tp_truth": tp_t, "fn": fn, "tp_query": tp_q, "fp": fp,
            "recall": recall, "recall_lo": r_lo, "recall_hi": r_hi,
            "recall_boot_lo": b_lo, "recall_boot_hi": b_hi,
            "precision": precision, "precision_lo": p_lo, "precision_hi": p_hi,
            "precision_boot_lo": pb_lo, "precision_boot_hi": pb_hi,
            "f1": f1,
            "dist_median": float(t["dist_edge"].median()) if len(t) else np.nan,
            "rel_median": float(t["rel_pos"].median()) if len(t) else np.nan,
        })
        rows.append(rec)
    tab = pd.DataFrame(rows)
    out = f"{out_prefix}_{label}.tsv"
    tab.to_csv(out, sep="\t", index=False)
    print(f"\nSaved: {out}")
    return tab


def show(tab, index_cols, col, value, title):
    print(f"\n--- {title} ---")
    print(tab.pivot_table(index=index_cols, columns=col, values=value,
                          observed=True).round(4).to_string())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lr-vcf")
    ap.add_argument("--sr-vcf")
    ap.add_argument("--tract-bed",
                    help="original low-mappability BED, NOT _highconf")
    ap.add_argument("--vtype", default="SNP", choices=["SNP", "INDEL", "ALL"])
    ap.add_argument("--coarse", action="store_true",
                    help="use the 4 display bands instead of the 7 analysis bands")
    ap.add_argument("--min-n", type=int, default=50,
                    help="minimum metric denominator: truth for recall, query for precision, both for F1")
    ap.add_argument("--tables-prefix", help="redraw saved count TSVs without reading VCFs")
    ap.add_argument("--region-label", default="", help="region name shown in figure titles")
    ap.add_argument("--no-plot", action="store_true",
                    help="write only the TSVs, skip the figures")
    ap.add_argument("--ci", action="store_true",
                    help="draw error bars; off by default because variants "
                         "inside a tract are not independent")
    ap.add_argument("--bootstrap", type=int, default=0,
                    help="resamples for the tract-level cluster bootstrap; "
                         "500 is plenty. Use this with --ci for an interval "
                         "that respects the clustering")
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()
    if args.min_n < 1:
        ap.error("--min-n must be positive")
    Path(args.out_prefix).parent.mkdir(parents=True, exist_ok=True)
    title = f"{args.region_label or Path(args.out_prefix).name} | {args.vtype}"
    lb, ll = ((COARSE_LEN_BINS, COARSE_LEN_LABELS) if args.coarse
              else (FINE_LEN_BINS, FINE_LEN_LABELS))
    if args.tables_prefix:
        tables = [pd.read_csv(f"{args.tables_prefix}_{args.vtype}_{suffix}.tsv", sep="\t")
                  for suffix in ("by_length", "by_length_x_distance", "by_length_x_relpos")]
        render_strata_tables(*tables, args.out_prefix, args.vtype, title,
                             args.min_n, args.ci, bands=ll)
        return
    if not all((args.lr_vcf, args.sr_vcf, args.tract_bed)):
        ap.error("provide --tables-prefix or all of --lr-vcf, --sr-vcf, --tract-bed")

    df = pd.concat([parse_both_sides(args.lr_vcf, "LR", args.vtype),
                    parse_both_sides(args.sr_vcf, "SR", args.vtype)], ignore_index=True)
    if args.vtype != "ALL":
        df = df[df["vtype"] == args.vtype].copy()
    print(f"\n{args.vtype}: {len(df)} records total")

    iv = load_merged_intervals(args.tract_bed)
    tl, de, ts = tract_geometry(iv, df["chrom"].to_numpy(),
                                df["pos0"].to_numpy(int))
    df["tract_len"], df["dist_edge"], df["tract_start"] = tl, de, ts
    n_out = df["tract_len"].isna().sum()
    print(f"Outside any tract: {n_out} ({n_out / len(df):.2%})")
    df = df.dropna(subset=["tract_len"]).copy()
    df["tract_id"] = df["chrom"] + ":" + df["tract_start"].astype(int).astype(str)
    df["rel_pos"] = df["dist_edge"] / (df["tract_len"] / 2)

    lb, ll = ((COARSE_LEN_BINS, COARSE_LEN_LABELS) if args.coarse
              else (FINE_LEN_BINS, FINE_LEN_LABELS))
    df["len_bin"] = pd.cut(df["tract_len"], lb, labels=ll)
    df["dist_bin"] = pd.cut(df["dist_edge"], DIST_BINS, labels=DIST_LABELS)
    df["rel_bin"] = pd.cut(df["rel_pos"], REL_BINS, labels=REL_LABELS)

    print("\nTruth variants and distinct tracts per length band:")
    print(df[df["side"] == "truth"].groupby(["len_bin", "tech"], observed=True)
          .agg(n=("bd", "size"), tracts=("tract_id", "nunique")).to_string())

    # ---- length band only ---------------------------------------------------
    t1 = metrics_table(df, ["tech", "len_bin"], f"{args.vtype}_by_length",
                       args.out_prefix, n_boot=args.bootstrap)
    for v in ("recall", "precision", "f1"):
        show(t1, ["len_bin"], "tech", v, f"{v} by tract length")

    # ---- length band x absolute distance ------------------------------------
    t2 = metrics_table(df, ["tech", "len_bin", "dist_bin"],
                       f"{args.vtype}_by_length_x_distance", args.out_prefix,
                       n_boot=args.bootstrap)
    keep = t2[t2["n_truth"] >= args.min_n]
    for tech in ("SR", "LR"):
        show(keep[keep["tech"] == tech], ["len_bin"], "dist_bin", "recall",
             f"{tech} recall, rows = tract length, cols = ABSOLUTE distance to edge")
    show(keep[keep["tech"] == "SR"], ["len_bin"], "dist_bin", "precision",
         "SR precision, same layout")
    print("\nRead the SR recall table down each column: if the values in one "
          "distance column are similar across length bands, reach explains the "
          "effect. If they keep falling as tract length grows at the same "
          "absolute distance, tract length carries something beyond reach.")

    # ---- length band x relative position ------------------------------------
    t3 = metrics_table(df, ["tech", "len_bin", "rel_bin"],
                       f"{args.vtype}_by_length_x_relpos", args.out_prefix,
                       n_boot=args.bootstrap)
    keep3 = t3[t3["n_truth"] >= args.min_n]
    for tech in ("SR", "LR"):
        show(keep3[keep3["tech"] == tech], ["len_bin"], "rel_bin", "recall",
             f"{tech} recall, rows = tract length, cols = RELATIVE position "
             f"(0 = edge, 1 = midpoint)")
    print("\nCompare the two layouts: whichever collapses the length bands "
          "onto a common curve is the scale the mechanism actually operates "
          "on.")

    if args.no_plot:
        return

    render_strata_tables(t1, t2, t3, args.out_prefix, args.vtype, title,
                         args.min_n, args.ci, bands=ll)


if __name__ == "__main__":
    main()
