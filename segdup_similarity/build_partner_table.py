#!/usr/bin/env python3

# Turn the minimap2 PAF of segdup windows into a per-window partner table:
# where else each window aligns, and how similar that place is.
# The self-hit is dropped. Identity = PAF col10 / col11 (gaps counted).
#   python build_partner_table.py --paf paf_1000.paf --size 1000 --out-windows w.pkl --out-partners p.pkl

import argparse

import pandas as pd
import numpy as np


# Config
MIN_QCOV = 0.5     # a hit covering under half the window is not a real partner


# Step 1: Parse PAF records, flag the self-hit
def parse(paf, size):
    rows = []
    with open(paf) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            qname, qlen, qs, qe, strand = f[0], int(f[1]), int(f[2]), int(f[3]), f[4]
            tname, ts, te = f[5], int(f[7]), int(f[8])
            nmatch, alen, mapq = int(f[9]), int(f[10]), int(f[11])
            qchrom, rng = qname.rsplit(":", 1)
            qstart, qend = (int(x) for x in rng.split("-"))
            # self-hit: lands back on the interval it was cut from
            if tname == qchrom and ts < qend and qstart < te:
                rows.append((qname, qchrom, qstart, qend, True,
                             tname, ts, te, strand, nmatch / max(1, alen),
                             (qe - qs) / qlen, mapq))
                continue
            rows.append((qname, qchrom, qstart, qend, False,
                         tname, ts, te, strand, nmatch / max(1, alen),
                         (qe - qs) / qlen, mapq))
    return pd.DataFrame(rows, columns=[
        "win", "chrom", "start", "end", "is_self",
        "t_chrom", "t_start", "t_end", "strand", "ident", "qcov", "mapq"])


# Main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--paf", required=True)
    ap.add_argument("--size", type=int, required=True)
    ap.add_argument("--out-windows", required=True)
    ap.add_argument("--out-partners", required=True)
    a = ap.parse_args()

    hits = parse(a.paf, a.size)
    print(f"{len(hits)} PAF records, {hits.win.nunique()} windows with any hit")

    part = hits[(~hits.is_self) & (hits.qcov >= MIN_QCOV)].copy()
    part.to_pickle(a.out_partners)
    print(f"{len(part)} partner hits (self-hit dropped, qcov>={MIN_QCOV})")

    g = part.groupby("win")
    w = pd.DataFrame({
        "n_partners": g.size(),
        "max_ident": g.ident.max(),
        "mean_ident": g.ident.mean(),
        "n_ge95": g.ident.apply(lambda s: (s >= 0.95).sum()),
        "n_ge99": g.ident.apply(lambda s: (s >= 0.99).sum()),
        "n_chrom": g.t_chrom.nunique(),
    })

    # every window that was submitted, including the ones with no partner at all
    allwin = hits[["win", "chrom", "start", "end"]].drop_duplicates("win").set_index("win")
    w = allwin.join(w)
    w["n_partners"] = w.n_partners.fillna(0).astype(int)
    for c in ("n_ge95", "n_ge99", "n_chrom"):
        w[c] = w[c].fillna(0).astype(int)
    # distinguishing positions per kb against the most confusable partner
    w["psv_per_kb"] = (1 - w.max_ident) * 1000
    w["win_size"] = a.size
    w = w.reset_index()
    w.to_pickle(a.out_windows)

    print(f"\nwindows: {len(w)}")
    print(f"  no partner at all      : {(w.n_partners == 0).sum():>8}  "
          f"({100*(w.n_partners==0).mean():.1f}%)")
    print(f"  >=1 partner            : {(w.n_partners >= 1).sum():>8}")
    print(f"  >=1 partner at >=99%   : {(w.n_ge99 >= 1).sum():>8}")
    print("\nmax_ident distribution over windows that have a partner:")
    print(w.loc[w.n_partners > 0, "max_ident"].describe().to_string())
    print("\nn_partners distribution:")
    print(w.n_partners.describe().to_string())


if __name__ == "__main__":
    main()
