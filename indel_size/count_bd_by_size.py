#!/usr/bin/env python3

# Count hap.py TP/FN/FP INDELs by length, for one caller.
# Truth side (TP/FN) gives recall, query side (TP/FP) gives precision.
# Deletions are negative, insertions positive. Writes one row per length.
#
#   python count_bd_by_size.py hap.vcf.gz out.tsv "LIB|CALLER"

import gzip
import sys
from collections import defaultdict


# Config
vcf_path, out_path, tag = sys.argv[1], sys.argv[2], sys.argv[3]


# Step 1: Indel length from the sample's GT, longest change wins
def size_of(gt, ref, alts):
    best = 0
    for a in gt.replace("|", "/").split("/"):
        if a in (".", "0"):
            continue
        try:
            alt = alts[int(a) - 1]
        except (ValueError, IndexError):
            continue
        d = len(alt) - len(ref)
        if abs(d) > abs(best):
            best = d
    return best


# Step 2: Count BD labels per length
cnt = defaultdict(lambda: defaultdict(int))
n_rows = 0
with gzip.open(vcf_path, "rt") as fh:
    for line in fh:
        if line[0] == "#":
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                i_t, i_q = cols.index("TRUTH"), cols.index("QUERY")
            continue
        f = line.rstrip("\n").split("\t")
        ref, alts, fmt = f[3], f[4].split(","), f[8].split(":")
        n_rows += 1
        for i, side in ((i_t, "truth"), (i_q, "query")):
            d = dict(zip(fmt, f[i].split(":")))
            bd = d.get("BD")
            if bd not in ("TP", "FN", "FP") or d.get("BVT") != "INDEL":
                continue
            if side == "truth" and bd == "FP":      # FP only exists on the query side
                continue
            if side == "query" and bd == "FN":
                continue
            s = size_of(d.get("GT", "."), ref, alts)
            cnt[s][f"{side}_{bd}"] += 1

# Step 3: Write counts
with open(out_path, "w") as out:
    out.write("caller\tsize\ttruth_TP\ttruth_FN\tquery_TP\tquery_FP\n")
    for s in sorted(cnt):
        c = cnt[s]
        out.write(f"{tag}\t{s}\t{c['truth_TP']}\t{c['truth_FN']}\t"
                  f"{c['query_TP']}\t{c['query_FP']}\n")
tot = {k: sum(c.get(k, 0) for c in cnt.values())
       for k in ("truth_TP", "truth_FN", "query_TP", "query_FP")}
print(f"{tag}: rows={n_rows:,} " + " ".join(f"{k}={v:,}" for k, v in tot.items()), flush=True)
