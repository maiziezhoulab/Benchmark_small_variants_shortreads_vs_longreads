#!/usr/bin/env python3

# Query-side TP/FP calls from a hap.py VCF, used to compute precision per stratum.
# Writes chrom, pos, bvt, bd, one row per query call with BD in {TP, FP}.
#
#   python extract_query_bd.py happy.vcf.gz out.tsv

import gzip
import sys

vcf_path, out_path = sys.argv[1], sys.argv[2]
n = {"TP": 0, "FP": 0}
with gzip.open(vcf_path, "rt") as fh, open(out_path, "w") as out:
    out.write("chrom\tpos\tbvt\tbd\n")
    for line in fh:
        if line[0] == "#":
            if line.startswith("#CHROM"):
                i_q = line.rstrip("\n").split("\t").index("QUERY")
            continue
        f = line.rstrip("\n").split("\t")
        d = dict(zip(f[8].split(":"), f[i_q].split(":")))
        bd = d.get("BD")
        if bd not in n:
            continue
        out.write(f"{f[0]}\t{f[1]}\t{d.get('BVT', '.')}\t{bd}\n")
        n[bd] += 1
print(f"{out_path}: TP={n['TP']:,} FP={n['FP']:,}", flush=True)
