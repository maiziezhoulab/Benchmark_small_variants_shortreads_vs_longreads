#!/usr/bin/env python3

# Per-truth-site TP/FN from a hap.py annotated VCF.
# A truth variant has BD=TP (found by the query) or BD=FN (missed); BVT is SNP or INDEL.
# Writes chrom, pos, ref, alt, bvt, bd, one row per truth variant.
#
#   python extract_truth_bd.py happy.vcf.gz out.tsv

import gzip
import sys

vcf_path, out_path = sys.argv[1], sys.argv[2]

n_out = 0
with gzip.open(vcf_path, "rt") as fh, open(out_path, "w") as out:
    out.write("chrom\tpos\tref\talt\tbvt\tbd\n")
    for line in fh:
        if line[0] == "#":
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                # find the TRUTH column instead of assuming its position
                i_truth = cols.index("TRUTH")
            continue
        f = line.rstrip("\n").split("\t")
        fmt = f[8].split(":")
        truth = f[i_truth].split(":")
        d = dict(zip(fmt, truth))
        bd = d.get("BD", ".")
        if bd not in ("TP", "FN"):        # '.' = no truth variant here (FP/UNK rows)
            continue
        out.write(f"{f[0]}\t{f[1]}\t{f[3]}\t{f[4]}\t{d.get('BVT','.')}\t{bd}\n")
        n_out += 1

print(f"{vcf_path} -> {out_path}: {n_out:,} truth variants", flush=True)
