#!/usr/bin/env python3

# Cut the GIAB segdup intervals into fixed-size windows and write them as FASTA.
# Window size stands in for read length: 150 bp ~ Illumina, 10 kb ~ HiFi.
# Read names are >chrom:start-end (0-based) so the self-hit can be dropped later.
#   python make_windows.py --bed segdups.bed --ref GRCh38.fa --size 1000 --out win_1000.fa

import argparse

import pysam


# Config
MAIN = set(f"chr{i}" for i in range(1, 23)) | {"chrX", "chrY"}


# Main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bed", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--size", type=int, required=True)
    ap.add_argument("--step", type=int, default=0, help="0 = non-overlapping")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    step = a.step or a.size

    fa = pysam.FastaFile(a.ref)
    n_win = n_skip = 0
    with open(a.out, "w") as out:
        for line in open(a.bed):
            if line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t")
            chrom, s, e = f[0], int(f[1]), int(f[2])
            if chrom not in MAIN:
                continue
            for ws in range(s, e - a.size + 1, step):
                we = ws + a.size
                seq = fa.fetch(chrom, ws, we).upper()
                # a window that is mostly N tells us nothing about paralogy
                if seq.count("N") > 0.1 * a.size:
                    n_skip += 1
                    continue
                out.write(f">{chrom}:{ws}-{we}\n{seq}\n")
                n_win += 1
    print(f"size={a.size} step={step}  windows={n_win}  skipped_N={n_skip}")


if __name__ == "__main__":
    main()
