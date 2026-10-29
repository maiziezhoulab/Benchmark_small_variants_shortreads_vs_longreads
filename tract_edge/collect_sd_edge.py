#!/usr/bin/env python3

# Sample reads at fixed distances from every segdup edge and record their MAPQ.
# For each merged GIAB segdup interval (>= 1 kb), seek at a ladder of distances
# inward and outward from both edges and keep the primary reads covering each point.
# Geometry is left to analyze_sd_edge.py. Writes one gzipped TSV of reads.
#
#   python collect_sd_edge.py --bam x.bam --segdup segdups.bed --out sd_edge_LR.tsv.gz

import argparse
import gzip
from multiprocessing import Pool

import numpy as np
import pysam


# Config
MAIN = set(f"chr{i}" for i in range(1, 23)) | {"chrX", "chrY"}
MIN_LEN = 1000                      # segdups are >= 1 kb by definition; shorter are BED debris
IN_D = [0, 50, 100, 200, 350, 500, 750, 1000, 1500, 2000, 3000, 5000, 7500,
        10000, 15000, 20000, 30000, 50000, 100000]
OUT_D = [50, 100, 200, 350, 500, 750, 1000, 1500, 2000, 3000, 5000, 7500,
         10000, 15000, 20000]
MAX_PER_SEEK = 60
SKIP = 0x4 | 0x100 | 0x200 | 0x400 | 0x800   # unmapped, secondary, qcfail, dup, supplementary


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
        out[c] = [(s, e) for s, e in m if e - s >= MIN_LEN]
    return out


# Step 2: Build seek points

# (pos, interval_start, interval_end, rung), rung > 0 inside, < 0 outside
def seeks_for(chrom, ivs):

    out = []
    for k, (s, e) in enumerate(ivs):
        L = e - s
        prev_e = ivs[k - 1][1] if k > 0 else -10**12
        next_s = ivs[k + 1][0] if k + 1 < len(ivs) else 10**12
        for d in IN_D:
            if d >= L / 2:
                break
            out.append((s + d, s, e, d))
            out.append((e - 1 - d, s, e, d))
        for d in OUT_D:
            p = s - 1 - d                                  # left flank
            if p > prev_e and (s - p) < (p - prev_e):
                out.append((p, s, e, -d))
            p = e + d                                      # right flank
            if p < next_s and (p - e) < (next_s - p):
                out.append((p, s, e, -d))
    return out


# Step 3: Worker function

# fetch primary reads at each seek point, dedup by read name within the chunk
def do_chunk(args):

    chrom, seeks, bam_path, seed = args
    rng = np.random.default_rng(seed)
    bam = pysam.AlignmentFile(bam_path, "rb")
    seen, rows = set(), []
    for pos, s, e, d in seeks:
        try:
            reads = [r for r in bam.fetch(chrom, pos, pos + 1)
                     if not (r.flag & SKIP) and r.reference_end]
        except ValueError:
            continue
        if len(reads) > MAX_PER_SEEK:                       # random, not the leftmost starters
            reads = [reads[i] for i in rng.choice(len(reads), MAX_PER_SEEK, replace=False)]
        for r in reads:
            key = (r.query_name, r.flag & 0xC0)
            if key in seen:
                continue
            seen.add(key)
            ct = r.cigartuples or []
            cl = ct[0][1] if ct and ct[0][0] in (4, 5) else 0
            cr = ct[-1][1] if ct and ct[-1][0] in (4, 5) else 0
            mq = r.get_tag("MQ") if r.has_tag("MQ") else -1
            rows.append((chrom, r.reference_start, r.reference_end, r.mapping_quality,
                         r.flag, mq, r.template_length, cl, cr, s, e, d))
    bam.close()
    return chrom, len(seeks), rows


# Main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--segdup", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--chunks", type=int, default=96)
    a = ap.parse_args()

    ivs = load_merged(a.segdup)
    allseeks = []
    for c in sorted(ivs):
        allseeks.extend((c, *x) for x in seeks_for(c, ivs[c]))
    allseeks.sort(key=lambda x: (x[0], x[1]))
    n_iv = sum(len(v) for v in ivs.values())
    print(f"{n_iv} intervals >= {MIN_LEN} bp, {len(allseeks):,} seek points", flush=True)

    # equal-size chunks that never cross a chromosome
    jobs, size = [], max(1, len(allseeks) // a.chunks)
    cur, cur_c = [], None
    for c, *rest in allseeks:
        if cur and (c != cur_c or len(cur) >= size):
            jobs.append((cur_c, cur, a.bam, len(jobs)))
            cur = []
        cur_c = c
        cur.append(tuple(rest))
    if cur:
        jobs.append((cur_c, cur, a.bam, len(jobs)))

    cols = ["chrom", "rs", "re", "mapq", "flag", "mate_mapq", "tlen",
            "clip_l", "clip_r", "iv_s", "iv_e", "rung"]
    n_rows = done = 0
    with gzip.open(a.out, "wt") as fo, Pool(a.threads) as pool:
        fo.write("\t".join(cols) + "\n")
        for chrom, ns, rows in pool.imap_unordered(do_chunk, jobs):
            for r in rows:
                fo.write("\t".join(map(str, r)) + "\n")
            n_rows += len(rows)
            done += 1
            print(f"  chunk {done}/{len(jobs)} {chrom}: {ns} seeks, {len(rows)} reads "
                  f"(total {n_rows:,})", flush=True)
    print("done", flush=True)


if __name__ == "__main__":
    main()
