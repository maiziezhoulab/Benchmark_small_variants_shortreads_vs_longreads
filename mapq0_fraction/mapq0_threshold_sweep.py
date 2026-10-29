#!/usr/bin/env python3

# Segdup overlap of SR MAPQ=0 loci and their XA/SA partner locations,
# swept over mapq0_fraction bands. Reads the SR lowmap TSV and SR BAM,
# writes per-partner, per-source and summary TSVs.
#
#   python mapq0_threshold_sweep.py [chrom|all]

import os
import statistics
import sys
from bisect import bisect_right

import pandas as pd
import pysam


# Config

TSV = "/data/maiziezhou_lab/Xiwen/Projects/SNP_INDEL_LR/Platinum-Pedigree/secondary_alignment/NA24385/lowmap_tsv/NA24385_DeepVariant_L1_SR_secalign_datapoint.lowmap.tsv"
BAM = "/lfs/archer.accre.vu/maiziezhou_lab/maiziezhou_lab/Yichen/Projects/SNP_INDEL_LR/ALN_GRCh38_no_alt/SR/SR_L1.bam"
SEGDUP = "/data/maiziezhou_lab/Yichen/Projects/SNP_INDEL_LR/Benchmark_GRCh38/Genome_stratification_v3.6/GRCh38_segdups.bed"
OUTDIR = "/data/maiziezhou_lab/Xiwen/Projects/SNP_INDEL_LR/Platinum-Pedigree/secondary_alignment/NA24385/mapq0_partners/threshold_sweep"
TARGET = sys.argv[1] if len(sys.argv) > 1 else "all"
MAIN = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
BANDS = [(0.25, 0.5), (0.5, 0.75), (0.75, 0.9), (0.9, 1.01)]
GAP = 500   # loci closer than this are merged into one hotspot


# Step 1: Segdup lookup

# merged, sorted intervals per chrom
def load_bed(path):
    iv = {}
    for line in open(path):
        c, s, e = line.split("\t")[:3]
        iv.setdefault(c, []).append((int(s), int(e)))
    idx = {}
    for c, v in iv.items():
        v.sort()
        m = []
        for s, e in v:
            if m and s <= m[-1][1]:
                m[-1][1] = max(m[-1][1], e)
            else:
                m.append([s, e])
        idx[c] = ([x[0] for x in m], m)
    return idx


# is pos0 inside an interval
def hit(idx, c, pos0):
    if c not in idx:
        return False
    st, iv = idx[c]
    i = bisect_right(st, pos0) - 1
    return i >= 0 and iv[i][0] <= pos0 < iv[i][1]


# Step 2: Load loci
os.makedirs(OUTDIR, exist_ok=True)
seg = load_bed(SEGDUP)
chroms = MAIN if TARGET == "all" else [TARGET]
df = pd.read_csv(TSV, sep="\t", usecols=["chrom", "pos", "total_mapped_depth", "mapq0_fraction"])
df = df[(df.total_mapped_depth >= 20) & (df.chrom.isin(chroms))]
bam = pysam.AlignmentFile(BAM, "rb")

# Step 3: Walk each band, collect XA/SA partners of MAPQ=0 reads
long_rows, src_rows, summary = [], [], []
for lo, hi in BANDS:
    band = f"[{lo:.2f},{hi:.2f})"
    bsub = df[(df.mapq0_fraction >= lo) & (df.mapq0_fraction < hi)]
    n_loci, n_reg, s_in, perread = 0, 0, 0, {}
    xa_n, xa_in = 0, 0
    for chrom in chroms:
        sub = sorted(bsub.pos[bsub.chrom == chrom])
        for p in sub:
            src_rows.append((band, chrom, p, int(hit(seg, chrom, p - 1))))
        s_in += sum(hit(seg, chrom, p - 1) for p in sub)
        n_loci += len(sub)
        regions = []
        for pos in sub:
            if regions and pos - regions[-1][1] <= GAP:
                regions[-1][1] = pos
            else:
                regions.append([pos, pos])
        n_reg += len(regions)
        for s, e in regions:
            hot = f"{chrom}:{s}-{e}"
            for r in bam.fetch(chrom, s - 1, e):
                if r.mapping_quality != 0 or r.is_secondary or r.is_supplementary:
                    continue
                for tag, srct in (("XA", "XA"), ("SA", "SA")):
                    if not r.has_tag(tag):
                        continue
                    hs = r.get_tag(tag).rstrip(";").split(";")
                    if srct == "XA":
                        perread[r.query_name] = len(hs)
                    for h in hs:
                        f = h.split(",")
                        c, p = f[0], abs(int(f[1]))
                        nm = f[3] if srct == "XA" else f[5]
                        d_seg = int(hit(seg, c, p - 1))
                        long_rows.append((band, hot, r.query_name, c, p, int(nm), srct, d_seg))
                        if srct == "XA":
                            xa_n += 1
                            xa_in += d_seg
    med = statistics.median(perread.values()) if perread else 0
    summary.append((band, n_loci, n_reg, s_in, n_loci and round(100 * s_in / n_loci, 1),
                    xa_n, xa_in, xa_n and round(100 * xa_in / xa_n, 1), med))
bam.close()

# Step 4: Write out
pd.DataFrame(long_rows, columns=["band", "hotspot", "read", "dest_chrom", "dest_pos", "NM", "source", "dest_in_segdup"]
             ).to_csv(f"{OUTDIR}/sweep_long.tsv", sep="\t", index=False)
pd.DataFrame(src_rows, columns=["band", "chrom", "pos", "src_in_segdup"]
             ).to_csv(f"{OUTDIR}/sweep_sources.tsv", sep="\t", index=False)
sm = pd.DataFrame(summary, columns=["band", "n_loci", "n_hotspots", "src_in_segdup_n", "src_in_segdup_pct",
                                    "n_XA_dests", "dest_in_segdup_n", "dest_in_segdup_pct", "median_XA_per_read"])
sm.to_csv(f"{OUTDIR}/sweep_summary.tsv", sep="\t", index=False)
print(f"scope={TARGET}")
print(sm.to_string(index=False))
print(f"\nwrote -> {OUTDIR}/")
