#!/usr/bin/env python3

# Join LR/SR truth calls with the repeat structure around each site.
# One row per GIAB truth variant: whether HiFi and Illumina recovered it, plus
# segdup escape distance / paralogue identity and tandem repeat period / copy number.
# Writes master.pkl.

import numpy as np
import pandas as pd


# Config
A = "/data/maiziezhou_lab/Xiwen/Projects/SNP_INDEL_LR/Platinum-Pedigree/segdup_structure"


# Step 1: Load truth calls
lr = pd.read_csv(f"{A}/data/truth_LR_HifiL1_DV.tsv", sep="\t")
sr = pd.read_csv(f"{A}/data/truth_SR_L1_DV.tsv", sep="\t")

# both hap.py runs share the truth VCF, so row order should match
assert (lr.chrom.values == sr.chrom.values).all() and (lr.pos.values == sr.pos.values).all()

df = lr[["chrom", "pos", "ref", "alt", "bvt"]].copy()
df["LR"] = (lr.bd == "TP").values
df["SR"] = (sr.bd == "TP").values
df["idx"] = np.arange(1, len(df) + 1)          # matches column 4 of sites.bed


# Step 2: Segdup, one row per (site, overlapping UCSC segdup record)
sd = pd.read_csv(f"{A}/data/sites_x_segdup.tsv", sep="\t", header=None,
                 names=["c", "s", "pos0", "idx", "dc", "ds", "de",
                        "alignL", "fracMatch", "fracMatchIndel",
                        "oc", "os", "oe"])
# distance to the nearest edge of the duplicated block
sd["escape"] = np.minimum(sd.pos0 - sd.ds, sd.de - sd.pos0)
sd["blocklen"] = sd.de - sd.ds

g = sd.groupby("idx")
seg = pd.DataFrame({
    "sd_n":        g.size(),                       # paralogue pairs at this locus
    "sd_escape":   g.escape.max(),                 # worst-case walk to unique sequence
    "sd_blocklen": g.blocklen.max(),
    "sd_id_max":   g.fracMatch.max(),
})
# identity of the record that sets the escape distance
hardest = sd.loc[g.escape.idxmax(), ["idx", "fracMatch", "blocklen"]] \
            .rename(columns={"fracMatch": "sd_id_hard", "blocklen": "sd_len_hard"}) \
            .set_index("idx")
seg = seg.join(hardest)


# Step 3: Tandem repeats
tr = pd.read_csv(f"{A}/data/sites_x_trf.tsv", sep="\t", header=None,
                 names=["c", "s", "pos0", "idx", "dc", "ds", "de",
                        "period", "copyNum", "perMatch", "entropy"])
tr["arraylen"] = tr.de - tr.ds
# keep the longest array when several TRF records overlap
tk = tr.loc[tr.groupby("idx").arraylen.idxmax(),
            ["idx", "period", "copyNum", "perMatch", "arraylen"]].set_index("idx")
tk.columns = ["tr_period", "tr_copynum", "tr_permatch", "tr_arraylen"]


# Step 4: Join and write
df = df.set_index("idx").join(seg).join(tk).reset_index()
df["in_segdup"] = df.sd_n.notna()
df["in_tr"] = df.tr_period.notna()

df.to_pickle(f"{A}/data/master.pkl")

n = len(df)
print(f"truth variants: {n:,}")
print(f"  in UCSC segdup: {df.in_segdup.sum():,} ({100*df.in_segdup.mean():.2f}%)")
print(f"  in UCSC TRF:    {df.in_tr.sum():,} ({100*df.in_tr.mean():.2f}%)")
for v in ("SNP", "INDEL"):
    d = df[df.bvt == v]
    ds = d[d.in_segdup]
    print(f"{v}: n={len(d):,}  recall LR={d.LR.mean():.4f} SR={d.SR.mean():.4f} | "
          f"segdup n={len(ds):,} recall LR={ds.LR.mean():.4f} SR={ds.SR.mean():.4f}")
