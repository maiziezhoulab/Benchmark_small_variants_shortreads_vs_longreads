#!/bin/bash
#SBATCH --job-name=lmedge_LR
#SBATCH --partition=batch
#SBATCH --account=maiziezhou_lab
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/data/maiziezhou_lab/Xiwen/Projects/SNP_INDEL_LR/Platinum-Pedigree/segdup_map/logs/lowmap_edge_LR.out
#SBATCH --error=/data/maiziezhou_lab/Xiwen/Projects/SNP_INDEL_LR/Platinum-Pedigree/segdup_map/logs/lowmap_edge_LR.err
set -e
# /lfs is a FUSE mount that can lag on a fresh compute node
[ -r "/lfs/archer.accre.vu/maiziezhou_lab/maiziezhou_lab/Yichen/Projects/SNP_INDEL_LR/ALN_GRCh38_no_alt/Hifi/Hifi_L1.bam" ] || { scontrol requeue $SLURM_JOB_ID; exit 1; }
/home/zhanx88/.conda/envs/radar2/bin/python /data/maiziezhou_lab/Xiwen/Projects/SNP_INDEL_LR/Platinum-Pedigree/segdup_map/scripts/collect_sd_edge.py \
  --bam /lfs/archer.accre.vu/maiziezhou_lab/maiziezhou_lab/Yichen/Projects/SNP_INDEL_LR/ALN_GRCh38_no_alt/Hifi/Hifi_L1.bam \
  --segdup /data/maiziezhou_lab/Yichen/Projects/SNP_INDEL_LR/Benchmark_GRCh38/Genome_stratification_v3.6/GRCh38_lowmappabilityall.bed \
  --out /data/maiziezhou_lab/Xiwen/Projects/SNP_INDEL_LR/Platinum-Pedigree/segdup_map/data/lowmap_edge_LR.tsv.gz \
  --threads 16
