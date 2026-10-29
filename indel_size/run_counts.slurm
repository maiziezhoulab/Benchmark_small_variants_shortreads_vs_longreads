#!/bin/bash
#SBATCH --job-name=indel_size
#SBATCH --partition=batch
#SBATCH --account=maiziezhou_lab
#SBATCH --cpus-per-task=12
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=/data/maiziezhou_lab/Xiwen/Projects/SNP_INDEL_LR/Platinum-Pedigree/indel_size/logs/counts.out
#SBATCH --error=/data/maiziezhou_lab/Xiwen/Projects/SNP_INDEL_LR/Platinum-Pedigree/indel_size/logs/counts.err
set -e
# count indel sizes for every hap.py VCF in jobs.txt, 12 at a time
B=/data/maiziezhou_lab/Xiwen/Projects/SNP_INDEL_LR/Platinum-Pedigree/indel_size
PY=/home/zhanx88/.conda/envs/radar2/bin/python
while IFS=$'\t' read -r vcf tag; do
  lib=${tag%%|*}; cal=${tag##*|}
  echo "$vcf|$lib|$cal"
done < $B/data/jobs.txt | xargs -P 12 -I{} bash -c '
  IFS="|" read -r vcf lib cal <<< "{}"
  '"$PY"' '"$B"'/scripts/count_bd_by_size.py "$vcf" '"$B"'/data/size_${lib}_${cal}.tsv "${lib}|${cal}"'
echo "all done"
