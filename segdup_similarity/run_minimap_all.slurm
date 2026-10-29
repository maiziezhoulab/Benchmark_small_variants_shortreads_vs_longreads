#!/bin/bash
#SBATCH --job-name=segdup_realign_all
#SBATCH --partition=batch
#SBATCH --account=maiziezhou_lab
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --export=ALL
#SBATCH --output=/data/maiziezhou_lab/Xiwen/Projects/SNP_INDEL_LR/Platinum-Pedigree/segdup_map/logs/realign_all.out
#SBATCH --error=/data/maiziezhou_lab/Xiwen/Projects/SNP_INDEL_LR/Platinum-Pedigree/segdup_map/logs/realign_all.err
set -e
module load StdEnv/2023
module load minimap2/2.30
BASE=/data/maiziezhou_lab/Xiwen/Projects/SNP_INDEL_LR/Platinum-Pedigree/segdup_map
REF=/data/maiziezhou_lab/Softwares/refdata-GRCh38_no_alt/GRCh38.no_alt_analysis_set.fa

# keep up to 50 secondary hits (the sr preset turns them off by default)
echo "=== 150 (Illumina scale) $(date) ==="
minimap2 -x sr      --secondary=yes -N 50 -p 0.1 -c -t 16 "$REF" "$BASE/data/win_150.fa"   > "$BASE/data/paf_150.paf"
echo "=== 1000 $(date) ==="
minimap2 -x asm10   --secondary=yes -N 50 -p 0.1 -c -t 16 "$REF" "$BASE/data/win_1000.fa"  > "$BASE/data/paf_1000.paf"
echo "=== 10000 (HiFi scale) $(date) ==="
minimap2 -x asm10   --secondary=yes -N 50 -p 0.1 -c -t 16 "$REF" "$BASE/data/win_10000.fa" > "$BASE/data/paf_10000.paf"
wc -l "$BASE"/data/paf_*.paf
echo "done $(date)"
