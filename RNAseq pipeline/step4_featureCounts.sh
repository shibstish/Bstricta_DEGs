#!/bin/bash
#SBATCH --job-name=feat_counts
#SBATCH --partition-cas,kamiak
#SBATCH --output=/data/lab/busch/stisinai/bstricta/step4_featureCounts/fc.out
#SBATCH --error=/data/lab/busch/stisinai/bstricta/step4_featureCounts/fc.out
#SBATCH --workdir=/data/lab/busch/stisinai/bstricta/
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1-10:00:00
#SBATCH --mail=END,FAIL
#SBATCH --mail-user=

# load condo environment and activate featureCounts

module load anaconda3
source activate featureCounts_env

featureCounts /data/lab/busch/stisinai/bstricta/step3_sam2bam/*.sorted.bam -p -O -T 16 -a /data/lab/busch/stisinai/ref_genome/TAIR10/Arabidopsis_thaliana.TAIR10.gtf -o /data/lab/busch/stisinai/step4_featureCounts/featCounts_mx.txt