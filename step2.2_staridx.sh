#!/bin/bash
#SBATCH --partition=cas,kamiak
#SBATCH --job-name=star_idx
#SBATCH --output=/data/lab/busch/stisinai/bstricta/step2_star/star_idx.out
#SBATCH --error=/data/lab/busch/stisinai/bstricta/step2_star/star_idx.err
#SBATCH --workdir=/data/lab/busch/stisinai/bstricta/
#SBATCH --mem=75G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1-10:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=

module load star/2.7.6a

mkdir -p /data/lab/busch/stisinai/bstricta/step2_star_index_genome

STAR --runThreadN 16 --runMode genomeGenerate --genomeDir /data/lab/busch/stisinai/bstricta/step2_star/index_genome/ --genomeFastaFiles /data/lab/busch/stisinai/ref_genome/TAIR10/aThaliana.fa --sjdbGTFfile /data/lab/busch/stisinai/ref_genome/TAIR10/Arabidopsis_thaliana.TAIR10.gtf

echo "index genome generated"