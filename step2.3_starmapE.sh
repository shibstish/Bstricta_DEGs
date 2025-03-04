#!/bin/bash
#SBATCH --partition=cas,kamiak
#SBATCH --job-name=star_mapE
#SBATCH --output=/data/lab/busch/stisinai/bstricta/step2_star/star_mapE.out
#SBATCH --error=/data/lab/busch/stisinai/bstricta/step2_star/star_mapE.err
#SBATCH --workdir=/data/lab/busch/stisinai/bstricta/
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1-10:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=

module load star2.7.6a

cd /data/lab/busch/stisinai/bstricta/step1_trimgalore/trimd

#E_7
STAR --genomeDir /data/lab/busch/stisinai/bstricta/step2_star/index_genome --runThreadN 20 --readFilesIn E_7_1_1_val_1.fq.gz E_7_1_2_val_2.fq.gz --readFilesCommand zcat --outSAMstrandField intronMotif --outMulimapperOrder Random --outSAMmultNmax 1 --outFilterMulimapNmax 10000 --outFileNamePrefix /data/lab/busch/stisinai/bstricta/step2_star/mapped/E_7_1

STAR --genomeDir /data/lab/busch/stisinai/bstricta/step2_star/index_genome --runThreadN 20 --readFilesIn E_7_2_1_val_1.fq.gz E_7_2_2_val_2.fq.gz --readFilesCommand zcat --outSAMstrandField intronMotif --outMulimapperOrder Random --outSAMmultNmax 1 --outFilterMulimapNmax 10000 --outFileNamePrefix /data/lab/busch/stisinai/bstricta/step2_star/mapped/E_7_2

STAR --genomeDir /data/lab/busch/stisinai/bstricta/step2_star/index_genome --runThreadN 20 --readFilesIn E_7_3_1_val_1.fq.gz E_7_3_2_val_2.fq.gz --readFilesCommand zcat --outSAMstrandField intronMotif --outMulimapperOrder Random --outSAMmultNmax 1 --outFilterMulimapNmax 10000 --outFileNamePrefix /data/lab/busch/stisinai/bstricta/step2_star/mapped/E_7_3

#E_13
STAR --genomeDir /data/lab/busch/stisinai/bstricta/step2_star/index_genome --runThreadN 20 --readFilesIn E_13_1_1_val_1.fq.gz E_13_1_2_val_2.fq.gz --readFilesCommand zcat --outSAMstrandField intronMotif --outMulimapperOrder Random --outSAMmultNmax 1 --outFilterMulimapNmax 10000 --outFileNamePrefix /data/lab/busch/stisinai/bstricta/step2_star/mapped/E_13_1

STAR --genomeDir /data/lab/busch/stisinai/bstricta/step2_star/index_genome --runThreadN 20 --readFilesIn E_13_2_1_val_1.fq.gz E_13_2_2_val_2.fq.gz --readFilesCommand zcat --outSAMstrandField intronMotif --outMulimapperOrder Random --outSAMmultNmax 1 --outFilterMulimapNmax 10000 --outFileNamePrefix /data/lab/busch/stisinai/bstricta/step2_star/mapped/E_13_2

STAR --genomeDir /data/lab/busch/stisinai/bstricta/step2_star/index_genome --runThreadN 20 --readFilesIn E_13_3_1_val_1.fq.gz E_13_3_2_val_2.fq.gz --readFilesCommand zcat --outSAMstrandField intronMotif --outMulimapperOrder Random --outSAMmultNmax 1 --outFilterMulimapNmax 10000 --outFileNamePrefix /data/lab/busch/stisinai/bstricta/step2_star/mapped/E_13_3
