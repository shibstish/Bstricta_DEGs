#!/bin/bash
#SBATCH --partition=cas,kamiak
#SBATCH --job-name=star_mapA
#SBATCH --output=/data/lab/busch/stisinai/bstricta/step2_star/star_mapA.out
#SBATCH --error=/data/lab/busch/stisinai/bstricta/step2_star/star_mapA.err
#SBATCH --workdir=/data/lab/busch/stisinai/bstricta/
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1-10:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=

module load star2.7.6a

cd /data/lab/busch/stisinai/bstricta/step1_trimgalore/trimd

#A_10
STAR --genomeDir /data/lab/busch/stisinai/bstricta/step2_star/index_genome --runThreadN 20 --readFilesIn A_10_1_1_val_1.fq.gz A_10_1_2_val_2.fq.gz --readFilesCommand zcat --outSAMstrandField intronMotif --outMulimapperOrder Random --outSAMmultNmax 1 --outFilterMulimapNmax 10000 --outFileNamePrefix /data/lab/busch/stisinai/bstricta/step2_star/mapped/A_10_1

STAR --genomeDir /data/lab/busch/stisinai/bstricta/step2_star/index_genome --runThreadN 20 --readFilesIn A_10_2_1_val_1.fq.gz A_10_2_2_val_2.fq.gz --readFilesCommand zcat --outSAMstrandField intronMotif --outMulimapperOrder Random --outSAMmultNmax 1 --outFilterMulimapNmax 10000 --outFileNamePrefix /data/lab/busch/stisinai/bstricta/step2_star/mapped/A_10_2

STAR --genomeDir /data/lab/busch/stisinai/bstricta/step2_star/index_genome --runThreadN 20 --readFilesIn A_10_3_1_val_1.fq.gz A_10_3_2_val_2.fq.gz --readFilesCommand zcat --outSAMstrandField intronMotif --outMulimapperOrder Random --outSAMmultNmax 1 --outFilterMulimapNmax 10000 --outFileNamePrefix /data/lab/busch/stisinai/bstricta/step2_star/mapped/A_10_3

#A_26
STAR --genomeDir /data/lab/busch/stisinai/bstricta/step2_star/index_genome --runThreadN 20 --readFilesIn A_26_1_1_val_1.fq.gz A_26_1_2_val_2.fq.gz --readFilesCommand zcat --outSAMstrandField intronMotif --outMulimapperOrder Random --outSAMmultNmax 1 --outFilterMulimapNmax 10000 --outFileNamePrefix /data/lab/busch/stisinai/bstricta/step2_star/mapped/A_26_1

STAR --genomeDir /data/lab/busch/stisinai/bstricta/step2_star/index_genome --runThreadN 20 --readFilesIn A_26_2_1_val_1.fq.gz A_26_2_2_val_2.fq.gz --readFilesCommand zcat --outSAMstrandField intronMotif --outMulimapperOrder Random --outSAMmultNmax 1 --outFilterMulimapNmax 10000 --outFileNamePrefix /data/lab/busch/stisinai/bstricta/step2_star/mapped/A_26_2

STAR --genomeDir /data/lab/busch/stisinai/bstricta/step2_star/index_genome --runThreadN 20 --readFilesIn A_26_3_1_val_1.fq.gz A_26_3_2_val_2.fq.gz --readFilesCommand zcat --outSAMstrandField intronMotif --outMulimapperOrder Random --outSAMmultNmax 1 --outFilterMulimapNmax 10000 --outFileNamePrefix /data/lab/busch/stisinai/bstricta/step2_star/mapped/A_26_3


#A_5
STAR --genomeDir /data/lab/busch/stisinai/bstricta/step2_star/index_genome --runThreadN 20 --readFilesIn A_5_1_1_val_1.fq.gz A_5_1_2_val_2.fq.gz --readFilesCommand zcat --outSAMstrandField intronMotif --outMulimapperOrder Random --outSAMmultNmax 1 --outFilterMulimapNmax 10000 --outFileNamePrefix /data/lab/busch/stisinai/bstricta/step2_star/mapped/A_5_1

STAR --genomeDir /data/lab/busch/stisinai/bstricta/step2_star/index_genome --runThreadN 20 --readFilesIn A_5_2_1_val_1.fq.gz A_5_2_2_val_2.fq.gz --readFilesCommand zcat --outSAMstrandField intronMotif --outMulimapperOrder Random --outSAMmultNmax 1 --outFilterMulimapNmax 10000 --outFileNamePrefix /data/lab/busch/stisinai/bstricta/step2_star/mapped/A_5_2

STAR --genomeDir /data/lab/busch/stisinai/bstricta/step2_star/index_genome --runThreadN 20 --readFilesIn A_5_3_1_val_1.fq.gz A_5_3_2_val_2.fq.gz --readFilesCommand zcat --outSAMstrandField intronMotif --outMulimapperOrder Random --outSAMmultNmax 1 --outFilterMulimapNmax 10000 --outFileNamePrefix /data/lab/busch/stisinai/bstricta/step2_star/mapped/A_5_3