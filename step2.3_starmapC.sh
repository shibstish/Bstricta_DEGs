#!/bin/bash
#SBATCH --partition=cas,kamiak
#SBATCH --job-name=star_mapC
#SBATCH --output=/data/lab/busch/stisinai/bstricta/step2_star/star_mapC.out
#SBATCH --error=/data/lab/busch/stisinai/bstricta/step2_star/star_mapC.err
#SBATCH --workdir=/data/lab/busch/stisinai/bstricta/
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1-10:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=

module load star2.7.6a

cd /data/lab/busch/stisinai/bstricta/step1_trimgalore/trimd

#C_10
STAR --genomeDir /data/lab/busch/stisinai/bstricta/step2_star/index_genome --runThreadN 20 --readFilesIn C_10_1_1_val_1.fq.gz C_10_1_2_val_2.fq.gz --readFilesCommand zcat --outSAMstrandField intronMotif --outMulimapperOrder Random --outSAMmultNmax 1 --outFilterMulimapNmax 10000 --outFileNamePrefix /data/lab/busch/stisinai/bstricta/step2_star/mapped/C_10_1

STAR --genomeDir /data/lab/busch/stisinai/bstricta/step2_star/index_genome --runThreadN 20 --readFilesIn C_10_2_1_val_1.fq.gz C_10_2_2_val_2.fq.gz --readFilesCommand zcat --outSAMstrandField intronMotif --outMulimapperOrder Random --outSAMmultNmax 1 --outFilterMulimapNmax 10000 --outFileNamePrefix /data/lab/busch/stisinai/bstricta/step2_star/mapped/C_10_2

STAR --genomeDir /data/lab/busch/stisinai/bstricta/step2_star/index_genome --runThreadN 20 --readFilesIn C_10_3_1_val_1.fq.gz C_10_3_2_val_2.fq.gz --readFilesCommand zcat --outSAMstrandField intronMotif --outMulimapperOrder Random --outSAMmultNmax 1 --outFilterMulimapNmax 10000 --outFileNamePrefix /data/lab/busch/stisinai/bstricta/step2_star/mapped/C_10_3

#C_11
STAR --genomeDir /data/lab/busch/stisinai/bstricta/step2_star/index_genome --runThreadN 20 --readFilesIn C_11_1_1_val_1.fq.gz C_11_1_2_val_2.fq.gz --readFilesCommand zcat --outSAMstrandField intronMotif --outMulimapperOrder Random --outSAMmultNmax 1 --outFilterMulimapNmax 10000 --outFileNamePrefix /data/lab/busch/stisinai/bstricta/step2_star/mapped/C_11_1

STAR --genomeDir /data/lab/busch/stisinai/bstricta/step2_star/index_genome --runThreadN 20 --readFilesIn C_11_2_1_val_1.fq.gz C_11_2_2_val_2.fq.gz --readFilesCommand zcat --outSAMstrandField intronMotif --outMulimapperOrder Random --outSAMmultNmax 1 --outFilterMulimapNmax 10000 --outFileNamePrefix /data/lab/busch/stisinai/bstricta/step2_star/mapped/C_11_2

STAR --genomeDir /data/lab/busch/stisinai/bstricta/step2_star/index_genome --runThreadN 20 --readFilesIn C_11_3_1_val_1.fq.gz C_11_3_2_val_2.fq.gz --readFilesCommand zcat --outSAMstrandField intronMotif --outMulimapperOrder Random --outSAMmultNmax 1 --outFilterMulimapNmax 10000 --outFileNamePrefix /data/lab/busch/stisinai/bstricta/step2_star/mapped/C_11_3


#C_22
STAR --genomeDir /data/lab/busch/stisinai/bstricta/step2_star/index_genome --runThreadN 20 --readFilesIn C_22_1_1_val_1.fq.gz C_22_1_2_val_2.fq.gz --readFilesCommand zcat --outSAMstrandField intronMotif --outMulimapperOrder Random --outSAMmultNmax 1 --outFilterMulimapNmax 10000 --outFileNamePrefix /data/lab/busch/stisinai/bstricta/step2_star/mapped/C_22_1

STAR --genomeDir /data/lab/busch/stisinai/bstricta/step2_star/index_genome --runThreadN 20 --readFilesIn C_22_2_1_val_1.fq.gz C_22_2_2_val_2.fq.gz --readFilesCommand zcat --outSAMstrandField intronMotif --outMulimapperOrder Random --outSAMmultNmax 1 --outFilterMulimapNmax 10000 --outFileNamePrefix /data/lab/busch/stisinai/bstricta/step2_star/mapped/C_22_2

STAR --genomeDir /data/lab/busch/stisinai/bstricta/step2_star/index_genome --runThreadN 20 --readFilesInC_22_3_1_val_1.fq.gz C_22_3_2_val_2.fq.gz --readFilesCommand zcat --outSAMstrandField intronMotif --outMulimapperOrder Random --outSAMmultNmax 1 --outFilterMulimapNmax 10000 --outFileNamePrefix /data/lab/busch/stisinai/bstricta/step2_star/mapped/C_22_3