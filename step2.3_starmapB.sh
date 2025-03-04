#!/bin/bash
#SBATCH --partition=cas,kamiak
#SBATCH --job-name=star_mapB
#SBATCH --output=/data/lab/busch/stisinai/bstricta/step2_star/star_mapB.out
#SBATCH --error=/data/lab/busch/stisinai/bstricta/step2_star/star_mapB.err
#SBATCH --workdir=/data/lab/busch/stisinai/bstricta/
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1-10:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=

module load star2.7.6a

cd /data/lab/busch/stisinai/bstricta/step1_trimgalore/trimd

#B_19
STAR --genomeDir /data/lab/busch/stisinai/bstricta/step2_star/index_genome --runThreadN 20 --readFilesIn B_19_1_1_val_1.fq.gz B_19_1_2_val_2.fq.gz --readFilesCommand zcat --outSAMstrandField intronMotif --outMulimapperOrder Random --outSAMmultNmax 1 --outFilterMulimapNmax 10000 --outFileNamePrefix /data/lab/busch/stisinai/bstricta/step2_star/mapped/B_19_1

STAR --genomeDir /data/lab/busch/stisinai/bstricta/step2_star/index_genome --runThreadN 20 --readFilesIn B_19_2_1_val_1.fq.gz B_19_2_2_val_2.fq.gz --readFilesCommand zcat --outSAMstrandField intronMotif --outMulimapperOrder Random --outSAMmultNmax 1 --outFilterMulimapNmax 10000 --outFileNamePrefix /data/lab/busch/stisinai/bstricta/step2_star/mapped/B_19_2

STAR --genomeDir /data/lab/busch/stisinai/bstricta/step2_star/index_genome --runThreadN 20 --readFilesIn B_19_3_1_val_1.fq.gz B_19_3_2_val_2.fq.gz --readFilesCommand zcat --outSAMstrandField intronMotif --outMulimapperOrder Random --outSAMmultNmax 1 --outFilterMulimapNmax 10000 --outFileNamePrefix /data/lab/busch/stisinai/bstricta/step2_star/mapped/B_19_3

#B_28
STAR --genomeDir /data/lab/busch/stisinai/bstricta/step2_star/index_genome --runThreadN 20 --readFilesIn B_28_1_1_val_1.fq.gz B_28_1_2_val_2.fq.gz --readFilesCommand zcat --outSAMstrandField intronMotif --outMulimapperOrder Random --outSAMmultNmax 1 --outFilterMulimapNmax 10000 --outFileNamePrefix /data/lab/busch/stisinai/bstricta/step2_star/mapped/B_28_1

STAR --genomeDir /data/lab/busch/stisinai/bstricta/step2_star/index_genome --runThreadN 20 --readFilesIn B_28_2_1_val_1.fq.gz B_28_2_2_val_2.fq.gz --readFilesCommand zcat --outSAMstrandField intronMotif --outMulimapperOrder Random --outSAMmultNmax 1 --outFilterMulimapNmax 10000 --outFileNamePrefix /data/lab/busch/stisinai/bstricta/step2_star/mapped/B_28_2

STAR --genomeDir /data/lab/busch/stisinai/bstricta/step2_star/index_genome --runThreadN 20 --readFilesIn B_28_3_1_val_1.fq.gz B_28_3_2_val_2.fq.gz --readFilesCommand zcat --outSAMstrandField intronMotif --outMulimapperOrder Random --outSAMmultNmax 1 --outFilterMulimapNmax 10000 --outFileNamePrefix /data/lab/busch/stisinai/bstricta/step2_star/mapped/B_28_3


#B_30
STAR --genomeDir /data/lab/busch/stisinai/bstricta/step2_star/index_genome --runThreadN 20 --readFilesIn B_30_1_1_val_1.fq.gz B_30_1_2_val_2.fq.gz --readFilesCommand zcat --outSAMstrandField intronMotif --outMulimapperOrder Random --outSAMmultNmax 1 --outFilterMulimapNmax 10000 --outFileNamePrefix /data/lab/busch/stisinai/bstricta/step2_star/mapped/B_30_1

STAR --genomeDir /data/lab/busch/stisinai/bstricta/step2_star/index_genome --runThreadN 20 --readFilesIn B_30_2_1_val_1.fq.gz B_30_2_2_val_2.fq.gz --readFilesCommand zcat --outSAMstrandField intronMotif --outMulimapperOrder Random --outSAMmultNmax 1 --outFilterMulimapNmax 10000 --outFileNamePrefix /data/lab/busch/stisinai/bstricta/step2_star/mapped/B_30_2

STAR --genomeDir /data/lab/busch/stisinai/bstricta/step2_star/index_genome --runThreadN 20 --readFilesIn B_30_3_1_val_1.fq.gz B_30_3_2_val_2.fq.gz --readFilesCommand zcat --outSAMstrandField intronMotif --outMulimapperOrder Random --outSAMmultNmax 1 --outFilterMulimapNmax 10000 --outFileNamePrefix /data/lab/busch/stisinai/bstricta/step2_star/mapped/B_30_3