#!/bin/bash
#SBATCH --partition=cas,kamiak
#SBATCH --job-name=star_mapD
#SBATCH --output=/data/lab/busch/stisinai/bstricta/step2_star/star_mapD.out
#SBATCH --error=/data/lab/busch/stisinai/bstricta/step2_star/star_mapD.err
#SBATCH --workdir=/data/lab/busch/stisinai/bstricta/
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1-10:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=

module load star2.7.6a

cd /data/lab/busch/stisinai/bstricta/step1_trimgalore/trimd

#D_14
STAR --genomeDir /data/lab/busch/stisinai/bstricta/step2_star/index_genome --runThreadN 20 --readFilesIn D_14_1_1_val_1.fq.gz D_14_1_2_val_2.fq.gz --readFilesCommand zcat --outSAMstrandField intronMotif --outMulimapperOrder Random --outSAMmultNmax 1 --outFilterMulimapNmax 10000 --outFileNamePrefix /data/lab/busch/stisinai/bstricta/step2_star/mapped/D_14_1

STAR --genomeDir /data/lab/busch/stisinai/bstricta/step2_star/index_genome --runThreadN 20 --readFilesIn D_14_2_1_val_1.fq.gz D_14_2_2_val_2.fq.gz --readFilesCommand zcat --outSAMstrandField intronMotif --outMulimapperOrder Random --outSAMmultNmax 1 --outFilterMulimapNmax 10000 --outFileNamePrefix /data/lab/busch/stisinai/bstricta/step2_star/mapped/D_14_2

STAR --genomeDir /data/lab/busch/stisinai/bstricta/step2_star/index_genome --runThreadN 20 --readFilesIn D_14_3_1_val_1.fq.gz D_14_3_2_val_2.fq.gz --readFilesCommand zcat --outSAMstrandField intronMotif --outMulimapperOrder Random --outSAMmultNmax 1 --outFilterMulimapNmax 10000 --outFileNamePrefix /data/lab/busch/stisinai/bstricta/step2_star/mapped/D_14_3

#D_16
STAR --genomeDir /data/lab/busch/stisinai/bstricta/step2_star/index_genome --runThreadN 20 --readFilesIn D_16_1_1_val_1.fq.gz D_16_1_2_val_2.fq.gz --readFilesCommand zcat --outSAMstrandField intronMotif --outMulimapperOrder Random --outSAMmultNmax 1 --outFilterMulimapNmax 10000 --outFileNamePrefix /data/lab/busch/stisinai/bstricta/step2_star/mapped/D_16_1

STAR --genomeDir /data/lab/busch/stisinai/bstricta/step2_star/index_genome --runThreadN 20 --readFilesIn D_16_2_1_val_1.fq.gz D_16_2_2_val_2.fq.gz --readFilesCommand zcat --outSAMstrandField intronMotif --outMulimapperOrder Random --outSAMmultNmax 1 --outFilterMulimapNmax 10000 --outFileNamePrefix /data/lab/busch/stisinai/bstricta/step2_star/mapped/D_16_2

STAR --genomeDir /data/lab/busch/stisinai/bstricta/step2_star/index_genome --runThreadN 20 --readFilesIn D_16_3_1_val_1.fq.gz D_16_3_2_val_2.fq.gz --readFilesCommand zcat --outSAMstrandField intronMotif --outMulimapperOrder Random --outSAMmultNmax 1 --outFilterMulimapNmax 10000 --outFileNamePrefix /data/lab/busch/stisinai/bstricta/step2_star/mapped/D_16_3


#D_20
STAR --genomeDir /data/lab/busch/stisinai/bstricta/step2_star/index_genome --runThreadN 20 --readFilesIn D_20_1_1_val_1.fq.gz D_20_1_2_val_2.fq.gz --readFilesCommand zcat --outSAMstrandField intronMotif --outMulimapperOrder Random --outSAMmultNmax 1 --outFilterMulimapNmax 10000 --outFileNamePrefix /data/lab/busch/stisinai/bstricta/step2_star/mapped/D_20_1

STAR --genomeDir /data/lab/busch/stisinai/bstricta/step2_star/index_genome --runThreadN 20 --readFilesIn D_20_2_1_val_1.fq.gz D_20_2_2_val_2.fq.gz --readFilesCommand zcat --outSAMstrandField intronMotif --outMulimapperOrder Random --outSAMmultNmax 1 --outFilterMulimapNmax 10000 --outFileNamePrefix /data/lab/busch/stisinai/bstricta/step2_star/mapped/D_20_2

STAR --genomeDir /data/lab/busch/stisinai/bstricta/step2_star/index_genome --runThreadN 20 --readFilesIn D_20_3_1_val_1.fq.gz D_20_3_2_val_2.fq.gz --readFilesCommand zcat --outSAMstrandField intronMotif --outMulimapperOrder Random --outSAMmultNmax 1 --outFilterMulimapNmax 10000 --outFileNamePrefix /data/lab/busch/stisinai/bstricta/step2_star/mapped/D_20_3