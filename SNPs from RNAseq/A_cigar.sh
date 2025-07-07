#!/bin/bash
#SBATCH --partition=cas
#SBATCH --job-name=A_cigar
#SBATCH --output=/data/lab/busch/stisinai/bstricta/step4_gatk/cigar.out
#SBATCH --error=/data/lab/busch/stisinai/bstricta/step4_gatk/cigar.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=70G
#SBATCH --time=10:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=shelby.tisinai@wsu.edu

module load gatk/4.1.8.1

gatk SplitNCigarReads \
      -R /data/lab/busch/stisinai/bstricta/reference/Bstricta/assembly/Bstricta_278_v1.fa \
      -I /data/lab/busch/stisinai/bstricta/step4_gatk/markd/A_10_3_marked.bam \
      -O /data/lab/busch/stisinai/bstricta/step4_gatk/cigar/A_10_3.bam

gatk SplitNCigarReads \
      -R /data/lab/busch/stisinai/bstricta/reference/Bstricta/assembly/Bstricta_278_v1.fa \
      -I /data/lab/busch/stisinai/bstricta/step4_gatk/markd/A_26_2_marked.bam \
      -O /data/lab/busch/stisinai/bstricta/step4_gatk/cigar/A_26_2.bam

gatk SplitNCigarReads \
      -R /data/lab/busch/stisinai/bstricta/reference/Bstricta/assembly/Bstricta_278_v1.fa \
      -I /data/lab/busch/stisinai/bstricta/step4_gatk/markd/A_5_1_marked.bam \
      -O /data/lab/busch/stisinai/bstricta/step4_gatk/cigar/A_5_1.bam
