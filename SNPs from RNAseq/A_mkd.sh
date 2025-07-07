#!/bin/bash
#SBATCH --partition=cas          
#SBATCH --job-name=A_MarkDupes          
#SBATCH --output=/data/lab/busch/stisinai/bstricta/step4_gatk/markd/A_mkd.out
#SBATCH --error=/data/lab/busch/stisinai/bstricta/step4_gatk/markd/A_mkd.err
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --ntasks-per-node=1
#SBATCH --mail-type=END,FAIL  
#SBATCH --mail-user=shelby.tisinai@wsu.edu

module load java/oracle_1.8.0_92
module load picard/2.21.4

picard MarkDuplicates \
I=/data/lab/busch/stisinai/bstricta/step3_s2b/A_10_3.Aligned.out.sam.sorted.bam \
M=/data/lab/busch/stisinai/bstricta/step3_s2b/A_10_3_metrics.txt \
O=/data/lab/busch/stisinai/bstricta/step4_gatk/markd/A_10_3_marked.bam \
REMOVE_DUPLICATES=true \
VALIDATION_STRINGENCY=SILENT

picard MarkDuplicates \
I=/data/lab/busch/stisinai/bstricta/step3_s2b/A_26_2.Aligned.out.sam.sorted.bam \
M=/data/lab/busch/stisinai/bstricta/step3_s2b/A_26_2_metrics.txt \
O=/data/lab/busch/stisinai/bstricta/step4_gatk/markd/A_26_2_marked.bam \
REMOVE_DUPLICATES=true \
VALIDATION_STRINGENCY=SILENT

picard MarkDuplicates \
I=/data/lab/busch/stisinai/bstricta/step3_s2b/A_5_1.Aligned.out.sam.sorted.bam \
M=/data/lab/busch/stisinai/bstricta/step3_s2b/A_5_1_metrics.txt \
O=/data/lab/busch/stisinai/bstricta/step4_gatk/markd/A_5_1_marked.bam \
REMOVE_DUPLICATES=true \
VALIDATION_STRINGENCY=SILENT
