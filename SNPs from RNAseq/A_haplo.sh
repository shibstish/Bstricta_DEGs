#!/bin/bash
#SBATCH --partition=cas
#SBATCH --job-name=A_haplo
#SBATCH --output=/data/lab/busch/stisinai/bstricta/step4_gatk/Ahaplo.out
#SBATCH --error=/data/lab/busch/stisinai/bstricta/step4_gatk/Ahaplo.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=70G
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=shelby.tisinai@wsu.edu

module load gatk/4.1.4.1
module load samtools/1.9

#After putting read groups on the bam files AND verifying that you have created the dictionary file, push the bam files through HaplotypeCaller, outputting gziped vcfs.

#gatk HaplotypeCaller -R /data/lab/busch/stisinai/bstricta/reference/Bstricta/assembly/Bstricta_278_v1.fa -I /data/lab/busch/stisinai/bstricta/step4_gatk/cigar/A_10_RG.bam -O /data/lab/busch/stisinai/bstricta/step4_gatk/haplo3/A_10_vcf.gz -ERC GVCF --sample-name A10

gatk HaplotypeCaller -R /data/lab/busch/stisinai/bstricta/reference/Bstricta/assembly/Bstricta_278_v1.fa -I /data/lab/busch/stisinai/bstricta/step4_gatk/cigar/A_26_RG.bam -O /data/lab/busch/stisinai/bstricta/step4_gatk/haplo3/A_26_vcf.gz -ERC GVCF --sample-name A26

#gatk HaplotypeCaller -R /data/lab/busch/stisinai/bstricta/reference/Bstricta/assembly/Bstricta_278_v1.fa -I /data/lab/busch/stisinai/bstricta/step4_gatk/cigar/A_5_RG.bam -O /data/lab/busch/stisinai/bstricta/step4_gatk/haplo3/A_5_vcf.gz -ERC GVCF --sample-name A5
