#!/bin/bash
#SBATCH --partition=cas
#SBATCH --job-name=GenotypeGVCF
#SBATCH --output=/data/lab/busch/stisinai/bstricta/geno.out
#SBATCH --error=/data/lab/busch/stisinai/bstricta/geno.err
#SBATCH --workdir=/data/lab/busch/stisinai/bstricta/step4_gatk
#SBATCH --nodes=1
#SBATCH --mem=64G
#SBATCH --time=70:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=shelby.tisinai@wsu.edu

module load gatk/4.1.4.1
module load samtools/1.9

#This script is to be used after merging the haplotyped sample vcfs together to form one vcf for each population
#genotype the population

gatk GenotypeGVCFs \
   -R /data/lab/busch/stisinai/bstricta/reference/Bstricta/assembly/Bstricta_278_v1.fa \
   -V /data/lab/busch/stisinai/bstricta/step4_gatk/merged/bstricta_F_vcf \
   -O /data/lab/busch/stisinai/bstricta/step4_gatk/merged/bstricta_G.vcf.gz
