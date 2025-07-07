#!/bin/bash
#SBATCH --job-name=ReadGroups
#SBATCH --partition=cas    
#SBATCH --mem-per-cpu=40gb           
#SBATCH --mail-type=END,FAIL          
#SBATCH --mail-user=shelby.tisinai@wsu.edu
#SBATCH --ntasks=1                   
#SBATCH --cpus-per-task=1            
#SBATCH --output=/data/lab/busch/stisinai/bstricta/A_RG.out
#SBATCH --error=/data/lab/busch/stisinai/bstricta/A_RG.err


module load java/oracle_1.8.0_92
module load picard/2.21.4
module load samtools/1.9

#Add the read groups. Required ReadGroup flags below. Other flags available: see https://gatk.broadinstitute.org/hc/en-us/articles/360037226472-AddOrReplaceReadGroups-Picard- 
#RGLB = library
#RGPL = sequencing platform (e.g., illumina)
#RGPU = sequencing platform unit/machine (e.g., unit1)
#RGSM = sample name

picard AddOrReplaceReadGroups I=/data/lab/busch/stisinai/bstricta/step4_gatk/cigar/A_10_3.bam O=/data/lab/busch/stisinai/bstricta/step4_gatk/cigar/A_10_RG.bam  RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=A10

#Now, need to re-index the bamfile output from AddOrReplaceReadGroups
samtools index /data/lab/busch/stisinai/bstricta/step4_gatk/cigar/A_10_RG.bam


picard AddOrReplaceReadGroups I=/data/lab/busch/stisinai/bstricta/step4_gatk/cigar/A_26_2.bam O=/data/lab/busch/stisinai/bstricta/step4_gatk/cigar/A_26_RG.bam  RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=A26

#Now, need to re-index the bamfile output from AddOrReplaceReadGroups
samtools index /data/lab/busch/stisinai/bstricta/step4_gatk/cigar/A_26_RG.bam


picard AddOrReplaceReadGroups I=/data/lab/busch/stisinai/bstricta/step4_gatk/cigar/A_5_1.bam O=/data/lab/busch/stisinai/bstricta/step4_gatk/cigar/A_5_RG.bam  RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=A5

#Now, need to re-index the bamfile output from AddOrReplaceReadGroups
samtools index /data/lab/busch/stisinai/bstricta/step4_gatk/cigar/A_5_RG.bam
