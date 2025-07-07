#!/bin/bash
#SBATCH --partition=cas
#SBATCH --job-name=HardFilter
#SBATCH --partition=cas
#SBATCH --output=/data/lab/busch/stisinai/bstricta/filter.out
#SBATCH --error=/data/lab/busch/stisinai/bstricta/filter.err
#SBATCH --nodes=1
#SBATCH --mem=64G
#SBATCH --time=70:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=shelby.tisinai@wsu.edu

module load gatk/4.2.6.1
module load samtools/1.9

#after genotyping the merged gvcf file, apply hard filters. Anything starting with the flag "-filter" is an industry standard flag. The "--genotype-filter-expression/-genotype-filter-name" and "--set-filtered-genotype-to-no-call" are additional filters used to filter on the format fields (e.g., GT, dp, ect.) of the vcf, instead of the usual "filter" field

#gatk VariantFiltration \
#    	-R /data/lab/busch/stisinai/bstricta/reference/Bstricta/assembly/Bstricta_278_v1.fa \
#    	-V /data/lab/busch/stisinai/bstricta/step4_gatk/merged/bstricta_vcf.gz \
#    	-O /data/lab/busch/stisinai/bstricta/step4_gatk/merged/bstricta_F_vcf.gz \
#    	-filter-expresson "QD < 2.0" --filter-name "QD2" \
#    	-filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
#    	-filter-expression "SOR > 3.0" --filter-name "SOR3" \
#    	-filter-expression "FS > 60.0" --filter-name "FS60" \
#	-filter-expression "MQ < 50.0" --filter-name "MQ50" \
#    	-filter-expression "MAPQ=255" --filter-name "MQ255" \
#    	-filter-expression "MQRankSum < -2.5" --filter-name "MQRankSum-2.5" \
#    	-filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
#    	-filter-expression "ExcessHet > 11" --filter-name "ExcessHet" \
#    	--genotype-filter-expression "DP < 10" --genotype-filter-name "DP10" \
#    	--set-filtered-genotype-to-no-call true

gatk VariantFiltration \
	-R /data/lab/busch/stisinai/bstricta/reference/Bstricta/assembly/Bstricta_278_v1.fa \
        -V /data/lab/busch/stisinai/bstricta/step4_gatk/merged/bstricta_vcf.gz \
	--window 35 \
	--cluster 3 \
	--filter-name "FS" \
	--filter "FS > 30.0" \
	--filter-name "QD" \
	--filter "QD < 2.0" \
	-O /data/lab/busch/stisinai/bstricta/step4_gatk/merged/bstricta_F_vcf.gz
