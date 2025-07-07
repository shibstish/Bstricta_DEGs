#!/bin/bash
#SBATCH --job-name=s2b
#SBATCH --partition-cas,kamiak
#SBATCH --output=/data/lab/busch/stisinai/bstricta/step3_sam2bam/s2b.out
#SBATCH --error=/data/lab/busch/stisinai/bstricta/step3_sam2bam/s2b.err
#SBATCH --workdir=/data/lab/busch/stisinai/bstricta/
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1-10:00:00
#SBATCH --mail=END,FAIL
#SBATCH --mail-user=

module load samtools

cd /data/lab/busch/stisinai/bstricta/step2_star/mapped

for i in *.sam
Do
samtools view -S -b $i>/data/lab/busch/stisinai/bstricta/step3_sam2bam/$i.bam

samtools sort /data/lab/busch/stisinai/bstricta/step3_sam2bam/$i.bam -o /data/lab/busch/stisinai/bstricta/step3_sam2bam/$i.sorted.bam
done