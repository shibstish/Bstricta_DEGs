#!/bin/bash
#SBATCH --partition=cas
#SBATCH --job-name=gff2gtf
#SBATCH --output=/data/lab/busch/stisinai/log_folder/out/gff2gtf.out
#SBATCH --error=/data/lab/busch/stisinai/log_folder/err/gff2gtf.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=70G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=

gffread="/data/lab/busch/stisinai/programs/gffread-0.12.7/gffread/gffread"

$gffread -E /data/lab/busch/stisinai/ref_genome/TAIR10/Arabidopsis_thaliana.TAIR10.59.gff3 -T -o /data/lab/busch/stisinai/ref_genome/TAIR10/Arabidopsis_thaliana.TAIR10.gtf | more 