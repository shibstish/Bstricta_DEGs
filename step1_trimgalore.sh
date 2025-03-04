#!/bin/bash
#SBATCH --partition=cas
#SBATCH --job-name=Trimgalore
#SBATCH --output=/data/lab/busch/stisinai/bstricta/step1_trimgalore/trimd.out
#SBATCH --error=/data/lab/busch/stisinai/bstricta/step1_trimgalore/trimd.err
#SBATCH --workdir=/data/lab/busch/stisinai/bstricta
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=75G
#SBATCH --time=10:00:00
#SBATCH --mail-type=END,FAIL     
#SBATCH --mail-user=

module load fastqc/
module load cutadapt/
module load python/2.7.10


trimg="/data/busch/stisinai/programs/TrimGalore-0.6.6/trim_galore"
pw_data="/data/lab/busch/stisinai/bstricta/RawData"
OUTDIR="/data/lab/busch/stisinai/bstricta/step1_trimgalore"

for i in `cat /data/lab/busch/stisinai/bstricta/samplist`; do
$trimg $pw_data/$i/$i*_1.fq.gz $pw_data/$i/$i*_2.fq.gz --fastqc --paired --illumina --gzip --trim-n --output_dir $OUTDIR --fastqc_args "--nogroup --outdir $OUTDIR" --quality 28
done