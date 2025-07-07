#!/bin/bash
#SBATCH --partition=cas
#SBATCH --job-name=faidx
#SBATCH --output=/data/lab/busch/stisinai/bstricta/index.out
#SBATCH --error=/data/lab/busch/stisinai/bstricta/index.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=70G
#SBATCH --time=10:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=shelby.tisinai@wsu.edu

module load samtools

samtools faidx /data/lab/busch/stisinai/bstricta/reference/Bstricta/assembly/Bstricta_278_v1.fa
