#!/bin/bash
#SBATCH --job-name=dictionary
#SBATCH --output=/data/lab/busch/stisinai/bstricta/dict.out
#SBATCH --error=/data/lab/busch/stisinai/bstricta/dict.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=70G
#SBATCH --time=10:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=shelby.tisinai@wsu.edu

module load gatk/4.1.8.1

gatk CreateSequenceDictionary -R /data/lab/busch/stisinai/bstricta/reference/Bstricta/assembly/Bstricta_278_v1.fa
