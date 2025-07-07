#!/bin/bash
#SBATCH --job-name=angsd_fst
#SBATCH --partition=batch
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=150gb
#SBATCH --time=100:00:00
#SBATCH --output=/scratch/dd66718/logs/Log.%j
#SBATCH --mail-user=dd66718@uga.edu
#SBATCH --mail-type=END,FAIL

### This script was prepared by Derek Denney, Sept 2024
### This uses ANGSD to generate saf files, which are then used to estimate Fst for each Boechera stricta garden population
### Custom R scripts were used to summarize global Fst and per-locus Fst
### See https://www.popgen.dk/angsd/index.php/Fst for more information

INDIR='/scratch/dd66718/gwas_gea/dupe_filtered/'
OUTDIR='/scratch/dd66718/diversity/angsd/results'
REF='/work/jtalab/Bstr_genome/SAD12v2.2_withgap.fasta'
cores='4'

ml angsd/0.940

### Generate SAF files using ANGSD
cd ${INDIR}
while read garden
do
	angsd -bam $INDIR/${garden}.bamlist \
	-ref ${REF} \
	-anc ${REF} \
	-out ${OUTDIR}/${garden} \
	-minQ 20 \
	-dosaf 1 \
	-doCounts 1 \
	-GL 1 \
	-baq 2 \
	-doGLF 2 \
	-minMapQ 30 \
	-doMajorMinor 4 \
	-doPost 1 \
	-doGeno 2 \
	-P ${cores} \
	-doMaf 1
done < gardens.list

cores='4'
cd $OUTDIR

### Use realSFS to calculate folded SFS from saf files

for x in *.saf.idx
do
	realSFS -seed 1234 -P ${cores} -maxIter 1000000 -fold 1 ${x} > ${x%%.saf.idx}.sfs
	realSFS saf2theta -seed 1234 -P ${cores} -sfs ${x%%.saf.idx}.sfs ${x} -outname ${x%%.saf.idx} -fold 1
	thetaStat do_stat ${x%%.saf.idx}.thetas.idx
done

### Estimate Fst from SFS for each pairwise population

mkdir -p ${OUTDIR}/FST

realSFS fst index 250.saf.idx 273.saf.idx -sfs FST/NB_250-273.ml -fstout FST/250_273 -whichFst 1
realSFS fst index 255.saf.idx 270.saf.idx -sfs FST/NB_270-255.ml -fstout FST/255_270 -whichFst 1
realSFS fst index 255.saf.idx 273.saf.idx -sfs FST/NB_273-255.ml -fstout FST/255_273 -whichFst 1
realSFS fst index 255.saf.idx 273.saf.idx -sfs FST/NB_273-255.ml -fstout FST/255_273 -whichFst 1
realSFS fst index 255.saf.idx 283.saf.idx -sfs FST/NB_283-255.ml -fstout FST/255_283 -whichFst 1
realSFS fst index 270.saf.idx 273.saf.idx -sfs FST/NB_273-270.ml -fstout FST/270_273 -whichFst 1
realSFS fst index 270.saf.idx 283.saf.idx -sfs FST/NB_283-270.ml -fstout FST/270_283 -whichFst 1
realSFS fst index 273.saf.idx 283.saf.idx -sfs FST/NB_283-273.ml -fstout FST/273_283 -whichFst 1

cd ${OUTDIR}/FST

### Summarize Fst using Rscript
ml R
for x in *fst.idx
do
	realSFS fst print ${x} > ${x%%.sfs.idx}.print
	Rscript fst_print.R ${x%%.sfs.idx}
done

Rscript compile_Fst.R
