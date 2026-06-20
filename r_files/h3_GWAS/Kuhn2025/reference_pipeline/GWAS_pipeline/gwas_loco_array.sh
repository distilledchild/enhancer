#!/bin/bash

#PBS -N gwas
#PBS -S /bin/bash
#PBS -l walltime=4:00:00
#PBS -l nodes=1:ppn=5
#PBS -j oe
#PBS -q hotel
#PBS -o log/$PBS_JOBNAME.out
#PBS -e log/$PBS_JOBNAME.err


#give current working directory to the array script
echo PBS: working directory is $PBS_O_WORKDIR

cd $PBS_O_WORKDIR


pheno=`head -$PBS_ARRAYID data/phenotypes.txt | tail -1`



/oasis/tscc/scratch/aschitre/gcta/gcta64 --thread-num 4 --pheno data/pheno/$pheno.txt --bfile genotypes/genotypes --mlma-loco --out loco/$pheno
