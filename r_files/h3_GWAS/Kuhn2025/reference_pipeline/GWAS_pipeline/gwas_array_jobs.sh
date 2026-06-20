#!/bin/bash

#PBS -N gcta_chr_pheno
#PBS -S /bin/bash
#PBS -l walltime=4:00:00
#PBS -l nodes=1:ppn=5
#PBS -j oe
#PBS -q hotel



#give current working directory to the array script
echo PBS: working directory is $PBS_O_WORKDIR

cd $PBS_O_WORKDIR

#create pheno and chr file and save it in temp directory




#depending on max jobs split the following file into separate files
#then qsub an array job for these file parts
#19 separate array jobs
#submit 1, depend afterok the first one

pheno_jobs=`head -$PBS_ARRAYID data/gwas_jobs/gwas.jobs.$pheno_file.txt | tail -1`

chr=$(echo $pheno_jobs | cut -f1 -d' ')
pheno=$(echo $pheno_jobs | cut -f2 -d' ')




/oasis/tscc/scratch/aschitre/gcta/gcta64 --thread-num 1 --pheno data/pheno/$pheno.txt --bfile genotypes/genotypes  --grm grm/genotypes --mlma-subtract-grm grm/$chr.genotypes --mlma  --out results/gwas/$chr.$pheno
