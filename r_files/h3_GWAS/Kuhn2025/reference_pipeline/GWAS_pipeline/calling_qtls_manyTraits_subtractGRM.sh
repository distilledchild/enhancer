#!/bin/bash

#PBS -S /bin/bash
#PBS -l walltime=1:00:00
#PBS -l nodes=1:ppn=8
#PBS -j oe
#PBS -q condo
#PBS -o log/$PBS_JOBNAME.out
#PBS -e log/$PBS_JOBNAME.err



export R_LIBS=/home/aschitre/R_libs:$R_LIBS
module load R

cd $PBS_O_WORKDIR/


pheno=`head -$PBS_ARRAYID data/phenotypes.txt | tail -1`

Rscript /projects/ps-palmer/apurva/genetic_analysis/code/calling_qtls_manyTraits_subtractGRM.R $pheno $gwas_threshold
