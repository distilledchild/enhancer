#!/bin/bash

#PBS -S /bin/bash
#PBS -l walltime=1:00:00
#PBS -l nodes=1:ppn=3
#PBS -j oe
#PBS -q condo
#PBS -o log/$PBS_JOBNAME.out
#PBS -e log/$PBS_JOBNAME.err


export R_LIBS=/home/aschitre/R_libs:$R_LIBS
module load R

cd $PBS_O_WORKDIR/

Rscript /projects/ps-palmer/apurva/genetic_analysis/code/phewas_gbs.R
