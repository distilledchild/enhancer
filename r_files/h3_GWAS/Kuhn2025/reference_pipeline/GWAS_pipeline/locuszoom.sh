#!/bin/bash

#PBS -S /bin/bash
#PBS -l walltime=6:00:00
#PBS -l nodes=1:ppn=3
#PBS -j oe
#PBS -q condo

export R_LIBS=/home/aschitre/R_libs:$R_LIBS
module load R
module load python

cd $PBS_O_WORKDIR/

Rscript /projects/ps-palmer/apurva/genetic_analysis/code/locuszoom.R
