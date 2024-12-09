#!/bin/bash
#
#$ -S /bin/bash
#$ -o ./Rout  ## where to put standard output (to screen)
#$ -e ./Rout  ## where to put error messages
#$ -cwd            #-- start in current working directory
#$ -l h_data=5G
#$ -l time=18:20:00       #-- runtime limit (this requests 24 hours)
#$ -l highp
#$ -t 1:2950

export OMP_NUM_THREADS=1 

Rscript gwas.R ${SGE_TASK_ID} 
Rscript gxewas.R ${SGE_TASK_ID} 

## module load R/3.6.0; module load plink; qsub -cwd -V gxewas.sh
