#!/bin/bash
#
#$ -S /bin/bash
#$ -o ./Rout  ## where to put standard output (to screen)
#$ -e ./Rout  ## where to put error messages
#$ -cwd            #-- start in current working directory
#$ -l h_data=4G
#$ -l time=18:20:00       #-- runtime limit (this requests 24 hours)
#$ -l highp
#$ -t 1:10

export OMP_NUM_THREADS=1 

Rscript gwas.R -11
Rscript gwas.R -12
Rscript gwas.R -13
Rscript gwas.R -14
Rscript gwas.R -15
Rscript gwas.R -16
Rscript gwas.R -17
Rscript gxewas.R -11
Rscript gxewas.R -12
Rscript gxewas.R -13
Rscript gxewas.R -14
Rscript gxewas.R -15
Rscript gxewas.R -16
Rscript gxewas.R -17

## module load R/3.6.0; module load plink; qsub -cwd -V targeted_gxewas.sh
