#!/bin/bash
#
#$ -S /bin/bash
#$ -o ./Rout  ## where to put standard output (to screen)
#$ -e ./Rout  ## where to put error messages
#$ -cwd            #-- start in current working directory
#$ -l h_data=10G      #### 5G good for 2K samples,
#$ -l time=22:20:00       #-- runtime limit (this requests 24 hours)
#$ -l highp
#$ -t 1:20

export OMP_NUM_THREADS=1 

Rscript ce.R

## module load R/3.6.1; module load plink; qsub -cwd -V ce.sh 
