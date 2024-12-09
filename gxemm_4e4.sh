#!/bin/bash
#
#$ -S /bin/bash
#$ -o ./Rout  ## where to put standard output (to screen)
#$ -e ./Rout  ## where to put error messages
#$ -cwd            #-- start in current working directory
#$ -l h_data=190G      #### 5G good for 2K samples,
#$ -l time=84:20:00       #-- runtime limit (this requests 24 hours)
#$ -l highp
#$ -t 1:2


export OMP_NUM_THREADS=1 

Rscript gxemm.R 4e4

## module load R/3.6.0; module load plink; qsub -cwd -V gxemm_4e4.sh
