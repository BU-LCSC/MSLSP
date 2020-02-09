#!/bin/bash -l
#$ -j y

#SCC
#Load packages
#######
#Load required modules
module purge
module load R/3.6.0
module load rstudio/1.2.1335
module load python2/2.7.16
module load awscli

tile=$1
outDir=$2
imgStartYr=$3
imgEndYr=$4

Rscript download_S10.r $tile $outDir $imgStartYr $imgEndYr &
wait

