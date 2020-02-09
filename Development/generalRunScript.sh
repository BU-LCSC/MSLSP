#!/bin/bash -l
#$ -pe omp 28
#$ -l h_rt=12:00:00
#$ -N runTile

#Load required modules
module purge
module load R/3.6.0
module load geos/3.7.0
module load rstudio/1.2.1335
module load qgis/3.4.2
module load gdal/2.3.2

#Run the script
######################
Rscript /usr2/postdoc/dbolt/CodeFolder/git_repos/MSLSP/supporting/run_phenology_for_CDL_anom.r 



