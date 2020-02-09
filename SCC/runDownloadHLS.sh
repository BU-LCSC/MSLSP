#!/bin/bash -l
#$ -j y


tile=$1
baseDir=$2
imgStartYr=$3
imgEndYr=$4

imgDir="${baseDir}${tile}/images/"


#Get Images from web
##########
for y in `seq $imgStartYr $imgEndYr`;
do
   downloadHLS.sh -t $tile -y $y $imgDir 
done


