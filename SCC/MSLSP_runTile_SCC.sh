#!/bin/bash -l
#$ -j y

#SCC
#Load packages
#######
#Load required modules
module purge
module load jq
module load python3/3.7.7
module load gdal
module load grass
module load sqlite3/3.37.2
module load R/4.2.1
module load rstudio/2022.07.2-576


tile=$1
parameters=$2
timeStamp=$3


rScript=$( jq --raw-output .SCC.rScript $parameters )
dataDir=$( jq --raw-output .SCC.dataDir $parameters )
workDir=$( jq --raw-output .SCC.workDir $parameters )
imgDir=$( jq --raw-output .dirs.imgDir $parameters )
logDir=$( jq --raw-output .SCC.logDir $parameters ) 
numCores=$( jq --raw-output .SCC.numCores $parameters ) 



tempParams=${parameters}".tmp"

#Add directories to json file that are specific to the tile 
#Will vary depending on if we are doing S10 or HLS30
###############
if $(jq .SCC.runS10 $parameters )
then
  fmask10m="${dataDir}S10/${tile}/fmask/"
	mkdir -p ${fmask10m}	
  cp ${parameters} ${tempParams}; jq --arg fmask10m "$fmask10m" '.dirs.fmask10m = $fmask10m' ${tempParams}>${parameters}
fi


fmaskDir="${dataDir}HLS30/${tile}/fmask/"
tempDir="${workDir}${tile}/temp/"	
chunkDir="${workDir}${tile}/imageChunks/"
phenDir="${workDir}${tile}/phenoMetrics/"

mkdir -p $fmaskDir
mkdir -p $tempDir
mkdir -p $chunkDir
mkdir -p $phenDir

cp ${parameters} ${tempParams}; jq --arg fmaskDir "$fmaskDir" '.dirs.fmaskDir = $fmaskDir' ${tempParams}>${parameters}
cp ${parameters} ${tempParams}; jq --arg tempDir "$tempDir" '.dirs.tempDir = $tempDir' ${tempParams}>${parameters}
cp ${parameters} ${tempParams}; jq --arg chunkDir "$chunkDir" '.dirs.chunkDir = $chunkDir' ${tempParams}>${parameters}
cp ${parameters} ${tempParams}; jq --arg phenDir "$phenDir" '.dirs.phenDir = $phenDir' ${tempParams}>${parameters}

cp ${parameters} ${tempParams}; jq --arg AWS_or_SCC "SCC" '.setup.AWS_or_SCC = "SCC"' ${tempParams}>${parameters}

rm -f $tempParams


#If we are going to preprocess imagery, then empty the chunk directory in case there was a previous run
#It is critical to do this, to prevent any chunks from a previous run from remaining
#Deletes all directories starting with letter "c" in the chunkDir
#This is the safer option to avoid accidently deleting other things if chunkDir is not correctly specified
if $( jq .setup.preprocessImagery $parameters )
then
  rm -rf ${chunkDir}c*
fi


#If S10, need to resample water and topographic files
#Otherwise just copy over
######################
if $( jq .SCC.runS10 $parameters )
then
  gdalwarp -overwrite -r average -ts 10980 10980 -of GTiff /projectnb/modislc/data/water/hls_tiles/water_${tile}.tif ${imgDir}water_${tile}.tif
  gdalwarp -overwrite -r average -ts 10980 10980 -of GTiff /projectnb/modislc/data/dem/usgs_ned/hls_tiles/dem/dem_${tile}.tif ${imgDir}dem_${tile}.tif
  gdalwarp -overwrite -r average -ts 10980 10980 -of GTiff /projectnb/modislc/data/dem/usgs_ned/hls_tiles/slope/slope_${tile}.tif ${imgDir}slope_${tile}.tif
  gdalwarp -overwrite -r average -ts 10980 10980 -of GTiff /projectnb/modislc/data/dem/usgs_ned/hls_tiles/aspect/aspect_${tile}.tif ${imgDir}aspect_${tile}.tif
else
  cp /projectnb/modislc/data/water/hls_tiles/water_${tile}.tif $imgDir
  cp /projectnb/modislc/data/dem/usgs_ned/hls_tiles/dem/dem_${tile}.tif $imgDir
  cp /projectnb/modislc/data/dem/usgs_ned/hls_tiles/slope/slope_${tile}.tif $imgDir
  cp /projectnb/modislc/data/dem/usgs_ned/hls_tiles/aspect/aspect_${tile}.tif $imgDir
fi



runLog="${logDir}${tile}_instanceInfo_${timeStamp}.txt"
errorLog="${logDir}${tile}_errorLog_${timeStamp}.txt"

echo "tile:${tile}"  >"${runLog}"
echo "time:${timeStamp}"  >>"${runLog}"
echo "run:SCC"  >>"${runLog}"
echo "num-cores:${numCores}"  >>"${runLog}"

#Run the script
######################
Rscript $rScript $tile $parameters $runLog $errorLog &
wait


#rm -rf ${inDir}
#rm -rf ${outDir}temp

