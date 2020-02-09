#!/bin/bash -l
#$ -j y


tile=$1
parameters=$2
timeStamp=$3
masterID=$4
logFile=$5


#For now, need to load environmental variables that allows aws commands to work (look into this later)
source /shared/code/MSLSP/AWS/envVars/env_variables_aws.sh


rScript=$( jq --raw-output .AWS.rScript $parameters )
dataDir=$( jq --raw-output .AWS.dataDir $parameters )
workDir=$( jq --raw-output .AWS.workDir $parameters )
s3Dir=$( jq --raw-output .AWS.s3Dir $parameters )
logDir=$( jq --raw-output .AWS.logDir $parameters ) 
numCores=$( jq --raw-output .AWS.numCores $parameters ) 

imgStartYr=$(( $( jq .setup.imgStartYr $parameters ) - 1 ))  #Add one year buffer
imgEndYr=$(( $( jq .setup.imgEndYr $parameters ) + 1 ))      #Add one year buffer

mkdir -p $logDir

tempParams=${parameters}".tmp"

imgDir="${dataDir}${tile}/images/"
fmaskDir="${dataDir}${tile}/fmask/"
tempDir="${workDir}${tile}/temp/"	
chunkDir="${workDir}${tile}/imageChunks/"
phenDir="${workDir}${tile}/phenoMetrics/"

mkdir -p $imgDir
mkdir -p $fmaskDir
mkdir -p $tempDir
mkdir -p $chunkDir
mkdir -p $phenDir

cp ${parameters} ${tempParams}; jq --arg imgDir "$imgDir" '.dirs.imgDir = $imgDir' ${tempParams}>${parameters}
cp ${parameters} ${tempParams}; jq --arg fmaskDir "$fmaskDir" '.dirs.fmaskDir = $fmaskDir' ${tempParams}>${parameters}
cp ${parameters} ${tempParams}; jq --arg tempDir "$tempDir" '.dirs.tempDir = $tempDir' ${tempParams}>${parameters}
cp ${parameters} ${tempParams}; jq --arg chunkDir "$chunkDir" '.dirs.chunkDir = $chunkDir' ${tempParams}>${parameters}
cp ${parameters} ${tempParams}; jq --arg phenDir "$phenDir" '.dirs.phenDir = $phenDir' ${tempParams}>${parameters}

cp ${parameters} ${tempParams}; jq --arg AWS_or_SCC "AWS" '.setup.AWS_or_SCC = "AWS"' ${tempParams}>${parameters}
rm -f $tempParams





#Set up time log
########
id="$(wget -qO- http://instance-data/latest/meta-data/instance-id)"

#Log info on instance
runLog="${logDir}${tile}_instanceInfo_${timeStamp}.txt"
echo "tile:${tile}"  >"${runLog}"
echo "time:${timeStamp}"  >>"${runLog}"
echo "run:AWS"  >>"${runLog}"
echo "num-cores:${numCores}"  >>"${runLog}"
echo "instance-id:${id}"  >>"${runLog}"
echo "instance-type:$(wget -qO- http://instance-data/latest/meta-data/instance-type)"  >>"${runLog}"
echo "reservation-id:$(wget -qO- http://instance-data/latest/meta-data/reservation-id)"  >>"${runLog}"
echo "ami-id:$(wget -qO- http://instance-data/latest/meta-data/ami-id)"  >>"${runLog}"
echo "master-id:${masterID}"  >>"${runLog}"
aws s3 cp ${runLog} ${s3Dir}logs/instanceInfo/ --quiet 

errorLog="${logDir}${tile}_errorLog_${timeStamp}_${id}.txt"
touch ${errorLog}



#Get Images from S3
##########
for y in `seq $imgStartYr $imgEndYr`;
do
l30_path="s3://hlsp1/HLS_LOCAL_IO.v1.4/L30/${y}/${tile:0:2}/${tile:2:1}/${tile:3:1}/${tile:4:1}/" 
s30_path="s3://hlsp1/HLS_LOCAL_IO.v1.4/S30/${y}/${tile:0:2}/${tile:2:1}/${tile:3:1}/${tile:4:1}/"  
aws s3 cp ${l30_path} ${imgDir} --recursive --quiet --request-payer requester
aws s3 cp ${s30_path} ${imgDir} --recursive --quiet --request-payer requester
done



#Get Ancillary from S3
##########
aws s3 cp s3://mslsp/inputs/water/water_${tile}.tif ${imgDir} --quiet 
aws s3 cp s3://mslsp/inputs/dem/dem_${tile}.tif ${imgDir} --quiet 
aws s3 cp s3://mslsp/inputs/slope/slope_${tile}.tif ${imgDir} --quiet 
aws s3 cp s3://mslsp/inputs/aspect/aspect_${tile}.tif ${imgDir} --quiet 


#Run the script
######################
#Load environmental variables that allows matlab to run
source /shared/code/MSLSP/AWS/envVars/env_variables.sh
Rscript $rScript $tile $parameters $runLog $errorLog &
wait

#After, we need to copy results to S3, including logs and assessment layers
#Load environmental variables that allows aws to run again
source /shared/code/MSLSP/AWS/envVars/env_variables_aws.sh

copyTo="${s3Dir}product/${tile}/"
aws s3 cp ${phenDir} ${copyTo} --recursive --quiet            #Copy all result layers
aws s3 cp ${runLog} ${copyTo} --quiet                        #Copy run log to output folder
aws s3 cp ${parameters} "${copyTo}parameters_${tile}_${timeStamp}.json" --quiet   #Copy parameters that were used


if $(jq .AWS.keepChunks $parameters )
then
    aws s3 cp ${chunkDir} ${s3Dir}imageChunks/${tile}/ --recursive --quiet 
fi



if [ -s ${errorLog} ]
then
    aws s3 cp ${errorLog} ${s3Dir}logs/errorLogs/ --quiet 
fi

aws s3 cp ${runLog} ${s3Dir}logs/instanceInfo/ --quiet            #Copy logs
aws s3 cp ${logFile} ${s3Dir}logs/completeLogs/ --quiet 








