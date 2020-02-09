
tileList=$1

parameters="/shared/code/MSLSP/MSLSP_Parameters.json"

numCores=$( jq .AWS.numCores $parameters )
logDir=$( jq --raw-output .AWS.logDir $parameters )

masterID="$(wget -qO- http://instance-data/latest/meta-data/instance-id)"

jobTime="$(date +%Y_%m_%d_%H_%M_%S)"


while read -r tile
do 
    paramName="${logDir}parameters_${tile}_${jobTime}.json"
    cp $parameters $paramName
    nameArg="-N T${tile}"
    logFile="${logDir}${tile}_completeLog_${jobTime}.txt"
    logArg="-o ${logFile}"
    coreArg="-pe smp ${numCores}"
    qsub ${nameArg} ${logArg} ${coreArg} MSLSP_runTile_AWS.sh $tile $paramName $jobTime $masterID $logFile
done < $tileList 
