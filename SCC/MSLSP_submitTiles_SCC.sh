

parameters="/usr3/graduate/mkmoon/GitHub/MSLSP/MSLSP_Parameters.json"
tileList="tileLists/europe.txt"


numCores=$( jq .SCC.numCores $parameters )
workDir=$( jq --raw-output .SCC.workDir $parameters )
dataDir=$( jq --raw-output .SCC.dataDir $parameters )

jobTime="$(date +%Y_%m_%d_%H_%M_%S)"

logDir=$( jq --raw-output .SCC.logDir $parameters ) 
mkdir -p $logDir

if $(jq .SCC.runS10 $parameters )
then
    baseDir="${dataDir}S10/"
    nodeArgs="-l h_rt=24:00:00 -l mem_per_core=16G -pe omp ${numCores}"
else
    baseDir="${dataDir}HLS30/"
    nodeArgs="-l h_rt=12:00:00 -pe omp ${numCores}"
fi


#Set up initial folders and copy parameter file 
while read -r tile
    do
        tileDir="${workDir}${tile}/"
        imgDir="${baseDir}${tile}/images/"
        mkdir -p $tileDir
        mkdir -p $imgDir
        
        paramName="${tileDir}parameters_${jobTime}.json"
        jq --arg imgDir "$imgDir" '.dirs.imgDir = $imgDir' ${parameters}>${paramName}
        
    done < $tileList


#If download is set to true, the code will check for new HLS images and only download those
if [ $(jq .setup.downloadImagery $parameters ) == true ] && [ $(jq .SCC.runS10 $parameters ) == false ]
then
    while read -r tile
    do
        nameArg="-N DL_${tile}"
        logArg_download="-o ${logDir}Download_${tile}_${jobTime}.txt"
        downloadArg="-l download"
        imgStartYr=$(( $( jq .setup.imgStartYr $parameters ) - 1 ))  #Add one year buffer
        imgEndYr=$(( $( jq .setup.imgEndYr $parameters ) + 1 ))      #Add one year buffer
        qsub $nameArg $logArg_download $downloadArg runDownloadHLS.sh $tile $baseDir $imgStartYr $imgEndYr
    done < $tileList

    while read -r tile
    do
        tileDir="${workDir}${tile}/"
        paramName="${tileDir}parameters_${jobTime}.json"
        nameArg="-N R_${tile}"
        logArg="-o ${logDir}Run_${tile}_${jobTime}.txt"
        holdArg="-hold_jid DL_${tile}"
        qsub $nameArg $logArg $nodeArgs $holdArg MSLSP_runTile_SCC.sh $tile $paramName $jobTime
    done < $tileList
else
    while read -r tile
    do
        tileDir="${workDir}${tile}/"
        paramName="${tileDir}parameters_${jobTime}.json"
        nameArg="-N R_${tile}"
        logArg="-o ${workDir}Run_${tile}_${jobTime}.txt"
        qsub $nameArg $logArg $nodeArgs MSLSP_runTile_SCC.sh $tile $paramName $jobTime
    done < $tileList
fi 

