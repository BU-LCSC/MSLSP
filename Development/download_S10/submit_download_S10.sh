

tileList="tileLists/S10_starting_set.txt"
imgStartYr=2015
imgEndYr=2020

baseDir='/projectnb/modislc/projects/landsat_sentinel/v1_4/S10/'

jobTime="$(date +%Y_%m_%d_%H_%M_%S)"

while read -r tile
do
    outDir="${baseDir}${tile}/images/"
    mkdir -p ${outDir}
    nameArg="-N S3_${tile}"
    downloadArg="-l download"
    logArg="-o "logs/"Download_${tile}_${jobTime}.txt"
    qsub $nameArg $downloadArg $logArg download_S10.sh $tile $outDir $imgStartYr $imgEndYr
done < $tileList
