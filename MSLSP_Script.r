#Douglas Bolton, Boston University
#Main script for running HLS Land Surface Phenology
#######################################################################################


start_time <- Sys.time()



#Load required libraries
##################################
print('----------------------------------------------------------------------------------------------')
print('Loading packages')
print('')
print('')

library(raster)
library(rgdal)
library(gdalUtils)
library(rgeos)
library(imager)   #needed for efficient distance to snow calculate

library(iterators)
library(foreach)
library(doMC)

library(matrixStats)   
library(WGCNA)
library(zoo)



print('----------------------------------------------------------------------------------------------')
print('Start processing')
print('')
print('')


#Read in arguments
args <- commandArgs(trailingOnly=T)
print(args)
tile <- args[1]
inDir <- args[2]
outDir <- args[3]
outLog <- args[4]
errorLog <- args[5]
functionPath <- args[6]
fmaskFunction <- args[7]
phenStartYr <- as.numeric(args[8]) 
phenEndYr <- as.numeric(args[9]) 
numCores <- as.numeric(args[10])
numChunks <- as.numeric(args[11])
topoMethod <- args[12]
dataFill <- as.logical(args[13])
stepsToDo <- as.numeric(args[14])

#Load functions
source(file=functionPath)

#Register the parallel backend
registerDoMC(cores=numCores)

#Get default parameters
pheno_pars <- DefaultPhenoParameters_HLS()



time_step1 <- 0
time_step2 <- 0

#Set up data
#####################################################

##Make required directories
########################3
tempFold <- paste0(outDir,"temp/")
phenoFold <- outDir      #Just to avoid a long folder structure, we will put phenology year folders right in out directory

#Make fmask folder, chunk folders, temporary output folders, and final output folders
####################
if (!dir.exists(paste0(tempFold,"fmask"))) {dir.create(paste0(tempFold,"fmask"),recursive=T)}

for (i in 1:numChunks) {dirName <- paste0(tempFold,"c",i,'/')
  if (!dir.exists(dirName)) {dir.create(dirName,recursive=T)}}

for (y in seq(phenStartYr,phenEndYr)) {for (i in seq(1,pheno_pars$numLyrs)) {dirName <- paste0(tempFold,'outputs/y',y,'/lyr',i,'/')
  if (!dir.exists(dirName)) {dir.create(dirName,recursive=T)}}}

for (y in seq(phenStartYr,phenEndYr)) {dirName <- paste0(phenoFold,'y',y)
  if (!dir.exists(dirName)) {dir.create(dirName,recursive=T)}}
  
  
#Get full image list
###########
imgList <- list.files(path=inDir, pattern=glob2rx("HLS*.hdf"), full.names=T)

#Get the year and doy of each image
##########################
yrdoy = as.numeric(matrix(NA,length(imgList),1))
for (i in 1:length(imgList)) {
  imgName_strip = tail(unlist(strsplit(imgList[i],'/')),n = 1)
  yrdoy[i]= unlist(strsplit(imgName_strip,'.',fixed = T))[4]}
doys <- as.numeric(format(as.Date(strptime(yrdoy, format="%Y%j")),'%j'))
years <- as.numeric(format(as.Date(strptime(yrdoy, format="%Y%j")),'%Y'))
imgYrs <- sort(unique(years))  #What years do we actually have imagery for?

#Get raster information from first image
##########################
qaName  <-  paste0('HDF4_EOS:EOS_GRID:',imgList[1], ':Grid:QA')
ref_info <- gdalinfo(qaName,proj4=TRUE,raw_output=FALSE)
ts <- c(ref_info$columns,ref_info$rows) #Get image rows and columns
numPix  <-  ts[1] * ts[2]   #Get total number of pixels
baseImage  <-  raster(qaName) #Set up base image that we'll use for outputs

#Sort out chunk boundaries
#################
boundaries <- chunkBoundaries(numPix,numChunks) #Chunk image for processing
chunkStart <- boundaries[[1]]   #Get pixel boundaries for each chunk
chunkEnd <- boundaries[[2]]

#Read in water mask
water <- readGDAL(paste0(inDir,'water_',tile,'.tif'),silent=T)$band1
waterMask <- water == 2 | water == 0    #Mask water and zero (zero = ocean far from shore)

#STEP 1 - Preprocess images
#####################################################
#####################################################
if (stepsToDo == 1 | stepsToDo == 3) {    #1 = Do preprocessing steps only. 3 = Do both image preprocessing and phenology runs
  
  #Apply mask and write out chunked images. Currently using the QA_and_Fmask version  (QA for Landsat, Fmask 4.0 for Sentinel)
  imgLog <- foreach(j=1:length(imgList),.combine=c) %dopar% {
                log <- try({ApplyMask_QA_and_Fmask(imgList[j], tile, waterMask, tempFold, inDir, fmaskFunction, chunkStart, chunkEnd, topoMethod, pheno_pars, deleteInput=F)},silent=T)
                if (inherits(log, 'try-error')) {cat(paste('ApplyMask_QA_and_Fmask: Error for', imgList[j],'\n'), file=errorLog, append=T)}  #If there's an error, keep going, but write to error log
                }   
  
  #Topographic Correction of images
  #TopoMethod is either set to VI or None. If set to None, this step is skipped.
  ###########################
  if (topoMethod == 'VI') {
    #In this approach, pixels will be clustered (kmeans) according to specified vegetation indices
    #Once clusters are determined, each cluster/band combo will be topographically corrected separately.
  
    #Read in slope and aspect rasters
    #IMPORTANT: Code expects units of slope and aspect to be radians * 10000
    slope <- raster(paste0(inDir,'slope_',tile,'.tif')) #Keeping this slope raster as a template for other temporary outputs
    slopeVals <- readGDAL(paste0(inDir,'slope_',tile,'.tif'),silent=T)$band1
    slopeVals[slopeVals == 65534] = NA
    slopeVals = slopeVals / 10000
    
    aspectVals = readGDAL(paste0(inDir,'aspect_',tile,'.tif'),silent=T)$band1
    aspectVals[aspectVals == 65534] = NA
    aspectVals = aspectVals / 10000
    
    for (yr in imgYrs) {
      subList <- imgList[years==yr]
      subDOY  <- doys[years==yr]
      
      if (length(subList) > 0) {
        #Check if the first year has enough data (images back to DOY 60 or whatever specified). If it doesn't, then use the next year's VIs for kmeans
        if (yr == imgYrs[1] & min(subDOY) > pheno_pars$requiredDoyStart) {yrPull <- yr+1} else {yrPull <- yr}
        #Check if the final year has enough data (images past DOY 300 or whatever specified). If it doesn't, then use the previous year's VIs for kmeans
        if (yr == imgYrs[length(imgYrs)] & max(subDOY) < pheno_pars$requiredDoyEnd) {yrPull <- yr-1} else {yrPull <- yr}
  
        #Calculate the VIs by reading in image chunks, and calcuating percentiles
        indexImg <- foreach(j=1:numChunks,.combine=rbind) %dopar% {
          numPix <- chunkEnd[j] - chunkStart[j] + 1  #need to know number of pixels in a chunk, in case there are no pixels and we need to fill
          getIndexQuantile(j, numPix, tempFold, yrPull, errorLog, pheno_pars)}
      
        #Convert VIs to z-scores prior to kmeans
        for (i in 1:dim(indexImg)[2]) {
              col <- indexImg[,i]
              indexImg[,i] <- (col - mean(col,na.rm=T)) / sd(col,na.rm=T)
        }
        
        #Perform kmeans. Must first remove NA values (otherwise kmeans fails). Topo correction will be performed for each class
        goodPix <- !is.na(rowMeans(indexImg))
        kClust <- kmeans(indexImg[goodPix,], pheno_pars$kmeansClasses, iter.max = pheno_pars$kmeansIterations)
        groups <- matrix(0,numPix)
        groups[goodPix] <-  kClust$cluster
        groups[is.na(groups) | is.na(slopeVals) | is.na(aspectVals)]  <- 0    #Assign zero value to groups if slope, aspect, or group is NA
  
        #Now that we have kmeans classes, run the topographic correction function. 
        imgLog <- foreach(j=1:length(subList),.combine=c) %dopar% {
                 log <- try({runTopoCorrection(subList[j], tempFold, groups, slopeVals, aspectVals, chunkStart,chunkEnd, errorLog, pheno_pars, writeImages=FALSE, slope=NA)}, silent=T)
                 if (inherits(log, 'try-error')) {cat(paste('RunTopoCorrection: Error for', subList[j],'\n'), file=errorLog, append=T)} 
        }
      }
    }
  }
  
  time_step1 <- as.numeric(difftime(Sys.time(),start_time,units="mins"))
  print('----------------------------------------------------------------------------------------------')
  print(paste('Finished preprocessing images in',round(time_step1,1),'minutes'))
  print('')
  print('')
  
}




#STEP 2 - Run Phenology
#####################################################
#####################################################
if (stepsToDo == 2 | stepsToDo == 3) {   #2 = Do phenology step only. 3 = Do both image preprocessing and phenology runs

  #Run phenology code for each image chunk
  imgLog <- foreach(j=1:numChunks) %dopar% {
    try({runPhenoChunk(j, tempFold, phenStartYr, phenEndYr, errorLog, pheno_pars, dataFill)},silent=T)
    if (inherits(log, 'try-error')) {cat(paste('RunPhenoChunk: Error for chunk', j,'\n'), file=errorLog, append=T)}
  }
  
  
  #Output Rasters
  ##################

  
  #Set up raster names
  #########
  names <- c('OGI','50PCGI','OGMx','OGD','50PCGD','OGMn',
             paste0('OGI_b',1:6),paste0('50PCGI_b',1:6),paste0('OGMx_b',1:6),
             paste0('OGD_b',1:6),paste0('50PCGD_b',1:6),paste0('OGMn_b',1:6),
             'EVImin','EVImax','EVIarea',
             'gupRsq','gupNumObs','gupMaxGap','gupMaxGapFilled',
             'gdownRsq','gdownNumObs','gdownMaxGap','gdownMaxGapFilled')
  
  lyrNames <- c('NumCycles',names,paste0(names,'_2')) #Duplicate columns for the second cycle
  
  #What columns need to be converted to day of year?
  doyColumns <- c('OGI','50PCGI','OGMx','OGD','50PCGD','OGMn','OGI_2','50PCGI_2','OGMx_2','OGD_2','50PCGD_2','OGMn_2')
  
  
  #Loop through the years processed
  yrs <- phenStartYr:phenEndYr  
  for (yr in yrs) {
    startDate <- as.numeric(as.Date(paste0(yr-1,'-12-31')))   #Determine date to subtract from each raster to get DOY
    foreach(i = 1:pheno_pars$numLyrs) %dopar% {     #Loop through each layer name
      lyrName <- lyrNames[i]
      mat <- matrix(NA,numPix,1)
      
      #Reconstruct image from chunks
      for (n in 1:numChunks) {
        fileName <- paste0(tempFold,'outputs/y',yr,'/lyr',i,'/c',n,'.Rds')
        matSub <- try(readRDS(fileName),silent=T)
        if (inherits(matSub, 'try-error')) {next} 
        mat[chunkStart[n]:chunkEnd[n]] <- matSub
      }
      
      if (lyrNames[i] %in% doyColumns) {mat <- mat - startDate}   #If layer needs to be DOY, then convert to DOY
      
      mat <- as.integer(round(mat))    #Convert to integer
      rast <- setValues(baseImage,mat)  #Convert to raster using base image
      
      #Determine data type and fill value (only NumCycles layer will be 8-bit)
      if (lyrNames[i] == 'NumCycles') {dType <- 'INT1U';naVal <- 255    
      } else {dType <- 'INT2S';naVal <- 32767}
      
      names(rast) <- paste0('LSP_',tile,'_',yr,'_',lyrNames[i])    #Assign a layer name to the raster
      outName <- paste0(phenoFold,'y',yr,'/LSP_',tile,'_',yr,'_',lyrNames[i],'.tif')     #Create an output name
      writeRaster(rast,filename=outName,format='GTiff',overwrite=TRUE,datatype=dType,NAflag=naVal)  #Write the raster to disk
    }
  }
  
  
  #Create 8-bit QA layers (Work in progress!)
  #foreach(yr = yrs) %dopar% {CreateQA(tile,yr,phenoFold,waterMask,baseImage,pheno_pars)}
    
  
  total_time <- as.numeric(difftime(Sys.time(),start_time,units="mins"))
  time_step2 <- total_time - time_step1
  print('----------------------------------------------------------------------------------------------')
  print(paste('Finished detecting phenology in',round(time_step2,1),'minutes'))
  print('')
  print('')
  
}


#Write to log
######################
total_time <- as.numeric(difftime(Sys.time(),start_time,units="mins"))
line <- paste(tile, time_step1, time_step2, total_time, sep=',')
line <- paste0(line,'\n')
cat(line,file=outLog,append=T)






