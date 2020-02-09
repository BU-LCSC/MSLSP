#Load required libraries
library(raster)
library(rgdal)
library(gdalUtils)
library(rgeos)
library(foreach)
library(iterators)
library(doMC)

library(matrixStats)
library(WGCNA)
library(zoo)

library(RcppRoll)

library(rjson)

library(ncdf4)



#Inputs
############

#Tile 
tileList <-'/usr2/postdoc/dbolt/CodeFolder/git_repos/MSLSP/SCC/tileLists/anom_tiles.txt'
tileList <-'/usr2/postdoc/dbolt/CodeFolder/git_repos/MSLSP/SCC/tileLists/kansas_tiles.txt'

numCores <- 28

imgYrs <- 2016:2019
phenYrs <- 2017:2019

#What data folder
dataBase <- '/projectnb/modislc/projects/landsat_sentinel/MSLSP_HLS30/'
imgBase <-   '/projectnb/modislc/projects/landsat_sentinel/v1_4/HLS30/'


#Define path to function
functions <- '/usr2/postdoc/dbolt/CodeFolder/git_repos/MSLSP/MSLSP_Functions.r'
functions_diagnostics <- '/usr2/postdoc/dbolt/CodeFolder/git_repos/MSLSP/supporting/MSLSP_Diagnostic_Functions_V1_0.r'


#Json file where phenology parameters are defined. Will just use the phenology paramaters from this 
jsonFile <- "/projectnb/modislc/users/dbolt/Crop_Analysis/MSLSP_Parameters_CDL_Runs.json"

cdlFold <- '/projectnb/modislc/data/lc_database/regional/united_states/cropland_data_layer/hls_tiles/'

cropTable <- read.csv('/projectnb/modislc/users/dbolt/Crop_Analysis/cdl_classes.csv',header=F,stringsAsFactors = F)

outFold <- '/projectnb/modislc/users/dbolt/Crop_Analysis/anomTiles/'
outFold <- '/projectnb/modislc/users/dbolt/Crop_Analysis/kansasTiles/'

numSamp <- 1000

#Load functions 
source(file=functions)
source(file=functions_diagnostics)


#Get default parameters
params <- fromJSON(file=jsonFile)


registerDoMC(cores=numCores)

#
tiles <- unlist(read.table(tileList,header=F,stringsAsFactors=F),use.names=F)


#Alter the jsonFile?
params$phenology_parameters$splineSpar=0.4      #Alter the smoothing parameter


crops <- c('Corn', 'Soybeans', 'Winter Wheat') 
crops <- c('Corn', 'Soybeans', 'Winter Wheat', 'Rye', 'Sorghum', 'Cotton', 'Canola', 'Millet', 'Alfalfa') 

print(tileList)
print(outFold)
print(crops)

for (crop in crops) {
  cropID <- cropTable[cropTable[,2]==crop,1]

for (y in 1:length(phenYrs)) {
  yr <- phenYrs[y]

  pixList <- c()

  for (t in 1:length(tiles)) {
    tile <- tiles[t]
    
    #First read in QA and only take mod/high quality pixels
    #dataFile <- paste0(dataBase,tile,'/phenoMetrics/y',yr,'/MS-LSP_',tile,'_',yr,'.nc')
    #ncFile <- nc_open(dataFile)
    #qa <- ncvar_get(ncFile, 'overallQA')
    #qa <- matrix(qa,nrow(qa)*ncol(qa))
    
    cdl <- readGDAL(paste0(cdlFold,yr,'/cdl_',tile,'.tif'),silent=T)$band1
    pixIDs <- 1:length(cdl)
    check <- pixIDs[cdl == cropID]
    
    #check <- pixIDs[cdl == cropID & qa < 4]
    
    tab <- cbind(rep(t,length(check)),check)
    
    pixList <- rbind(pixList,tab)
  }
  
  if (numSamp > length(pixList[,1])) {numSamp2 <- length(pixList[,1])} else {numSamp2 <- numSamp}
  samp <- sample(1:length(pixList[,1]),numSamp2)

  chosen <- pixList[samp,]

  uniq <- unique(chosen[,1])

  numDaysFit <- 365 + params$phenology_parameters$splineBuffer*2
  
  smoothTable <- matrix(NA,numDaysFit,nrow(chosen)+1)
  phenTable <- matrix(NA,nrow(chosen),9)
  count=0
  for (t in uniq) {
    tile <- tiles[t]
    #Where is the data?
    imgDir <- paste0(imgBase,tile,'/images/')
    chunkDir <- paste0(dataBase,tile,'/imageChunks/')
    
    pixelNumbers <- chosen[chosen[,1]==t,2]
    theTable <- Extract_Timeseries(tile, imgDir, chunkDir, imgYrs, phenYrs, numCores, params, pixelNumbers=pixelNumbers)
    
    #Just extract smooth curve for the correct year. Append to running lst
    for (i in 1:length(theTable)) {
      log <- try({
      subTab <- theTable[[i]]
      subTab <- subTab[[paste0('y',yr)]]
      smooth <- subTab$smoothed_vi
      if (length(subTab$phenDates) >= 7) {    #if length = 7, that means phenology was detected
        count <- count+1
        phenTable[count,] <- c(subTab$phenDates[1:7],subTab$stats$gup_maxgap_frac,subTab$stats$gdown_maxgap_frac)
        smoothTable[,count+1] <- subTab$smoothed_vi
        if (count ==1) {smoothTable[,1] <- subTab$smoothed_dates}
      }
      },silent=TRUE)
    }
  }
  
  phenTable <- data.frame(phenTable)
  colnames(phenTable) <- c('OGI','PCGI','OGMx','Peak','OGD','PCGD','OGMn','gupGap','gdownGap')
  outName <- paste0(outFold,'Phen_Dates_',crop,'_',yr,'.csv')
  write.table(phenTable,outName,sep=",",row.names=F)
  
  smoothTable <- data.frame(smoothTable)
  outName <- paste0(outFold,'Smoothed_timeseries_',crop,'_',yr,'.csv')
  write.table(smoothTable,outName,sep=",",row.names=F,col.names=F)
}
}
