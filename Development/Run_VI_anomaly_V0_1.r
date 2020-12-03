#######################################################################################

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
library(ncdf4)

library(iterators)
library(foreach)
library(doMC)

library(matrixStats)   
library(WGCNA)
library(zoo)

library(RcppRoll)

library(rjson)

#---------------------------------------------------------------------
#Calculate VI anomaly for each pixel
#Code adapted by Minkyu Moon
#Based on a function "DoPhenologyHLS"
#---------------------------------------------------------------------
DoAnomalyHLS <- function(b2, b4, vi, snowPix, dates, imgYrs, phenYrs, splineStart, splineEnd, numDaysFit, pheno_pars){
  
  
  #Despike, calculate dormant value, fill snow with dormant value, despike poorly fit snow values
  log <- try({
    #Despike
    spikes <- CheckSpike_MultiBand(b2/10000,b4/10000,vi, dates, pheno_pars)
    
    #Set despiked values to NA
    vi[spikes] <- NA
    
    dormIms <- dates >= pheno_pars$dormStart & dates <= pheno_pars$dormEnd
    
    vi_dorm <- quantile(vi[dormIms],probs=pheno_pars$dormantQuantile,na.rm=T)   #Calc vi dormant value
    snowPix <- Screen_SnowFills(vi,vi_dorm,snowPix,dates,pheno_pars)              #Screen poorly filled snow values
    
    vi[snowPix] <- vi_dorm   #Fill remaining snow values with dormant value for vi  
    
    #Determine gaps that require filling
    gDates <- dates[!is.na(vi)]  
    dDiff <- diff(gDates) > pheno_pars$gapLengthToFill     #Gaps greater than 30 days will be filled
    dStart <- gDates[c(dDiff,FALSE)]
    dEnd <- gDates[c(FALSE,dDiff)]
    
    
    #Locate gaps in date vector
    all_dates <- seq(min(splineStart), max(splineEnd), by="day")
    
    fill_locations <- matrix(FALSE,length(all_dates))
    for (d in 1:length(dStart)) {
      fill_locations[all_dates >= dStart[d] & all_dates < dEnd[d]] <- TRUE}
    
    fill_dates <- all_dates[fill_locations]
    
    yToDo <- 1:length(imgYrs)
    yToDo <- yToDo[imgYrs %in% phenYrs]
    yrsWithGaps <- c()
    for (y in yToDo) {
      pred_dates <- seq(splineStart[y], splineEnd[y], by="day")
      if (sum(pred_dates %in% fill_dates) > 0) {yrsWithGaps <- c(yrsWithGaps,imgYrs[y])}
    }
    
    
    #If there are gaps to be filled, then we will spline all years. 
    #If not, just spline product years
    if (length(yrsWithGaps) > 0) {
      yrs <- imgYrs
    } else {
      yrs <- phenYrs
      splineStart <- splineStart[imgYrs %in% phenYrs]
      splineEnd <- splineEnd[imgYrs %in% phenYrs]}
    
    numYrs <- length(yrs)
    daysVec <- 1:numDaysFit
    vecLength <- numDaysFit*numYrs
    
    
    
    
    
    #First, we will fit splines to each year invidually
    #To line up observations from each year, we will create a matrix for vi and each band (numDaysFit x numYears)
    
    smoothMat <- matrix(NA, numDaysFit, numYrs)
    maskMat <- matrix(0, numDaysFit, numYrs)
    fillMat <- smoothMat
    baseWeights <- maskMat
    
    for (y in 1:numYrs) {
      #Use try statement, because we don't want to stop processing if only an error in one year
      try({
        
        dateRange <- dates >= splineStart[y] & dates <= splineEnd[y] & !is.na(vi)   
        
        dateSub <- dates[dateRange]; viSub <- vi[dateRange]; snowSub <- snowPix[dateRange]
        
        #Get weights
        weights <- matrix(1,length(snowSub))
        weights[snowSub == 1] <- pheno_pars$snowWeight
        
        pred_dates <- seq(splineStart[y], splineEnd[y], by="day")
        
        #Assign weights and run cubic spline
        smoothed <- Smooth_VI(viSub, dateSub, pred_dates, weights, pheno_pars, vi_dorm)
        
        
        #Mask spline in gaps, and before/after first/last image
        maskMat[fill_locations[all_dates %in% pred_dates]] <- 1    #Mask spline in gaps
        maskMat[pred_dates < dateSub[1]] <- 1                      #Mask spline before first image and after last image
        maskMat[pred_dates > dateSub[length(dateSub)]] <- 1
        
        #Mask spline in the buffer years (only interested in comparing splines in target year)
        maskMat[format(pred_dates,'%Y') != yrs[y]]  <- 1
        
        fillDs <- pred_dates %in% dateSub
        
        smoothMat[,y] <- smoothed
        baseWeights[fillDs,y] <- weights
        fillMat[fillDs,y] <- viSub
        
      },silent=TRUE)
    }
    
    
    xs <- rep(daysVec,numYrs)
    ys <- matrix(fillMat,vecLength)
    ysGood <- !is.na(ys)
    baseW <- matrix(baseWeights,vecLength)   #Base Weights are 1=clear observation, 0.5=snow-filled
    
    smoothMat_Masked <- smoothMat
    smoothMat_Masked[maskMat] <- NA
    
    
    #Loop through years, compare spline to other years, weight each year based on similiarity, fit spline, calculate phenology
    #Just product years now
    yToDo <- 1:numYrs
    yToDo <- yToDo[yrs %in% phenYrs]
    
    
    weightArray <- calculateWeights(smoothMat_Masked, numDaysFit, numYrs, pheno_pars) 
    
    prevYear <- daysVec <= pheno_pars$splineBuffer
    inYear <- daysVec > pheno_pars$splineBuffer & daysVec <= (pheno_pars$splineBuffer+365)
    nextYear <- daysVec > (pheno_pars$splineBuffer+365)
    
  },silent=TRUE)
  #If there is an error despiking or other initial steps, return NAs
  if(inherits(log, "try-error")){return(matrix(NA,365*length(imgYrs)))}   
  
  if(length(imgYrs)==length(yToDo)){
    vi_mat <- matrix(NA,365,length(yToDo))
    for (y in yToDo) {
      log <- try({
        
        pred_dates <- seq(splineStart[y], splineEnd[y], by="day")
        
        
        if (yrs[y] %in% yrsWithGaps) {
          
          indPrev <- y-1; indPrev[indPrev<1] <- 1
          indNext <- y+1; indNext[indNext>numYrs] <- numYrs
          
          weights <- rbind(weightArray[prevYear,,indPrev],
                           weightArray[inYear,,y],
                           weightArray[nextYear,,indNext])
          
          #Where are the gaps?
          toFill <- fill_locations[all_dates %in% pred_dates]
          
          weights[!toFill,] <- 0     #Set weight to zero for observations that aren't in a gap
          weights[,y] <- 1           #Set weights in target year to 1
          
          
          #Now that we have weights, calculate phenology
          #######################
          weights <- matrix(weights,vecLength) * baseW   #Multiple weights by base weight (1=good,0.5=snow-filled)
          theInds <- ysGood & weights > 0
          xs_sub <- xs[theInds]; w_sub <- weights[theInds]
          smoothed_vi <- Smooth_VI(ys[theInds], xs_sub, daysVec, w_sub, pheno_pars, vi_dorm)  #Fit spline
          
        } else {
          
          #Variables needed for next steps if the above gap filling was not done
          theInds <- matrix(FALSE,length(ysGood))
          theInds[((y-1)*numDaysFit+1):(y*numDaysFit)] <- TRUE
          xs_sub <- xs[theInds]; w_sub <- baseW[theInds]
          
          smoothed_vi <- smoothMat[,y]   #if no gaps to fill, just use existing spline
        }
        
        vi_mat[,y] <- smoothed_vi[(pheno_pars$splineBuffer+1):(pheno_pars$splineBuffer+365)]
        
      },silent=TRUE)  #End of the try block
    }
    
    vi_mean <- apply(vi_mat,1,mean,na.rm=T)
    out <- vi_mat - vi_mean
    out <- matrix(out,365*length(yToDo))
  }else{
    out <- matrix(NA,365*length(imgYrs))
  }
  
  return(out)
}


#---------------------------------------------------------------------
#Run Phenology code for an image chunk
#Write phenology results for each chunk to disk
#Douglas Bolton
#---------------------------------------------------------------------
runAnomalyChunk <- function(chunk, numPix, imgYrs, phenYrs, errorLog, params) {
  
  pheno_pars <- params$phenology_parameters
  
  #Get all images to process
  ######################
  chunkFold <- paste0(params$dirs$chunkDir,'c',chunk,'/') 
  outFold <-   paste0(params$dirs$tempDir,'outputs/')  
  
  imgList <- list.files(path=chunkFold, pattern=glob2rx("HLS_*.Rds"), full.names=F)
  
  numImgs = length(imgList)
  
  yrdoy = as.numeric(matrix(NA,numImgs,1))
  for (i in 1:numImgs) {
    imName <- gsub('.Rds','',imgList[i])
    yrdoy[i]= as.numeric(unlist(strsplit(imName,'_',fixed = T))[4])}
  
  ord = order(yrdoy)  #Determine image order
  
  #Read in all imagery for chunk
  b2 <- matrix(NA,numPix,numImgs);  b3 <- matrix(NA,numPix,numImgs);  b4 <- matrix(NA,numPix,numImgs)
  b5 <- matrix(NA,numPix,numImgs);  b6 <- matrix(NA,numPix,numImgs);  b7 <- matrix(NA,numPix,numImgs)
  re1 <- matrix(NA,numPix,numImgs); re2 <- matrix(NA,numPix,numImgs); re3 <- matrix(NA,numPix,numImgs) 
  
  sensor <- as.character(rep('',numImgs))
  
  for (i in 1:length(ord)) {
    img <- imgList[ord[i]]
    imgData <- try(matrix(readRDS(paste0(chunkFold,img)),nrow=numPix),silent = TRUE)
    if (inherits(imgData, 'try-error')) {cat(paste('runPhenoChunk: Error for chunk',chunk,img), file=errorLog, append=T);next} 
    
    b2[,i] <- imgData[,1]; b3[,i] <- imgData[,2]; b4[,i] <- imgData[,3]
    b5[,i] <- imgData[,4]; b6[,i] <- imgData[,5]; b7[,i] <- imgData[,6]
    
    sensor[i] = unlist(strsplit(img,'_',fixed = T))[2]
    
    #If there are more than 7 bands, fill in the red edge bands
    if (sensor[i] == 'S30' | sensor[i] == 'S10') {
      re1[,i] <- imgData[,7]; re2[,i] <- imgData[,8]; re3[,i] <- imgData[,9]
    }
    
    remove(imgData)
  }
  
  dates = as.Date(strptime(yrdoy[ord], format="%Y%j"))  #Format as dates
  
  #Get snow pixels to fill 
  snowPix <- b2 == pheno_pars$snowFillVal
  snowPix[is.na(snowPix)] <- FALSE
  
  #Set snow pixels to NA
  b2[snowPix] <- NA; b3[snowPix] <- NA; b4[snowPix] <- NA
  b5[snowPix] <- NA; b6[snowPix] <- NA; b7[snowPix] <- NA 
  re1[snowPix] <- NA; re2[snowPix] <- NA; re3[snowPix] <- NA 
  
  
  #Calculate index, mask negative values
  vi <- calcIndex(blue=b2/10000,green=b3/10000,red=b4/10000,
                  nir=b5/10000,swirOne=b6/10000,swirTwo=b7/10000,
                  edge1=re1/10000,edge2=re2/10000,edge3=re3/10000,
                  whatIndex=pheno_pars$vegetation_index)
  if (pheno_pars$maskNegativeValues) {vi[vi < 0] <- NA}       #Option to mask negative values (defined in json file)
  
  #Average images that occur on the same day
  #Only keep the snow flag if both images on the day were classified as snow
  imgData <- averageDuplicates(list(b2, b3, b4, b5, b6, b7, vi), snowPix, sensor, dates)
  
  #Extract data from listd
  bands <- imgData[[1]]; snowPix <- imgData[[2]]; sensor <- imgData[[3]]; dates <- imgData[[4]] 
  b2 <- bands[[1]];  b4 <- bands[[3]]  
  vi <- bands[[7]]
  remove(imgData)
  remove(bands)
  
  
  #Define the start and end to the buffer period
  #If it is a leap year, there will be one more day in target year, one less day in the following buffer period
  splineStart <- as.Date(as.Date(paste0(imgYrs,'-01-01')) - pheno_pars$splineBuffer) 
  numDaysFit  <-  365 + (pheno_pars$splineBuffer * 2)    
  splineEnd <- splineStart+(numDaysFit-1)
  
  #Repeat data for the buffer period of the final year to ensure that the spline gets tied down
  #Determine the most recent image, and then fill BEYOND that date until the spline buffer is reached
  #Filling with images from the previous year
  #########################################################################
  startFillDate <- max(dates)+1   #Find the most recent image, and then subtract a year 
  endFillDate <- as.Date(paste0(max(imgYrs)+1,'-01-01')) + pheno_pars$splineBuffer #This is the    
  
  if (startFillDate < endFillDate) {
    vals <- dates > (startFillDate-365) & dates <= (endFillDate-365)   #Get images from year before and repeat them
    b2_add <- b2[,vals]; b4_add <- b4[,vals]
    vi_add <- vi[,vals]
    snow_add <- snowPix[,vals]
    dates_add <- dates[vals] + 365
    
    b2 <- cbind(b2,b2_add); b4 <- cbind(b4,b4_add)
    vi <- cbind(vi,vi_add);
    snowPix <- cbind(snowPix,snow_add);
    dates <- c(dates,dates_add)
  }
  
  
  #Loop through each pixel and estimate phenometrics
  anomaly_mat <- matrix(NA,numPix,365*length(phenYrs))
  for (i in 1:numPix) {anomaly_mat[i,] <- DoAnomalyHLS(b2[i,], b4[i,], vi[i,],
                                                       snowPix[i,],dates, imgYrs, phenYrs, splineStart, splineEnd, numDaysFit, pheno_pars)}
  return(anomaly_mat)
}



print('----------------------------------------------------------------------------------------------')
print('Start processing')
print('')
print('')



# #Read in arguments
# args <- commandArgs(trailingOnly=T)
# print(args)
# tile <- args[1]
# jsonFile <- args[2]
# runLog <- args[3] 
# errorLog <- args[4] 

tile <- '16SEJ'
jsonFile <- "/projectnb/modislc/projects/landsat_sentinel/MSLSP_HLS30/16SEJ/parameters_2020_07_31_16_36_19.json"
runLog <- "/projectnb/modislc/projects/landsat_sentinel/MSLSP_HLS30/runLogs/16SEJ_instanceInfo_2020_07_31_16_36_19.txt"


#Get default parameters
params <- fromJSON(file=jsonFile)

#Pull paths for AWS or SCC?
if (params$setup$AWS_or_SCC == 'SCC') {
  params$setup$rFunctions <- params$SCC$rFunctions
  params$setup$productTable <- params$SCC$productTable
  params$setup$fmaskFunction <- params$SCC$fmaskFunction 
  params$setup$numCores <- params$SCC$numCores
  params$setup$numChunks <- params$SCC$numChunks
} else {
  params$setup$rFunctions <- params$AWS$rFunctions
  params$setup$productTable <- params$AWS$productTable
  params$setup$fmaskFunction <- params$AWS$fmaskFunction
  params$setup$numCores <- params$AWS$numCores
  params$setup$numChunks <- params$AWS$numChunks}


#Load functions
source(file=params$setup$rFunctions)

#Register the parallel backend
registerDoMC(cores=params$setup$numCores)

time_step1 <- 0
time_step2 <- 0

#Set up data
#####################################################

#Pull a few variables
imgStartYr <- params$setup$imgStartYr
imgEndYr <- params$setup$imgEndYr
phenStartYr <- params$setup$phenStartYr
phenEndYr <- params$setup$phenEndYr
numChunks <- params$setup$numChunks

#
params$phenology_parameters$dormStart <- as.Date(params$phenology_parameters$dormStart)
params$phenology_parameters$dormEnd <- as.Date(params$phenology_parameters$dormEnd)

#Sort product table
############################
productTable <- read.csv(params$setup$productTable,header=T,stringsAsFactors = F)

params$phenology_parameters$numLyrs <- sum(!is.na(productTable$calc_lyr))
params$phenology_parameters$numLyrsCycle2 <- length(grep('_2',productTable$short_name[!is.na(productTable$calc_lyr)]))

#Location of specific layers for which we will product outputs for ALL pixels (regardless of having phenology)
params$phenology_parameters$loc_numCycles <- productTable$calc_lyr[productTable$short_name == 'NumCycles']
params$phenology_parameters$loc_max <- productTable$calc_lyr[productTable$short_name == 'EVImax']
params$phenology_parameters$loc_amp <- productTable$calc_lyr[productTable$short_name == 'EVIamp']
params$phenology_parameters$loc_numObs <- productTable$calc_lyr[productTable$short_name == 'numObs']
params$phenology_parameters$loc_numObs_count_snow <- productTable$calc_lyr[productTable$short_name == 'numObsCountSnow']
params$phenology_parameters$loc_maxGap <- productTable$calc_lyr[productTable$short_name == 'maxGap']
params$phenology_parameters$loc_maxGap_count_snow <- productTable$calc_lyr[productTable$short_name == 'maxGapCountSnow']

#If we are running S10, set resolution to 10m. Otherwise, doing things at 30m
if (params$setup$AWS_or_SCC == "SCC" & params$SCC$runS10) {params$setup$image_res <- 10} else {params$setup$image_res <- 30}


#Sort years
######################
imgYrs <- imgStartYr:imgEndYr        #What years will we consider for time-series filling?
phenYrs <- phenStartYr:phenEndYr     #What years will we calculate phenology for?


#Get full image list
###########
imgList <- list.files(path=params$dirs$imgDir, pattern=glob2rx("HLS*.hdf"), full.names=T)

#Get the year and doy of each image. Restrict to time period of interest
##########################
sensor <- matrix(NA,length(imgList),1)
yrdoy = as.numeric(matrix(NA,length(imgList),1))
for (i in 1:length(imgList)) {
  imgName_strip = tail(unlist(strsplit(imgList[i],'/')),n = 1)
  sensor[i] <- unlist(strsplit(imgName_strip,'.',fixed = T))[2]
  yrdoy[i]= unlist(strsplit(imgName_strip,'.',fixed = T))[4]}
doys <- as.numeric(format(as.Date(strptime(yrdoy, format="%Y%j")),'%j'))
years <- as.numeric(format(as.Date(strptime(yrdoy, format="%Y%j")),'%Y'))

keep <- (years >= (imgStartYr - 1)) & (years <= (imgEndYr + 1))  #Only keep imagery +/- 1 (need 6 month buffer)

#Deterime which data to keep
if (!params$setup$includeLandsat) {drop <- sensor == 'L30'; keep[drop] <- FALSE}
if (!params$setup$includeSentinel) {drop <- sensor == 'S30'; keep[drop] <- FALSE}  

imgList <- imgList[keep]
yrdoy <- yrdoy[keep]
doys <- doys[keep]
years <- years[keep]
sensor <- sensor[keep]


uniqueYrs <- sort(unique(years))  #What years do we actually have imagery for? (imgYears +/- 1 for buffer years)


#Record number of images for each sensor per year
for (y in uniqueYrs) {
  if (params$setup$AWS_or_SCC == "SCC" & params$SCC$runS10){
    nS10 <- sum(years == y & sensor == 'S10')
    cat(paste0('S10_',y,':',nS10,'\n'), file=runLog, append=T)
  } else {
    nL30 <- sum(years == y & sensor == 'L30')
    nS30 <- sum(years == y & sensor == 'S30')
    cat(paste0('L30_',y,':',nL30,'\n'), file=runLog, append=T)
    cat(paste0('S30_',y,':',nS30,'\n'), file=runLog, append=T)
  }
}



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
numPixPerChunk <- chunkEnd - chunkStart + 1  #Number of pixels in each chunk


#Read in water mask
water <- readGDAL(paste0(params$dirs$imgDir,'water_',tile,'.tif'),silent=T)$band1
waterMask <- water == 2 | water == 0    #Mask water and zero (zero = ocean far from shore)
remove(water)

#####################################################
#####################################################
#Run Anomaly code for each image chunk
# Read in arguments
args <- commandArgs()
print(args)

j <- as.numeric(args[3])

AnoChunk <- runAnomalyChunk(j, numPixPerChunk[j], imgYrs, phenYrs, errorLog, params)

setwd('/projectnb/modislc/users/mkmoon/MuSLI/VI_anomaly/16SEJ/')  
cn <- sprintf('%03d',j)
save(AnoChunk,file=paste(tile,'_',cn,'.rda',sep=''))

  

###############################################
path <- '/projectnb/modislc/users/mkmoon/MuSLI/VI_anomaly/16SEJ/'
sstr <- '16SEJ*.rda'
file <- list.files(path, pattern=glob2rx(sstr),full.names=T)
for(i in 1:length(file)){
  load(file[i])
  if(i==1){
    temp <- apply(AnoChunk[,214:305],1,mean,na.rm=T)
  }else{
    temp <- c(temp,apply(AnoChunk[,214:305],1,mean,na.rm=T))
  }
  print(i)
}

aa <- setValues(baseImage,temp)
plot(aa)

hist(temp,xlim=c(-0.07,0.07),breaks=seq(-10,10,0.005))

# setwd('/projectnb/modislc/users/mkmoon/MuSLI/VI_anomaly/18TYN/')
# writeRaster(aa,filename="test_2016.tif", format="GTiff", overwrite=TRUE)

lct <- raster('/projectnb/modislc/data/lc_database/regional/united_states/nlcd_2011_hls_tiles/18TYN_nlcd_2011.tif')
lct_val <- values(lct)

aa_f <- temp
aa_c <- temp

aa_f[which(lct_val!=41&lct_val!=42)] <- NA
aa_c[which(lct_val!=81&lct_val!=82&lct_val!=23&lct_val!=24)] <- NA

##################
library(RColorBrewer)

aa[aa>  0.07] <-  0.069999
aa[aa< -0.07] <- -0.069999

Pal <- colorRampPalette(c('darkmagenta','white','darkgreen'))
Col <- Pal(11)

lct_sub <- crop(lct,extent(705273,725700,4760000,4780000))
aa_sub <- crop(aa,extent(705273,725700,4760000,4780000))


setwd('/projectnb/modislc/users/mkmoon/MuSLI/figure/')
png(filename='anomaly_hist.png',width=12,height=5.5,units='in',res=600)

par(fig=c(0,0.5,0,1),oma=c(0,0,0,0),mar=c(0.5,0.5,0.5,0.5),mgp=c(3.5,1,0))
# plot(lct_sub)
plot(aa_sub,col=Col,axes=F,box=F,legend=F,colNA='grey50')

par(fig=c(0.5,1,0,1),oma=c(1,1,1,1),mar=c(3,3,1,1),mgp=c(3.5,1,0),new=T)
hist(temp,xlim=c(-0.1,0.1),breaks=seq(-10,10,0.005),
     main='EVI2 anomaly in 2016 summer',
     cex.main=1.5,cex.axis=1.2)
hist(aa_f,xlim=c(-0.1,0.1),breaks=seq(-10,10,0.005),col=rgb(0,100/256,0,0.9),add=T)
hist(aa_c,xlim=c(-0.1,0.1),breaks=seq(-10,10,0.005),col=rgb(139/256,0,139/256,0.9),add=T)

dev.off()

png(filename='anomaly_nlcd.png',width=6,height=5.5,units='in',res=600)
par(fig=c(0,0.5,0,1),oma=c(0,0,0,0),mar=c(0.5,0.5,0.5,0.5),mgp=c(3.5,1,0))
plot(lct_sub)
dev.off()

#############################
# Flux tower data
geog_crs = CRS("+proj=longlat +datum=WGS84")
site <- data.frame(1,-86.4131,39.3232)
colnames(site) <- c('id','lon','lat')
xy   <- site[,c(2,3)]
bb   <- SpatialPointsDataFrame(coords=xy,data=site,proj4string=geog_crs)
bb <- spTransform(bb,crs(baseImage))
plot(bb,add=T)

setwd('/projectnb/modislc/projects/landsat_sentinel/MSLSP_assessment/shps/')
writeOGR(bb, ".",paste(tile,'_pts',sep=''), driver="ESRI Shapefile",overwrite=T)

# # Load shapefiles
# shp_hls <- shapefile('/projectnb/modislc/projects/landsat_sentinel/shapefiles/sentinel2_tiles_world/sentinel2_tiles_world/sentinel2_tiles_north_america_Albers.shp')
# cpt     <- spTransform(bb,crs(shp_hls))
# shp_hls <- intersect(shp_hls,cpt)

load('/projectnb/modislc/users/mkmoon/autumn/sample_data/mms_gppFilled.rda')

gpp_day <- matrix(gpp,48,length(gpp)/48)
gpp_day <- apply(gpp_day,2,mean)*0.0864*44/2

# setwd('/projectnb/modislc/users/mkmoon/MuSLI/VI_anomaly/16SEJ/figure/')
# png(filename='gpp_daily.png',width=10,height=5,units='in',res=300)
par(oma=c(2,2,1,1),mar=c(4,5,4,4),mgp=c(2.5,1,0))
plot(gpp_day,xlab='Dates',type='o',
     ylab=expression(paste('g CO'[2],' ',m^-2,' ',day^-1)),
     ylim=c(-5,32),
     cex.lab=1.5,xaxt='n',cex.axis=1.5)
axis(1,at=c(1,366,366*2,366*3),c(2016,2017,2018,2019),cex.axis=1.5)
# title('US-MMS daily GPP',cex=2)
# dev.off()

path <- '/projectnb/modislc/projects/landsat_sentinel/MSLSP_HLS30/16SEJ/imageChunks/c1'
sstr <- '*16SEJ*.Rds'
file <- list.files(path, pattern=glob2rx(sstr),full.names=T)
load(file[100])
aa <- setValues(baseImage,temp)