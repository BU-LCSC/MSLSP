# Time-series and diagnostic functions
###############################




#---------------------------------------------------------------------
#Calculate pheno metrics for each pixel
#This version using alternate years to gap fill
#Code adapted by Douglas Bolton
#Based on MODIS C6 algorithm developed by Josh Gray
#---------------------------------------------------------------------
DoPhenologyHLS_Diagnostic <- function(b2, b3, b4, b5, b6, b7, vi, snowPix, dates, imgYrs, phenYrs, splineStart, splineEnd, numDaysFit, pheno_pars){
  
  vi_orig <- vi

  #Despike, calculate dormant value, fill snow with dormant value, despike poorly fit snow values
  log <- try({
    #Despike
    spikes <- CheckSpike_MultiBand(b2/10000,b4/10000,vi, dates, pheno_pars)
    
    #Set despiked values to NA
    vi[spikes] <- NA
    b2[spikes] <- NA; b3[spikes] <- NA; b4[spikes] <- NA
    b5[spikes] <- NA; b6[spikes] <- NA; b7[spikes] <- NA
    
    dormIms <- dates >= pheno_pars$dormStart & dates <= pheno_pars$dormEnd
    
    vi_dorm <- quantile(vi[dormIms],probs=pheno_pars$dormantQuantile,na.rm=T)   #Calc vi dormant value
    snowPix <- Screen_SnowFills(vi,vi_dorm,snowPix,dates,pheno_pars)              #Screen poorly filled snow values
    
    
    #now calculate dormancy values and fill individual bands
    dormObs <- dormIms & vi < vi_dorm    #Defining dormant observations for bands as median on dates when vi < vi_dorm
    b2_dorm <- median(b2[dormObs], na.rm=T); b2[snowPix] <- b2_dorm
    b3_dorm <- median(b3[dormObs], na.rm=T); b3[snowPix] <- b3_dorm
    b4_dorm <- median(b4[dormObs], na.rm=T); b4[snowPix] <- b4_dorm
    b5_dorm <- median(b5[dormObs], na.rm=T); b5[snowPix] <- b5_dorm
    b6_dorm <- median(b6[dormObs], na.rm=T); b6[snowPix] <- b6_dorm
    b7_dorm <- median(b7[dormObs], na.rm=T); b7[snowPix] <- b7_dorm
    
    vi[snowPix] <- vi_dorm   #Fill remaining snow values with dormant value for vi  
    
    
    outList <- list()  #Where to store outputs
    outList$original_VI <- vi_orig
    outList$filled_VI <- vi
    outList$snowPix <- snowPix
    outList$dates <- dates
    
    
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
    b2Mat <- smoothMat; b3Mat <- smoothMat; b4Mat <- smoothMat
    b5Mat <- smoothMat; b6Mat <- smoothMat; b7Mat <- smoothMat
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
        b2Mat[fillDs,y] <- b2[dateRange]; b3Mat[fillDs,y] <- b3[dateRange]; b4Mat[fillDs,y] <- b4[dateRange]
        b5Mat[fillDs,y] <- b5[dateRange]; b6Mat[fillDs,y] <- b6[dateRange]; b7Mat[fillDs,y] <- b7[dateRange]
        
      },silent=TRUE)
    }
    
    
    xs <- rep(daysVec,numYrs)
    ys <- matrix(fillMat,vecLength)
    ysGood <- !is.na(ys)
    baseW <- matrix(baseWeights,vecLength)   #Base Weights are 1=clear observation, 0.5=snow-filled
    
    b2Mat <- matrix(b2Mat,vecLength);b3Mat <- matrix(b3Mat,vecLength)
    b4Mat <- matrix(b4Mat,vecLength);b5Mat <- matrix(b5Mat,vecLength) 
    b6Mat <- matrix(b6Mat,vecLength);b7Mat <- matrix(b7Mat,vecLength)  
    
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
  if(inherits(log, "try-error")){return(NA)}   
  
  
  outAll=c()
  for (y in yToDo) {
    log <- try({
      
      yearList <- list()
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
      

      
      yearList$smoothed_vi <- smoothed_vi
      yearList$smoothed_vi_nofilling <- smoothMat[,y]
      yearList$smoothed_dates <- pred_dates
      yearList$smoothed_day <- daysVec
      yearList$filled_vi <- ys[theInds]
      yearList$filled_dates <- rep(pred_dates,length(yrs))[theInds]   #repeating the same dates for all years
      yearList$filled_day <- xs_sub
      yearList$filled_weigth <- w_sub
      yearList$filled_baseWeight <- baseW[theInds]

      
      
      # Fit spline to individual bands. Waiting until now so that we don't fit for pixels with no cycle
      yearList$bands$smooth_b2 <- Smooth_Bands(b2Mat[theInds], xs_sub, daysVec, w_sub, pheno_pars)
      yearList$bands$smooth_b3 <- Smooth_Bands(b3Mat[theInds], xs_sub, daysVec, w_sub, pheno_pars)
      yearList$bands$smooth_b4 <- Smooth_Bands(b4Mat[theInds], xs_sub, daysVec, w_sub, pheno_pars)
      yearList$bands$smooth_b5 <- Smooth_Bands(b5Mat[theInds], xs_sub, daysVec, w_sub, pheno_pars)
      yearList$bands$smooth_b6 <- Smooth_Bands(b6Mat[theInds], xs_sub, daysVec, w_sub, pheno_pars)
      yearList$bands$smooth_b7 <- Smooth_Bands(b7Mat[theInds], xs_sub, daysVec, w_sub, pheno_pars)
      
      yearList$bands$filled_b2 <- b2Mat[theInds]
      yearList$bands$filled_b3 <- b3Mat[theInds]
      yearList$bands$filled_b4 <- b4Mat[theInds]
      yearList$bands$filled_b5 <- b5Mat[theInds]
      yearList$bands$filled_b6 <- b6Mat[theInds]
      yearList$bands$filled_b7 <- b7Mat[theInds]
      
      
      
      #Fit phenology
      peaks <- FindPeaks(smoothed_vi)
      if (all(is.na(peaks))) {outList[[paste0('y',yrs[y])]] <- yearList;next}   #If no peaks, calc annual metrics, and move to next year
      
      #Find full segments
      full_segs <- GetSegs(peaks, smoothed_vi, pheno_pars)
      if (is.null(full_segs)) {outList[[paste0('y',yrs[y])]] <- yearList;next}  #If no valid segments, calc annual metrics, and move to next year
      
      #Only keep segments with peaks within year *****
      full_segs <- full_segs[inYear[sapply(full_segs, "[[", 2)] ]  #check if peaks are in the year
      
      #Get PhenoDates
      pheno_dates <- GetPhenoDates(full_segs, smoothed_vi, pred_dates, pheno_pars)
      phen <- unlist(pheno_dates, use.names=F)
      if (all(is.na(phen))) {outList[[paste0('y',yrs[y])]] <- yearList;next} #If no dates detected, calc annual metrics, and move to next year
      
      yearList$phenDates <- phen
      
      #Get metrics that describe the segments and the year
      
      
      #First, get metrics counting gap filled observations as "good" observations
      seg_metricsFill <- lapply(full_segs, GetSegMetricsLight, daysVec, sort(xs_sub))
      un <- unlist(seg_metricsFill, use.names=F)
      ln <- length(un)      
      yearList$gup_maxgap_frac_filled <- un[seq(1, ln, by=2)] * 100
      yearList$gdown_maxgap_frac_filled <- un[seq(2, ln, by=2)] * 100
      
      
      #Second, get segment metrics with snow observations counted as "good" observations
      filled_vi <- fillMat[,y]
      seg_metricsFill <- lapply(full_segs, GetSegMetricsLight, daysVec, daysVec[!is.na(filled_vi)])
      un <- unlist(seg_metricsFill, use.names=F)
      ln <- length(un)      
      yearList$gup_maxgap_frac_count_snow <- un[seq(1, ln, by=2)] * 100
      yearList$gdown_maxgap_frac_count_snow <- un[seq(2, ln, by=2)] * 100
      
      #And get calendar year metrics with snow counted as good
      yearList$numObs_count_snow <- sum(!is.na(filled_vi) & inYear)
      yearList$maxGap_annual_count_snow <- max(diff(c( pheno_pars$splineBuffer+1, daysVec[!is.na(filled_vi) & inYear], 365+pheno_pars$splineBuffer))) 
      
      
      #Now get the full segment metrics, not counting snow and not counting gap filled
      filled_vi[baseWeights[,y] < 1] <- NA    #If weight is less than 1, implies it is a snow-fill, and we don't want to count snow-filled as a valid observation. So set to NA.
      yearList$numObs <- sum(!is.na(filled_vi) & inYear)   #Number of observations in year
      yearList$maxGap_annual <- max(diff(c( pheno_pars$splineBuffer+1, daysVec[!is.na(filled_vi) & inYear], 365+pheno_pars$splineBuffer)))  #Max gap (in days) during year
      
      seg_metrics <- lapply(full_segs, GetSegMetrics, smoothed_vi, filled_vi[!is.na(filled_vi)], pred_dates, pred_dates[!is.na(filled_vi)]) #full segment metrics
      
      
      #Unlist and scale the seg metrics
      un <- unlist(seg_metrics, use.names=F)
      ln <- length(un)
      yearList$stats$seg_amp <- un[seq(1, ln, by=9)] * 10000
      yearList$stats$seg_max <- un[seq(2, ln, by=9)] * 10000
      yearList$stats$seg_int <- un[seq(3, ln, by=9)] * 100
      yearList$stats$gup_rsq <- un[seq(4, ln, by=9)] * 10000
      yearList$stats$gup_numObs <- un[seq(5, ln, by=9)]
      yearList$stats$gup_maxgap_frac <- un[seq(6, ln, by=9)] * 100
      yearList$stats$gdown_rsq <- un[seq(7, ln, by=9)] * 10000
      yearList$stats$gdown_numObs <- un[seq(8, ln, by=9)]
      yearList$stats$gdown_maxgap_frac <- un[seq(9, ln, by=9)] * 100
      
      yearList$stats$numRecords <- length(seg_amp)  #how many cycles were recorded
      
      naCheck <- is.na(seg_amp)
      yearList$stats$numCyc <- sum(naCheck == 0)  #how many cycles have good data (seg metrics has valid observations)
    

      
      if (numRecords == 1) {
        yearList$bands$composites <- as.vector(MakeComposites(phen,pred_dates,yearList$bands$smooth_b2, 
                                                              yearList$bands$smooth_b3, yearList$bands$smooth_b4, yearList$bands$smooth_b5, yearList$bands$smooth_b6, yearList$bands$smooth_b7))
   
        
      } else {
        #If there are multiple cycles, sort by amplitude and report two highest amplitudes (highest amplitude first)
        theOrd <- order(seg_amp,decreasing=T)   
        
        phen1 <- phen[seq(theOrd[1], length(phen), by = numRecords)]
        phen2 <- phen[seq(theOrd[2], length(phen), by = numRecords)]
        
        yearList$bands$composites <- as.vector(MakeComposites(phen1,pred_dates,yearList$bands$smooth_b2, 
                                                              yearList$bands$smooth_b3, yearList$bands$smooth_b4, yearList$bands$smooth_b5, yearList$bands$smooth_b6, yearList$bands$smooth_b7))
        yearList$bands$composites_2 <- as.vector(MakeComposites(phen2,pred_dates,yearList$bands$smooth_b2, 
                                                                yearList$bands$smooth_b3, yearList$bands$smooth_b4, yearList$bands$smooth_b5, yearList$bands$smooth_b6, yearList$bands$smooth_b7))
        
      }
      outList[[paste0('y',yrs[y])]] <- yearList
      
    },silent=TRUE)  #End of the try block
    if(inherits(log, "try-error")){outList[[paste0('y',yrs[y])]] <- yearList} 
  }
  
  return(outList)
}




#---------------------------------------------------------------------
#Run Phenology code for specific pixels in an image chunk
#Douglas Bolton
#---------------------------------------------------------------------
runPhenoChunk_Diagnostic <- function(chunk, pixels, ids, tile, numPix, chunkDir, imgYrs, phenYrs, params, codeVersion) {
  
  pheno_pars <- params$phenology_parameters
  
  chunkFold <- paste0(chunkDir,'c',chunk,'/') 
  
  #Get all images to process
  ######################
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
    if (inherits(imgData, 'try-error')) {next} 
    
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
  b2 <- bands[[1]];  b3 <- bands[[2]]; b4 <- bands[[3]]
  b5 <- bands[[4]]; b6 <- bands[[5]]; b7 <- bands[[6]];
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
    b2_add <- b2[,vals]; b3_add <- b3[,vals]; b4_add <- b4[,vals]
    b5_add <- b5[,vals]; b6_add <- b6[,vals]; b7_add <- b7[,vals]
    vi_add <- vi[,vals]
    snow_add <- snowPix[,vals]
    dates_add <- dates[vals] + 365
    
    b2 <- cbind(b2,b2_add); b3 <- cbind(b3,b3_add); b4 <- cbind(b4,b4_add)
    b5 <- cbind(b5,b5_add); b6 <- cbind(b6,b6_add); b7 <- cbind(b7,b7_add);
    vi <- cbind(vi,vi_add);
    snowPix <- cbind(snowPix,snow_add);
    dates <- c(dates,dates_add)
  }
  
  if (codeVersion == 'V1') {
  #Loop through each pixel and add to list
  result <- vector("list", length(pixels))
  for (i in 1:length(pixels)) {
    pix <- pixels[i]
    result[[i]] <- DoPhenologyHLS_Diagnostic(b2[pix,],  b3[pix,],  b4[pix,],  b5[pix,],  b6[pix,], b7[pix,],  vi[pix,],
                                                  snowPix[pix,],dates, imgYrs, phenYrs, splineStart, splineEnd, numDaysFit, pheno_pars)
  }
  } else {
    #Loop through each pixel and add to list
    result <- vector("list", length(pixels))
    for (i in 1:length(pixels)) {
      pix <- pixels[i]
      result[[i]] <- DoPhenologyHLS_Diagnostic_V0(b2[pix,],  b3[pix,],  b4[pix,],  b5[pix,],  b6[pix,], b7[pix,],  vi[pix,],
                                               snowPix[pix,],dates, imgYrs, splineStart, numDaysFit, pheno_pars)
    }
  }
  
  names(result) <- paste0(tile,'_',ids)
  return(result)
}









#---------------------------------------------------------------------
#Extract time series and run phenology code for point locations
#Douglas Bolton
#---------------------------------------------------------------------
Extract_Timeseries <- function(tile, imgDir, chunkDir, imgYrs, phenYrs, numCores, params, codeVersion='V1', pixelNumbers=NULL, shpName=NULL)  {
  
  
  #Sort out image chunks
  ######################
  
  #Determine number of image chunks based on how many chunk folders there are
  numChunks <- max(as.numeric(gsub("c","",(list.files(chunkDir, pattern=glob2rx("c*"), recursive=FALSE)))))
  
  #Use HLS image to get tile info. Make image of pixel ids
  imgList <- list.files(path=imgDir, pattern=glob2rx("HLS*.hdf"), full.names=T)
  baseImage <- raster(paste0('HDF4_EOS:EOS_GRID:',imgList[1], ':Grid:QA'))
  

  nrows <- dim(baseImage)[1]
  ncols <- dim(baseImage)[2]
  numPix <- nrows*ncols
  
  if (is.null(pixelNumbers)) {
    #Read in shapefile, reporject to match baseImage
    shp <- readOGR(shpName)
    shp$id <- as.numeric(paste0(shp$id))
    shp <- spTransform(shp, crs(baseImage))

    pixID <- 1:numPix
    pixRast <-setValues(baseImage,pixID)
    pixelNumbers <- extract(pixRast,shp)
  } else {shp <- list();shp$id=1:length(pixelNumbers)}
  
  #Sort out chunk boundaries
  boundaries <- chunkBoundaries(numPix,numChunks)
  chunkStart <- boundaries[[1]]
  chunkEnd <- boundaries[[2]]
  chunks <- 1:numChunks
  numPixPerChunk <- chunkEnd - chunkStart + 1  #Number of pixels in each chunk
  
  
  #Loop through points. Determine chunk # and pixel #
  theChunks <- matrix(0,length(shp$id))
  thePixels <- matrix(0,length(shp$id))
  theIDs <- matrix(0,length(shp$id))
  for (i in 1:length(shp$id)) {
    id <- shp$id[i]            #Get id of point (Shapefile needs id column!)
    
    pixelNumber <- pixelNumbers[shp$id==id]
    
    inds <- which(chunkEnd > pixelNumber)
    theChunk <- chunks[min(inds)]
    
    
    chunkPix <- chunkStart[theChunk]:chunkEnd[theChunk]
    
    thePixels[i] <- which(chunkPix == pixelNumber)  
    theChunks[i] <- theChunk
    theIDs[i] <- id
    
  }
  
  
  uniqChunks <- unique(theChunks)
  
  
  results <- foreach(i=1:length(uniqChunks),.combine='c',.multicombine = TRUE)  %dopar% {
    chunk <- uniqChunks[i]
    pixels <- thePixels[theChunks == chunk]
    ids <- theIDs[theChunks == chunk]
    numPix <- numPixPerChunk[chunk]
    
    out <- runPhenoChunk_Diagnostic(chunk, pixels, ids, tile, numPix, chunkDir, imgYrs, phenYrs, params,codeVersion=codeVersion)
    list(out)
  }
  return(unlist(results,recursive=F))
}
  
  
  
  
  
  






































#---------------------------------------------------------------------
#This is the originial version of the code V0
#Make sure the V0 functions are loaded when running this to get the true V0 result!
#---------------------------------------------------------------------
DoPhenologyHLS_Diagnostic_V0 <- function(b2, b3, b4, b5, b6, b7, vi, snowPix, dates, yrs, splineStart, numDaysFit, pheno_pars) {
  
  vi_orig <- vi
  
  numYrs <- length(yrs)
  daysVec <- 1:numDaysFit
  vecLength <- numDaysFit*numYrs
  
  
  #Despike, calculate dormant value, fill snow with dormant value, despike poorly fit snow values
  log <- try({
    #Despike
    spikes <- CheckSpike_MultiBand(b2/10000,b4/10000,vi, dates, pheno_pars)
    
    #Set despiked values to NA
    vi[spikes] <- NA
    b2[spikes] <- NA; b3[spikes] <- NA; b4[spikes] <- NA
    b5[spikes] <- NA; b6[spikes] <- NA; b7[spikes] <- NA
    
    dormIms <- dates >= pheno_pars$dormStart & dates <= pheno_pars$dormEnd
    
    vi_dorm <- quantile(vi[dormIms],probs=pheno_pars$dormantQuantile,na.rm=T)   #Calc vi dormant value
    snowPix <- Screen_SnowFills(vi,vi_dorm,snowPix,dates)              #Screen poorly filled snow values
    vi[snowPix] <- vi_dorm                                             #Fill remaining snow values with dormant value  
    
    #now calculate dormancy values and fill each band
    b2_dorm <- quantile(b2[dormIms],probs=pheno_pars$dormantQuantile,na.rm=T); b2[snowPix] <- b2_dorm
    b3_dorm <- quantile(b3[dormIms],probs=pheno_pars$dormantQuantile,na.rm=T); b3[snowPix] <- b3_dorm
    b4_dorm <- quantile(b4[dormIms],probs=pheno_pars$dormantQuantile,na.rm=T); b4[snowPix] <- b4_dorm
    b5_dorm <- quantile(b5[dormIms],probs=pheno_pars$dormantQuantile,na.rm=T); b5[snowPix] <- b5_dorm
    b6_dorm <- quantile(b6[dormIms],probs=pheno_pars$dormantQuantile,na.rm=T); b6[snowPix] <- b6_dorm
    b7_dorm <- quantile(b7[dormIms],probs=pheno_pars$dormantQuantile,na.rm=T); b7[snowPix] <- b7_dorm
    

    
    outList <- list()  #Where to store outputs
    outList$original_VI <- vi_orig
    outList$filled_VI <- vi
    outList$snowPix <- snowPix
    outList$dates <- dates
    
    

    
    #First, we will fit splines to each year invidually
    #To line up observations from each year, we will create a matrix for vi and each band (numDaysFit x numYears)
    
    #Assign base Weights
    baseW <- matrix(1,length(vi))
    baseW[snowPix == 1] <- pheno_pars$snowWeight
    
    #Option to not include observations from buffer period
    #climInds <- matrix(TRUE,numDaysFit)
    #climInds[pheno_pars$splineBuffer:(numDaysFit-pheno_pars$splineBuffer)] <- FALSE
    
    splineEnd <- splineStart+(numDaysFit-1)
    
    #Single spline through all years
    fullDates <- seq(splineStart[1], splineEnd[numYrs], by="day")
    smoothed <- Just_Smooth(vi, dates, fullDates, baseW, vi_dorm, pheno_pars)
    
    gDates <- dates[!is.na(vi)]
    
    #Mask spline where gap greater than 45 days
    dDiff <- diff(gDates) > pheno_pars$maxGapForSplineComparison
    dStart <- gDates[c(dDiff,FALSE)]
    dEnd <- gDates[c(FALSE,dDiff)]
    for (d in 1:length(dStart)) {smoothed[fullDates >= dStart[d] & fullDates < dEnd[d]] <- NA}
    
    #Mask spline before first image and after last image
    smoothed[fullDates < gDates[1]] <- NA
    smoothed[fullDates > gDates[length(gDates)]] <- NA
    
    smoothMat <- matrix(NA, numDaysFit, numYrs)
    fillMat <- smoothMat
    b2Mat <- smoothMat; b3Mat <- smoothMat; b4Mat <- smoothMat
    b5Mat <- smoothMat; b6Mat <- smoothMat; b7Mat <- smoothMat
    baseWeights <- smoothMat
    for (y in 1:numYrs) {
      dateRange <- dates >= splineStart[y] & dates <= splineEnd[y] & !is.na(vi)   
      fillDs <- as.numeric(dates[dateRange] - splineStart[y]) + 1
      
      smoothMat[,y] <- smoothed[fullDates >= splineStart[y] & fullDates <= splineEnd[y]]
      baseWeights[fillDs,y] <- baseW[dateRange]
      fillMat[fillDs,y] <- vi[dateRange]
      b2Mat[fillDs,y] <- b2[dateRange]; b3Mat[fillDs,y] <- b3[dateRange]; b4Mat[fillDs,y] <- b4[dateRange]
      b5Mat[fillDs,y] <- b5[dateRange]; b6Mat[fillDs,y] <- b6[dateRange]; b7Mat[fillDs,y] <- b7[dateRange] }
    
    
    xs <- rep(daysVec,numYrs)
    ys <- matrix(fillMat,vecLength)
    ysGood <- !is.na(ys)
    baseW <- matrix(baseWeights,vecLength)   #Base Weights are 1=clear observation, 0.5=snow-filled
    
    b2Mat <- matrix(b2Mat,vecLength); b3Mat <- matrix(b3Mat,vecLength) 
    b4Mat <- matrix(b4Mat,vecLength); b5Mat <- matrix(b5Mat,vecLength) 
    b6Mat <- matrix(b6Mat,vecLength); b7Mat <- matrix(b7Mat,vecLength) 
    
    

    
    
  },silent=TRUE)
  #If there is an error despiking or other initial steps, return NAs
  if(inherits(log, "try-error")){return(NA)}   
  
  
  outAll=c()
  for (y in 1:length(yrs)) {
    log <- try({
      
      yearList <- list()
      
      pred_dates <- seq(splineStart[y], splineEnd[y], by="day")  
      #Only compare on dates that have splined data in target year
      ind <- !is.na(smoothMat[,y])
      sub_vi <- smoothMat[ind,]
      numGoodDays <- colSums(!is.na(sub_vi))  #how many days actually have splined data?
      
      #What approach to use for weighting
      ######
      #Calculate euclidean distance
      eucl <- colSums((sub_vi - sub_vi[,y])^2,na.rm=T)^0.5 
      
      #Now calculate euclidean distance assuming the average through the year
      #Scale euculidean distances between this value and a perfect fit (0 to 1)
      theAvg <- matrix(mean(sub_vi[,y],na.rm=T),length(sub_vi[,y]),numYrs,byrow=T)
      theAvg[is.na(sub_vi)] <- NA   #only calculate for days that have data
      
      max_eucl <- colSums((theAvg - sub_vi[,y])^2,na.rm=T)^0.5    #calculate eucidean distance for this case
      scaled_eucl <- 1 - (eucl / max_eucl)
      scaled_eucl[scaled_eucl < 0] <- 0
      
      #Weigh as the scaled euclidean distance (0 = same/worse than assuming average, 1 = perfect fit)
      weight <-pheno_pars$maxWeight * scaled_eucl
      weight[numGoodDays < pheno_pars$minDaysForSplineComparison]  <- 0
      weight[is.na(weight)] <- 0
      weight[is.infinite(weight)] <- 0
      weight[y] <- 1  #Set weight of target year to 1
      weights <- matrix(weight,numDaysFit,numYrs,byrow=T)
      
      #weights[climInds,yrs!=yrs[y]] <- 0  #Option to only include WITHIN year observations for climatology (and not the buffer period)
      
      #Now that we have weights, calculate phenology
      #######################
      weights <- matrix(weights,vecLength) * baseW   #Multiple weights by base weight (1=good,0.5=snow-filled)
      theInds <- ysGood & weights > 0
      xs_sub <- xs[theInds]; w_sub <- weights[theInds]
      smooth_vi <- Just_Smooth(ys[theInds], xs_sub, daysVec, w_sub, vi_dorm, pheno_pars)  #Fit spline
      
      
      yearList$smoothed_vi <- smooth_vi
      yearList$smoothed_vi_nofilling <- smoothMat[,y]
      yearList$smoothed_dates <- pred_dates
      yearList$smoothed_day <- daysVec
      yearList$filled_vi <- ys[theInds]
      yearList$filled_dates <- rep(pred_dates,length(yrs))[theInds]   #repeating the same dates for all years
      yearList$filled_day <- xs_sub
      yearList$filled_weigth <- w_sub
      yearList$filled_baseWeight <- baseW[theInds]
      
      
      
      # Fit spline to individual bands. Waiting until now so that we don't fit for pixels with no cycle
      yearList$bands$smooth_b2 <- Just_Smooth(b2Mat[theInds], xs_sub, daysVec, w_sub, b2_dorm, pheno_pars)
      yearList$bands$smooth_b3 <- Just_Smooth(b3Mat[theInds], xs_sub, daysVec, w_sub, b3_dorm, pheno_pars)
      yearList$bands$smooth_b4 <- Just_Smooth(b4Mat[theInds], xs_sub, daysVec, w_sub, b4_dorm, pheno_pars)
      yearList$bands$smooth_b5 <- Just_Smooth(b5Mat[theInds], xs_sub, daysVec, w_sub, b5_dorm, pheno_pars)
      yearList$bands$smooth_b6 <- Just_Smooth(b6Mat[theInds], xs_sub, daysVec, w_sub, b6_dorm, pheno_pars)
      yearList$bands$smooth_b7 <- Just_Smooth(b7Mat[theInds], xs_sub, daysVec, w_sub, b7_dorm, pheno_pars)
      
      yearList$bands$filled_b2 <- b2Mat[theInds]
      yearList$bands$filled_b3 <- b3Mat[theInds]
      yearList$bands$filled_b4 <- b4Mat[theInds]
      yearList$bands$filled_b5 <- b5Mat[theInds]
      yearList$bands$filled_b6 <- b6Mat[theInds]
      yearList$bands$filled_b7 <- b7Mat[theInds]
      
      
      
      #Fit phenology
      peaks <- FindPeaks(smooth_vi)
      if (all(is.na(peaks))) {outList[[paste0('y',yrs[y])]] <- yearList;next}   #If no peaks, calc annual metrics, and move to next year
      
      #Find full segments
      full_segs <- GetSegs(peaks, smooth_vi, pheno_pars)
      if (is.null(full_segs)) {outList[[paste0('y',yrs[y])]] <- yearList;next}  #If no valid segments, calc annual metrics, and move to next year
      
      #Only keep segments with peaks within year *****
      kept_peaks <- unlist(full_segs, use.names=F)  
      kept_peaks  <- kept_peaks[seq(2, length(kept_peaks), by=3)]    #extract peaks. If there are no peaks, an error will thrown, and will move to next year
      full_segs <- full_segs[as.numeric(format(pred_dates[kept_peaks],'%Y')) == yrs[y]]  #check if peaks are in the year
      
      #Get PhenoDates
      pheno_dates <- GetPhenoDates(full_segs, smooth_vi, pred_dates, pheno_pars)
      phen <- unlist(pheno_dates, use.names=F)
      if (all(is.na(phen))) {outList[[paste0('y',yrs[y])]] <- yearList;next} #If no dates detected, calc annual metrics, and move to next year
      
      yearList$phenDates <- phen
      
      #Get metrics that describe the segments and the year
      
      
      #Get the segment metrics
      filled_vi <- fillMat[,y]
      filled_vi[baseWeights[,y] < 1] <- NA    #If weight is less than 1, implies it is a snow-fill, and we don't want to count snow-filled as a valid observation. So set to NA.
      seg_metrics <- lapply(full_segs, GetSegMetrics, smooth_vi, filled_vi[!is.na(filled_vi)], pred_dates, pred_dates[!is.na(filled_vi)])
      
      
      #Unlist and scale the seg metrics
      un <- unlist(seg_metrics, use.names=F)
      ln <- length(un)
      yearList$stats$seg_amp <- un[seq(1, ln, by=9)] * 10000
      yearList$stats$seg_max <- un[seq(2, ln, by=9)] * 10000
      yearList$stats$seg_int <- un[seq(3, ln, by=9)] * 100
      yearList$stats$gup_rsq <- un[seq(4, ln, by=9)] * 10000
      yearList$stats$gup_numObs <- un[seq(5, ln, by=9)]
      yearList$stats$gup_maxgap_frac <- un[seq(6, ln, by=9)] * 100
      yearList$stats$gdown_rsq <- un[seq(7, ln, by=9)] * 10000
      yearList$stats$gdown_numObs <- un[seq(8, ln, by=9)]
      yearList$stats$gdown_maxgap_frac <- un[seq(9, ln, by=9)] * 100
      
      yearList$stats$numRecords <- length(seg_amp)  #how many cycles were recorded
      
      naCheck <- is.na(seg_amp)
      yearList$stats$numCyc <- sum(naCheck == 0)  #how many cycles have good data (seg metrics has valid observations)
      
      
      
      if (numRecords == 1) {
        yearList$bands$composites <- as.vector(MakeComposites(phen,pred_dates,yearList$bands$smooth_b2, 
                                                              yearList$bands$smooth_b3, yearList$bands$smooth_b4, yearList$bands$smooth_b5, yearList$bands$smooth_b6, yearList$bands$smooth_b7))
        
        
      } else {
        #If there are multiple cycles, sort by amplitude and report two highest amplitudes (highest amplitude first)
        theOrd <- order(seg_amp,decreasing=T)   
        
        phen1 <- phen[seq(theOrd[1], length(phen), by = numRecords)]
        phen2 <- phen[seq(theOrd[2], length(phen), by = numRecords)]
        
        yearList$bands$composites <- as.vector(MakeComposites(phen1,pred_dates,yearList$bands$smooth_b2, 
                                                              yearList$bands$smooth_b3, yearList$bands$smooth_b4, yearList$bands$smooth_b5, yearList$bands$smooth_b6, yearList$bands$smooth_b7))
        yearList$bands$composites_2 <- as.vector(MakeComposites(phen2,pred_dates,yearList$bands$smooth_b2, 
                                                                yearList$bands$smooth_b3, yearList$bands$smooth_b4, yearList$bands$smooth_b5, yearList$bands$smooth_b6, yearList$bands$smooth_b7))
        
      }
      outList[[paste0('y',yrs[y])]] <- yearList
      
    },silent=TRUE)  #End of the try block
    if(inherits(log, "try-error")){outList[[paste0('y',yrs[y])]] <- yearList} 
  }
  
  return(outList)
}






