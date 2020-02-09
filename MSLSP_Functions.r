#Functions for MuSLI_LSP - Multisource Land Imaging of Land Surface Phenology

#Functions written by Josh Gray, Douglas Bolton, and Eli Melaas
###############



#---------------------------------------------------------------------
#Get Default algorithm parameters for MuSLI LSP
#Douglas Bolton and Josh Gray
#---------------------------------------------------------------------
DefaultPhenoParameters_HLS <- function(){

  pheno_pars <- list(
    
    ## valid segment definition ##
    min_seg_amplitude=0.1,    # absolute amplitude threshold
    min_increase_length=30,    #Minimum increase length set to a month   (removes spurious short cycles)
    max_increase_length=185,   #Maximum increase length set to half a year (removes spurious long cycles)
    min_decrease_length=30,    #Minimum decrease length set to a month   (removes spurious short cycles)
    max_decrease_length=185,   #Maximum decrease length set to half a year (removes spurious long cycles)
    rel_amp_frac=0.35, # segment amplitude must be >= rel_amp_frac * global_amp (max - min)
    rel_peak_frac=NA, # EVI value at peak must be >= rel_peak_frac * global max
    
    ## date thresholds ##
    gup_threshes=c(0.15,0.5,0.9),
    gdown_threshes=c(0.9,0.5,0.15),
    
    ## image and smoothing parameters ##        
    dormantQuantile=0.05,   #Quantile for the dormant value
    dormStart=as.Date('2016-01-01'), #Start date for getting the dormant value
    dormEnd=as.Date('2018-12-31'),   #End date for getting the dormant value
   
    spikeThresh=2,             #A spike has to deviate X times the difference between the first and last values
    minResid=0.1,             #A spike has to deviate at least X EVI2 in order to be removed
    maxDistance=45,            #Three observations must be within X days in order to assess a spike 
    maxDespikeIterations=50,   #Maximum iterations to look for spikes
    
    splineBuffer=185,                #Spline buffer in days (~6 months)
    splineSpar=0.55,                 #Spline smoothing parameter
    maxWeight=0.1,                   #Maximum spline weight for additional years
    maxGapForSplineComparison=45,    #If data gap is greater than X, spline fit will be poor and shouldn't be used for comparison against other years (set to 6 weeks)
    minDaysForSplineComparison=90,   #In order to compare splines from different years, the splines must overlap by at least X days 
    howToCompare='euclidean',        #euclidean or correlation
    
    
    ndmiSnowThresh=0.5,   
    distanceToSnow=5000,    #If snow is detected within this distance (in meters), we will perform ndmi test
    snowWindow=5,           #Number of pixels to filter snow before creating distance to snow search window
    snowFraction=0.5,       #Fraction of pixels in snow window that must be snow in order to create search window
    snowBuffer=100,         #Distance to buffer detected snow pixels (in meters). These pixels will no be labelled as snow, but will be masked
    
    snowWeight=0.50,        #What weight should we apply to snow observations in the spline?
    snowFillVal=32767,      #Value to fill in for snow pixels

    ## topo correction parameters specific to VI approach ##
    numILclass=5,
    numSamples=1000,
    topoVIs=c('ndvi','nbr'),       #What indices to use to stratify the landscape
    viQuantiles=0.9,          #What temporal quantiles to take for each VI (calculated annually)
    requiredDoyStart=60,       #If the images for the first year don't go back until at least DOY 60, then use the next year's VIs
    requiredDoyEnd=300,        #If the images for the most recent year don't go until at least DOY 300, then use the previous year's VIs
    kmeansClasses=5,
    kmeansIterations=500,
    
    r2_qa=0.75,            #Cuts offs for quality classes. If r2 of spline fit is less than 0.75, we are considering this a poor fit
    maxGap_qa=c(20,35),    #Cuts offs for quality classes (high quality, moderate quality). Maximum image gap as a percentage of segment length

    numLyrs=107,          #How many output layers are there?
    runFmaskS30=TRUE,   #Should we rerun Fmask on sentinel images?
    runOn='SCC'         #SCC or AWS. Needed for setting environmental variables for matlab
  )
  
  return(pheno_pars)
}




#---------------------------------------------------------------------
#Sort out image chunk boundaries
#Douglas Bolton
#---------------------------------------------------------------------
chunkBoundaries <-function(numPix, numChunks) {

lenChunk <- ceiling(numPix / numChunks) #Length of each chunk

chunkStart <- matrix(0,numChunks,1)
chunkEnd <- matrix(0,numChunks,1)
for (n in 1:numChunks) {
  chunkStart[n] <- lenChunk*(n-1) + 1
  chunkEnd[n] <- lenChunk*n}

chunkEnd[numChunks] <- numPix
boundaries <- list(chunkStart,chunkEnd)
return(boundaries)
}



#---------------------------------------------------------------------
#Derive image mask and snow mask from QA flags
#image mask includes cloud, cloud shadow, adjacent cloud, and snow pixels
#snow mask includes snow pixels only. 
#Douglas Bolton
#---------------------------------------------------------------------
getMasks <-function(qaName) {
  
  #Read in QA
  qaM <- readGDAL(qaName,silent=T)$band1

  #Make empty masks
  #Must have two separate masks, because pixels can have multiple labels (e.g., snow in fmask, cloud in lasrc)
  mask <- as.integer(matrix(0,length(qaM),1))
  snow <- as.integer(matrix(0,length(qaM),1))

  #Out of the image value for S30 is NA. Set mask to 1
  mask[is.na(qaM)] <- 1
  qaM[is.na(qaM)] <- 0
  
  #Out of the image value for L30 is 255. Set mask to 1
  mask[qaM == 255] <- 1
  qaM[qaM == 255] <- 0
  
  #If bits 6-7 are 11, 01 or 10 (high/average/low aerosol), don't label, just subtract values  (we aren't confident enough in these values yet)
  qaM[qaM >= 192] = qaM[qaM >= 192] - 192 
  qaM[qaM >= 128] = qaM[qaM >= 128] - 128 
  qaM[qaM >= 64] = qaM[qaM >= 64] - 64
  
  #check if water, but don't label mask! Could be a false water detection, and we don't want to remove those. We'll use ancillary data for water
  qaM[qaM >= 32] <- qaM[qaM >= 32] - 32
  
  #check if snow/ice, label snow layer as 1
  snow[(qaM >= 16)] <- 1
  qaM[qaM >= 16] <- qaM[qaM >= 16] - 16
  
  #if cloud shadow, adjacent to cloud, cloud, or cirrus, label as 1
  mask[(qaM >= 1)] <- 1
  
  mask <- mask == 1
  snow <- snow == 1
  
  out <- list(mask,snow)
  return(out)
}



#---------------------------------------------------------------------
#Apply QA mask to images
#TEMPORARY FIX where we are using Fmask 4.0 for Sentinel
#Chunk images and save to a temporary location as .Rds files
#If a chunk is empty (no data), don't write to disk. Huge I/O savings.
#Douglas Bolton
#---------------------------------------------------------------------
ApplyMask_QA_and_Fmask <-function(imgName, tile, waterMask, tempFold, inDir, fmaskFunction, chunkStart, chunkEnd, topoMethod, pheno_pars, deleteInput=F) {
  
  #Get sensor, names sorted
  imgName_strip = tail(unlist(strsplit(imgName,'/')),n = 1)
  sensor = unlist(strsplit(imgName_strip,'.',fixed = T))[2]
  tileTxt = unlist(strsplit(imgName_strip,'.',fixed = T))[3]
  date = unlist(strsplit(imgName_strip,'.',fixed = T))[4]
  outBase = paste0('HLS_',sensor,'_',tileTxt,'_',date)
  
  #Open the image bands
  ##############################################################################
  theBase <- paste0('HDF4_EOS:EOS_GRID:',imgName,':Grid:')
  if (sensor == 'L30') {
    blueName = paste0(theBase, 'band02'); greenName = paste0(theBase, 'band03')
    redName = paste0(theBase, 'band04'); nirName = paste0(theBase, 'band05')
    swirName = paste0(theBase, 'band06'); swir2Name = paste0(theBase, 'band07')
  } else {
    blueName = paste0(theBase, 'B02'); greenName = paste0(theBase, 'B03')
    redName = paste0(theBase, 'B04'); nirName = paste0(theBase, 'B8A')
    swirName = paste0(theBase, 'B11'); swir2Name = paste0(theBase, 'B12')} 
  
  blueVals = readGDAL(blueName,silent=T)$band1
  greenVals = readGDAL(greenName,silent=T)$band1
  redVals = readGDAL(redName,silent=T)$band1
  nirVals = readGDAL(nirName,silent=T)$band1
  swirVals = readGDAL(swirName,silent=T)$band1
  swir2Vals = readGDAL(swir2Name,silent=T)$band1
  
  
  #Mask for clouds and snow. Use QA flags for Landsat. Impliment Fmask 4.0 for Sentinel
  ################################################################################
  if (sensor == 'L30' | pheno_pars$runFmaskS30 == FALSE) {
    qaName = paste0('HDF4_EOS:EOS_GRID:',imgName, ':Grid:QA')
    maskList <- getMasks(qaName)
    mask <- maskList[[1]]
    snow <- maskList[[2]]
    snow[waterMask==1] <- 0   #Remove snow flags that are over water
  }
  
  if (sensor == 'S30' & pheno_pars$runFmaskS30 == TRUE) { 
    fmaskFold <- paste0(tempFold,'fmask/')
    #Execute system command to call fmask function
    RunFmask(imgName, tile, inDir, fmaskFold, fmaskFold, fmaskFunction,pheno_pars) 
    #Now read in that fmask layer
    
    fmaskName <- paste0(fmaskFold, gsub('hdf','',imgName_strip), "Fmask.tif")
    fmask <- readGDAL(fmaskName,silent=T)$band1     #If Fmask doesn't exist (didn't run), function will fail and error will be logged
    mask <- (fmask == 4) | (fmask == 2) | (fmask == 255)  #If cloud, cloud shadow, or no data
    snow  <- fmask == 3
    snow[waterMask==1] <- 0  #Remove snow flags that are over water (shouldn't be any though)
  }
  
  #If we aren't topo correcting, then we need to do the NDMI check now. If we are topo correcting, we will do this check after topo correction
  #If NDMI > 0.5 AND the pixel is within 5 km of detected snow, then mask
  if (topoMethod == 'None') {
    if (sum(snow) > 0) {
      ndmi <- calcIndex(nir=nirVals/10000,swirOne=swirVals/10000,whatIndex='ndmi')
      ndmiCheck = ndmi > pheno_pars$ndmiSnowThresh   #Find pixels that meet threshold (potential snow pixels)
  
      #We will only keep snow pixels when 50% of 5x5 pixel window is also snow (to remove speckle)
      imBlur  <- boxblur(as.cimg(matrix(snow,length(snow)^.5,length(snow)^.5)),pheno_pars$snowWindow,neumann = FALSE)    #neumann = FALSE assumes snow values of zero outside image edges
      snowCheck <- imBlur > pheno_pars$snowFraction & snow == 1
      
      #Calcuate the distance from each pixel to the nearest pixel detected as snow by fmask or lasrc
      sD <- matrix(as.matrix(distance_transform(snowCheck,1)),length(snow),1)
      sCheck <- (sD*30) < pheno_pars$distanceToSnow  #Find pixels where distance to snow is less than threshold (5 km)
      
      #Find pixels meeting ndmi rule that are within 5km of detected snow pixels (or pixels already detected as snow)
      toMask <- (ndmiCheck & sCheck) | snow == 1 
      
      #Buffer snow pixels by 100m
      mD <- matrix(as.matrix(distance_transform(as.cimg(matrix(toMask,length(snow)^.5,length(snow)^.5)),1)),length(snow),1)
      buffer <- (mD*30) < pheno_pars$snowBuffer
      
      #Mask these pixels. BUT we won't label as snow, as our confidence is lower, and we don't want to falsely fill dormant values. 
      mask[buffer]  <- TRUE
      
      #rast <- raster(blueName)
      #snowRast <- setValues(rast, toMask)
      #writeRaster(snowRast,filename='~/testSnow2.tif',format='GTiff',overwrite=TRUE)
    }   
  }
  
  
  bands <- cbind(blueVals,greenVals,redVals,nirVals,swirVals,swir2Vals)
  bands[bands < 0] <- NA    #Remove negative reflectance values
  bands[mask==1,]  <-  NA
  bands[waterMask==1,]  <-  NA
  
  #Mask pixels that have missing data in ANY band. These are likely problem pixels (and we need most bands for despiking, kmeans, snow detection, etc).
  check <- rowSums(is.na(bands)) > 0
  bands[check,] <- NA
  
  bands <- round(bands * 10000)    #multiply by 10000 so we can store as integer
  
  #We remove potentially snowy observations using NDMI, but we are only going to label pixels that were actually detected as snow as snow
  bands[snow==1,]  <-  pheno_pars$snowFillVal    

  for (n in 1:length(chunkStart)) {
    mat <- bands[chunkStart[n]:chunkEnd[n],]
    if (all(is.na(mat))) {next} #if there is no good data in the chunk, move to next chunk
    saveRDS(as.integer(mat),paste0(tempFold,'c',n,'/',outBase,'.Rds'))
  }
  #If requested, delete the input file once it is processed
  if (deleteInput == TRUE) {file.remove(imgName)}
}




RunFmask <-function(imgName, tile, auxFold, fmaskFold, outDir, fmaskFunction, pheno_pars) {
  
  #Different approach depending on if we are on SCC or AWS
  #If on SCC, need to set environmental variables 
  if (pheno_pars$runOn == 'SCC') {
    system(paste('MCRROOT="/share/pkg/matlab/2018b/install"',
                 'LD_LIBRARY_PATH=.:${MCRROOT}/runtime/glnxa64',
                 'LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64',
                 'LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64',
                 'LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/opengl/lib/glnxa64',
                 'export LD_LIBRARY_PATH',
                 paste('eval',fmaskFunction,tile,imgName,auxFold,fmaskFold),sep='\n'))
  
   } else {


    #Get sensor, names sorted
    imgName_strip <- tail(unlist(strsplit(imgName,'/')),n = 1)
    fileBase <- paste0(outDir, gsub('hdf','',imgName_strip))
    
    #Open the image bands
    ##############################################################################
    inBase <- paste0('HDF4_EOS:EOS_GRID:',imgName,':Grid:')
  
    bandNames <- c('B02','B03','B04','B8A','B11','B12','B10','B07','B08')
    for (b in bandNames) { 
      inName <- paste0(inBase, b)
      outName <- paste0(fileBase,b,'.tif')
      system(paste("gdal_translate -q",inName,outName))}
  
    
    #Get sun zenith and azimuth.
    #Some L8 scenes will have two values (two L8 images as inputs). We will take the average of these two values
    info = gdalinfo(imgName)
    ulx <- as.numeric(unlist(strsplit(info[pmatch('  ULX',info)],'='))[2]) + 15    #Adjust for matlab (needs center of pixel)
    uly <- as.numeric(unlist(strsplit(info[pmatch('  ULY',info)],'='))[2]) - 15    #Adjust for matlab (needs center of pixel)
    sza <- as.numeric(unlist(strsplit(info[pmatch('  MEAN_SUN_ZENITH_ANGLE(B01)',info)],'='))[2])
    saa <- as.numeric(unlist(strsplit(info[pmatch('  MEAN_SUN_AZIMUTH_ANGLE(B01)',info)],'='))[2])
    vza <- as.numeric(unlist(strsplit(info[pmatch('  MEAN_VIEW_ZENITH_ANGLE(B01)',info)],'='))[2])
    vaa <- as.numeric(unlist(strsplit(info[pmatch('  MEAN_VIEW_AZIMUTH_ANGLE(B01)',info)],'='))[2])  
  
    zone <- unlist(strsplit(info[pmatch('  HORIZONTAL_CS_NAME',info)],' '))
    zone <- zone[length(zone)]
    zone <- as.numeric(substr(zone,1,2))
  
    system(paste('eval',fmaskFunction,tile,fileBase,auxFold,fmaskFold,zone,ulx,uly,sza,saa,vza,vaa)) 
    
    #Delete the temporary images
    for (b in bandNames) {file.remove(paste0(fileBase,b,'.tif'))}
  }
}



#---------------------------------------------------------------------
#Calculate annual quantiles for spectral indices
#We will use these quantiles to create kmeans classes for topographic correction
#Douglas Bolton
#---------------------------------------------------------------------
getIndexQuantile <- function(chunk, numPix, tempFold, yr, errorLog, pheno_pars) {
  
  #Make empty output
  indexOut <- matrix(NA,numPix,length(pheno_pars$topoVIs)*length(pheno_pars$viQuantiles))
  
  #Get all images to process
  ######################
  chunkFold <- paste0(tempFold,'c',chunk,'/') 
  imgList <- list.files(path=chunkFold, pattern=glob2rx("HLS_*.Rds"), full.names=F)
  numImgs <- length(imgList)
  
  if (numImgs > 0) {
    yrdoy <- as.numeric(matrix(NA,numImgs,1))
    for (i in 1:numImgs) {
      imName <- gsub('.Rds','',imgList[i])
      yrdoy[i] <- as.numeric(unlist(strsplit(imName,'_',fixed = T))[4])}
    years <- as.numeric(format(as.Date(strptime(yrdoy, format="%Y%j")),'%Y'))
  
    imgList <- imgList[years == yr]   #Restrict to images in target year
    numImgs <- length(imgList)
    
    #Only run if more than one image in target year
    if (numImgs > 1) {
      b1 = matrix(NA,numPix,numImgs); b2 = matrix(NA,numPix,numImgs)
      b3 = matrix(NA,numPix,numImgs); b4 = matrix(NA,numPix,numImgs)
      b5 = matrix(NA,numPix,numImgs); b6 = matrix(NA,numPix,numImgs)
      for (i in 1:numImgs) {
        imgData <- try(matrix(readRDS(paste0(chunkFold,imgList[i])),numPix,6),silent = TRUE)
        if (inherits(imgData, 'try-error')) {cat(paste('getIndexQuantile: Error for chunk',chunk,imgList[i]), file=errorLog, append=T);next} 
        imgData[imgData == pheno_pars$snowFillVal] <- NA
        imgData <- imgData / 10000
        b1[,i] <- imgData[,1]; b2[,i] <- imgData[,2]; b3[,i] <- imgData[,3]
        b4[,i] <- imgData[,4]; b5[,i] <- imgData[,5]; b6[,i] <- imgData[,6]
      }
      
      #Loop through requested indices in pheno_pars$topoVIs
      #For each index, loop through requested quantiles in pheno_pars$viQuantiles
      column <- 0
      for (i in 1:length(pheno_pars$topoVIs)) { 
        ind <- calcIndex(blue=b1,green=b2,red=b3,nir=b4,swirOne=b5,swirTwo=b6,whatIndex=pheno_pars$topoVIs[i])
        ind[is.infinite(ind)] <- NA
        for (j in 1:length(pheno_pars$viQuantiles)) {
          column <- column + 1
          indexOut[,column] <- rowQuantileC(ind,pheno_pars$viQuantiles[j])
        }
      }
    }
  }
  return(indexOut)
}





#---------------------------------------------------------------------
#Run topographic correction
#Overwrite existing image chunks with new, corrected data
#Douglas Bolton
#---------------------------------------------------------------------
runTopoCorrection <- function(imgName, tempFold, groups, slopeVals, aspectVals, chunkStart,chunkEnd, errorLog, pheno_pars, writeImages=FALSE, slope=FALSE) {
  
  #Get sensor, names sorted
  imgName_strip = tail(unlist(strsplit(imgName,'/')),n = 1)
  sensor = unlist(strsplit(imgName_strip,'.',fixed = T))[2]
  tileTxt = unlist(strsplit(imgName_strip,'.',fixed = T))[3]
  date = unlist(strsplit(imgName_strip,'.',fixed = T))[4]
  outBase = paste0('HLS_',sensor,'_',tileTxt,'_',date)

  #Get sun zenith and azimuth.
  #Some L8 scenes will have two values (two L8 images as inputs). We will take the average of these two values
  info = gdalinfo(imgName)
  line <- info[pmatch('  MEAN_SUN_ZENITH_ANGLE',info)]
  line <- unlist(strsplit(line,'='))[2]
  zMean <- mean(as.numeric(unlist(strsplit(line,','))))
  sunzenith <- (pi/180) * zMean
  
  line = info[pmatch('  MEAN_SUN_AZIMUTH_ANGLE',info)]
  line <- unlist(strsplit(line,'='))[2]
  aMean <- mean(as.numeric(unlist(strsplit(line,','))))
  sunazimuth <- (pi/180) * aMean
  
  numChunks <- length(chunkEnd)
  bands <- matrix(NA,chunkEnd[numChunks],6)

  #Reconstruct each image from the chunked data
  for (n in 1:numChunks) {
    fName <- paste0(tempFold,'c',n,'/',outBase,'.Rds')
    if (file.exists(fName)) {
      numPix <- length(chunkStart[n]:chunkEnd[n])
      imgData <- try(matrix(readRDS(fName),numPix,6),silent = TRUE)
      if (inherits(imgData, 'try-error')) {cat(paste('runTopoCorrection: Error for chunk',n,outBase), file=errorLog, append=T);next} 
      bands[chunkStart[n]:chunkEnd[n],] <- imgData
      remove(imgData)
    }
  }
  
  
  if (!all(is.na(bands))) {                   #Only continue if there is actually good data
    
    snow <- bands[,1] == pheno_pars$snowFillVal     #Determine snow pixels
    snow[is.na(snow)] <- FALSE                      #Set NA snow values to FALSE (no snow)
    bands[bands == pheno_pars$snowFillVal] = NA     #Mask snow values
    bands  <- bands / 10000                         
    
    #If image outputs are requested, create a base name for plots
    plotBaseName <- NULL
    if (writeImages) {plotBaseName <- paste0(tempFold,'IL_plots/',outBase)}   #Base name for outputing scatterplots

    #Run topographic correction function by kmeans group
    corr <- topocorr_rotational_by_group(x=bands,groups=groups, 
                                           slope=slopeVals,aspect=aspectVals, 
                                           sunzenith = sunzenith, sunazimuth=sunazimuth,pheno_pars,plotBaseName)

    #NDMI check
    #Now that we have topo corrected the imagery, let's check for potentail snow....
    #If NDMI > 0.5 AND the pixel is within 5 km of detected snow, then mask pixel
    if (sum(snow) > 0) {
      ndmi <- calcIndex(nir=corr[,4],swirOne=corr[,5],whatIndex='ndmi')
      ndmiCheck = ndmi > pheno_pars$ndmiSnowThresh   #Find pixels that meet threshold (potential snow pixels)
      
      #We will only keep snow pixels when 50% of 5x5 pixel window is also snow (to remove speckle)
      imBlur  <- boxblur(as.cimg(matrix(snow,length(snow)^.5,length(snow)^.5)),pheno_pars$snowWindow,neumann = FALSE)    #neumann = FALSE assumes snow values of zero outside image edges
      snowCheck <- imBlur > pheno_pars$snowFraction & snow == 1
      
      #Calcuate the distance from each pixel to the nearest pixel detected as snow by fmask or lasrc
      sD <- matrix(as.matrix(distance_transform(snowCheck,1)),length(snow),1)
      sCheck <- (sD*30) < pheno_pars$distanceToSnow  #Find pixels where distance to snow is less than threshold (5 km)
      
      #Find pixels meeting ndmi rule that are within 5km of detected snow pixels (or pixels already detected as snow)
      toMask <- (ndmiCheck & sCheck) | snow == 1 
      
      #Buffer snow pixels by 100m
      mD <- matrix(as.matrix(distance_transform(as.cimg(matrix(toMask,length(snow)^.5,length(snow)^.5)),1)),length(snow),1)
      buffer <- (mD*30) < pheno_pars$snowBuffer
      
      #Mask these pixels
      corr[buffer,]  <- NA
    } 
    
    corr = round(corr * 10000)                  #Convert back to integer
    corr[corr < 0]  <-  NA                      #Remove negative reflectance values
    
    #Mask pixels that have missing data in ANY band. These are likely problem pixels (and we need most bands for despiking, kmeans, snow detection, etc).
    check <- rowSums(is.na(corr)) > 0
    corr[check,] <- NA
    
    #We removed potentially snowy observations using NDMI, but we are only going to label pixels that were actually detected as snow (this impacts dormant filling)
    corr[snow==1,]  <-  pheno_pars$snowFillVal   
    corr <- round(corr)
    
    #Overwrite existing chunks with new corrected data
    for (n in 1:numChunks) {
      fName <- paste0(tempFold,'c',n,'/',outBase,'.Rds')
      if (file.exists(fName)) {file.remove(fName)}  #Delete the existing file. We will write the topo corrected file to the same location
      mat <- corr[chunkStart[n]:chunkEnd[n],]
      if (all(is.na(mat))) {next} #if there is no good data in the chunk, move to next chunk
      saveRDS(as.integer(mat),fName)
    }
    
    ####Should we write out some images for viewing?
    if (writeImages) {
      dType <- 'INT2U'
      ##Write out a corrected version. Put snow back to NA
      corr[corr == pheno_pars$snowFillVal] = NA
      after <- stack(setValues(slope,corr[,5]),setValues(slope,corr[,4]),setValues(slope,corr[,3]))   #SWIR/NIR/RED composite
      afterName <- paste0(tempFold,'/sampleImages/',outBase,'_corrected.tif')
      writeRaster(after,filename=afterName,format='GTiff',overwrite=TRUE,datatype=dType)
      
      ##Write out an uncorrected version
      bands = round(bands * 10000)
      bands <- round(bands)
      before <- stack(setValues(slope,bands[,5]),setValues(slope,bands[,4]),setValues(slope,bands[,3])) #SWIR/NIR/RED composite
      beforeName <- paste0(tempFold,'/sampleImages/',outBase,'.tif')
      writeRaster(before,filename=beforeName,format='GTiff',overwrite=TRUE,datatype=dType) 
    }
  }   
}



#---------------------------------------------------------------------
#Run Phenology code for an image chunk
#Write phenology results for each chunk to disk
#Douglas Bolton
#---------------------------------------------------------------------
runPhenoChunk <- function(chunk, tempFold, startYear, endYear, errorLog, pheno_pars, dataFill) {

    #Get all images to process
    ######################
    chunkFold <- paste0(tempFold,'c',chunk,'/') 
    outFold <-   paste0(tempFold,'outputs/')  
    
    imgList <- list.files(path=chunkFold, pattern=glob2rx("HLS_*.Rds"), full.names=F)
    
    numImgs = length(imgList)
    numPix = length(readRDS(paste0(chunkFold,imgList[1]))) / 6
    
    yrdoy = as.numeric(matrix(NA,numImgs,1))
    for (i in 1:numImgs) {
      imName <- gsub('.Rds','',imgList[i])
      yrdoy[i]= as.numeric(unlist(strsplit(imName,'_',fixed = T))[4])}
    
    ord = order(yrdoy)  #Determine image order
    
    #Read in all imagery for chunk
    b1 = matrix(NA,numPix,numImgs); b2 = matrix(NA,numPix,numImgs); b3 = matrix(NA,numPix,numImgs)
    b4 = matrix(NA,numPix,numImgs); b5 = matrix(NA,numPix,numImgs); b6 = matrix(NA,numPix,numImgs)
    sensor = as.character(rep('',numImgs))
    
    for (i in 1:length(ord)) {
      img <- imgList[ord[i]]
      imgData <- try(matrix(readRDS(paste0(chunkFold,img)),numPix,6),silent = TRUE)
      if (inherits(imgData, 'try-error')) {cat(paste('runPhenoChunk: Error for chunk',chunk,img), file=errorLog, append=T);next} 
      
      b1[,i] <- imgData[,1]; b2[,i] <- imgData[,2]; b3[,i] <- imgData[,3]
      b4[,i] <- imgData[,4]; b5[,i] <- imgData[,5]; b6[,i] <- imgData[,6]
      sensor[i] = unlist(strsplit(img,'_',fixed = T))[2]
      remove(imgData)
    }
    
    dates = as.Date(strptime(yrdoy[ord], format="%Y%j"))  #Format as dates
    
    #Get snow pixels to fill 
    snowPix <- b1 == pheno_pars$snowFillVal
    snowPix[is.na(snowPix)] <- FALSE
    
    #Set snow pixels to NA
    b1[snowPix] <- NA; b2[snowPix] <- NA; b3[snowPix] <- NA
    b4[snowPix] <- NA; b5[snowPix] <- NA; b6[snowPix] <- NA
    
    
    
    
    #Temporary mask of late date 
    ###########################################################
    # test <- dates < as.Date('2019-04-15')
    # dates <- dates[test]
    # snowPix <- snowPix[,test]
    # b1 <- b1[,test]
    # b2 <- b2[,test]
    # b3 <- b3[,test]
    # b4 <- b4[,test]
    # b5 <- b5[,test]
    # b6 <- b6[,test]
    # sensor <- sensor[test]
    
    
    
    
    #Calculate EVI2, mask negative values
    evi2 <- calcIndex(red=b3/10000,nir=b4/10000,whatIndex='evi2')
    evi2[evi2 <= 0] <- NA   
    
    #Average images that occur on the same day
    #Only keep the snow flag if both images on the day were classified as snow
    imgData <- averageDuplicates(list(b1, b2, b3, b4, b5, b6, evi2), snowPix, sensor, dates)
    
    #Extract data from listd
    bands <- imgData[[1]]; snowPix <- imgData[[2]]; sensor <- imgData[[3]]; dates <- imgData[[4]] 
    b1 <- bands[[1]];  b2 <- bands[[2]];  b3 <- bands[[3]]; 
    b4 <- bands[[4]]; b5 <- bands[[5]]; b6 <- bands[[6]]
    evi2 <- bands[[7]]
    remove(imgData)
    remove(bands)
    
    
    yrs <- seq(startYear,endYear)
    splineStart <- as.Date(as.Date(paste0(yrs,'-01-01')) - pheno_pars$splineBuffer) 
    numDaysFit  <-  365 + (pheno_pars$splineBuffer * 2)
    pheno_pars$halfLyrs <- (pheno_pars$numLyrs-1)/2    #Number of layers that require filling if there is no second cycle. First cycle has 7 extra layers
    
    
    #Repeat data for the buffer period of the final year to ensure that the spline gets tied down
    #Determine the most recent image, and then fill BEYOND that date until the spline buffer is reached
    #Filling with images from the previous year
    #########################################################################
    startFillDate <- max(dates) - 365    #Find the most recent image, and then subtract a year
    endFillDate <- as.Date(paste0(endYear,'-01-01')) + pheno_pars$splineBuffer   #Find the upper bound that needs to be filled
    
    vals <- dates > startFillDate & dates < endFillDate  
    b1_add <- b1[,vals];b2_add <- b2[,vals];b3_add <- b3[,vals]
    b4_add <- b4[,vals];b5_add <- b5[,vals];b6_add <- b6[,vals]
    evi2_add <- evi2[,vals]
    snow_add <- snowPix[,vals]
    dates_add <- dates[vals] + 365
    
    b1 <- cbind(b1,b1_add);b2 <- cbind(b2,b2_add);b3 <- cbind(b3,b3_add);
    b4 <- cbind(b4,b4_add);b5 <- cbind(b5,b5_add);b6 <- cbind(b6,b6_add);
    evi2 <- cbind(evi2,evi2_add);
    snowPix <- cbind(snowPix,snow_add);
    dates <- c(dates,dates_add)
    
    
    #Loop through each pixel and estimate phenometrics
    #Run DoPhenologyHLS_Climatology if we want to fill gaps with other years. Otherwise run DoPhenologyHLS.
    pheno_mat <- matrix(NA,numPix,pheno_pars$numLyrs*length(yrs))
    if (dataFill) {
        for (i in 1:numPix) {pheno_mat[i,] <- DoPhenologyHLS_Climatology(b1[i,],  b2[i,],  b3[i,],  b4[i,],  b5[i,],  b6[i,],  evi2[i,],
                                                  snowPix[i,],dates, yrs, splineStart, numDaysFit, pheno_pars)}
    } else { 
      for (i in 1:numPix) {pheno_mat[i,] <- DoPhenologyHLS(b1[i,],  b2[i,],  b3[i,],  b4[i,],  b5[i,],  b6[i,],  evi2[i,],
                                        snowPix[i,],dates, yrs, splineStart, numDaysFit, pheno_pars)}
    }
    
    pheno_mat <- round(pheno_mat)
    
    #Write results to disk (.Rds files)
    for (y in 1:length(yrs)) {
      yr <- yrs[y]
      ind <- ((y-1)*pheno_pars$numLyrs+1):(y*pheno_pars$numLyrs)
      subTab <- pheno_mat[,ind]
      for (i in 1:pheno_pars$numLyrs) {saveRDS(as.integer(subTab[,i]),paste0(outFold,'y',yr,'/lyr',i,'/c',chunk,'.Rds'))}}
}



#---------------------------------------------------------------------
#Correct reflectance for topography
#Developed by Eli Melaas
#Rewritten by Douglas Bolton to only run rotational model (Tan et al. 2010) and to have a stratified sample of pixels based on IL
#---------------------------------------------------------------------
topocorr_rotational_by_group <-function(x, groups, slope, aspect, sunzenith, sunazimuth, pheno_pars, plotBaseName, IL.epsilon=0.000001) {
  # topographic correction for image x based on topography and sun location
  
  #IMPORTANT: slope, aspect, sunzenith, sunazimuth must be in radians! 
  
  #Calculate IL
  IL <- cos(slope) * cos(sunzenith) + sin(slope) * sin(sunzenith) * cos(sunazimuth - aspect)
  IL[IL <= 0] <- IL.epsilon # adding tiny increment eliminates Inf's later
  

  #Determine # of bands and create empty output matrix
  if (is.null(dim(x))) {numBands <- 1} else {numBands <- dim(x)[2]}
  
  xout <- x    #IMPORTANT: Set output to input. This means, if is no slope/aspect info, or too few pixels to correct, the ORIGINAL band values will be output instead of NA values
  #Remove NAs from groups
  groups[is.na(groups)] = 0
  
  #samples per class
  numPer <- round(pheno_pars$numSamples/pheno_pars$numILclass)
  
  #Loop through bands
  for (b in 1:numBands) {
    band <- x[,b]
    
    clSub <- groups
    clSub[is.na(band)] = 0
    
    clTypes <- unique(clSub)
    clTypes <- clTypes[clTypes != 0]
    
    #Loop through kmeans classes
    for (clType in clTypes) {
      clCheck <- clSub == clType
      
      if (sum(clCheck) < 5){next}   #If fewer than 5 pixels, don't correct. These pixels will be masked from the output
      
      bandSub <- band[clCheck]
      ilSub <- IL[clCheck]
      
      ##Break IL into groups. Sample within each group
      breaks <- seq(quantile(ilSub,prob=0.02),quantile(ilSub,prob=0.98),length.out=(pheno_pars$numILclass+1))    #Using 2 and 98% percentiles to reduce influece of outliers
      breaks[1] <- min(ilSub); breaks[length(breaks)] <- max(ilSub)
      ilGroup <- try({as.numeric(cut(ilSub,breaks=breaks,labels=1:pheno_pars$numILclass,include.lowest=T))},silent=T)  #Group pixels
      if (inherits(ilGroup, 'try-error')) {next}    #If pixels can't be grouped, move on
      pixID <- 1:length(ilGroup)
      pixToSample <- matrix(0,length(pixID))
      for (i in 1:pheno_pars$numILclass) {
        check <- ilGroup == i
        num <- sum(check)
        if (num == 0) {next}
        if (num >= numPer) {toSamp <- numPer} else {toSamp <- num}   
        pixSub <- sample(pixID[check],toSamp)    #Randomly sample pixels in each IL group
        pixToSample[pixSub] = 1
      }
      
      #Perform the correction
      band.lm <- lm(as.vector(bandSub[pixToSample == 1]) ~ as.vector(ilSub[pixToSample == 1]))
      A <- coefficients(band.lm)[[2]]
      xout[clCheck,b] <- bandSub - A * (ilSub - cos(sunzenith))    #Output the corrected band values
      
      #If image outputs were requested, then create some IL scatterplots
      if (!is.null(plotBaseName)) {
        try({
          if (b == 3 | b == 4) {   #only for RED and NIR channel
            
            output_name <- paste0(plotBaseName,'_b',b,'_c',clType,'.png')
            png(output_name,res = 300,width = 6,height = 3,units = "in",pointsize = 8)
            par(mai=c(.6,.6,.3,.1),fig=c(0.02,0.98,0.02,0.98))
            plot(as.vector(ilSub[pixToSample == 1]), as.vector(bandSub[pixToSample == 1]),pch=16,xlab='IL',ylab = paste('Band',b),main=paste('Class',clType))
            abline(band.lm)
            dev.off()
            
            #Output the data points to make figures later
            toWrite <- cbind(as.vector(ilSub[pixToSample == 1]), as.vector(bandSub[pixToSample == 1]))
            saveRDS(toWrite,paste0(plotBaseName,'_b',b,'_c',clType,'.Rds'))
          }
        })
      }
    }
  }
  return(xout)
}







#---------------------------------------------------------------------
#Generic function to calculate indicies
#Must provide individual bands and an index to calculate
#Douglas Bolton
#---------------------------------------------------------------------
calcIndex <- function(blue=NULL,green=NULL,red=NULL,edge1=NULL,edge2=NULL,edge3=NULL,nir=NULL,swirOne=NULL,swirTwo=NULL,whatIndex=NULL) {
  
  if (whatIndex == 'evi') {index <- 2.5*(nir - red) / (nir + 2.4*red - 7.5*blue)}
  if (whatIndex == 'evi2') {index <- 2.5*(nir - red) / (nir + 2.4*red + 1)}
  if (whatIndex == 'ndvi') {index <- (nir - red) / (nir + red)}
  if (whatIndex == 'ndmi') {index <- (nir - swirOne) / (swirOne + nir)}
  if (whatIndex == 'nbr') {index <- (nir - swirTwo) / (swirTwo + nir)}
  if (whatIndex == 'nbr2') {index <- (swirOne - swirTwo) / (swirOne + swirTwo)}
  if (whatIndex == 'savi') {index <- ((nir - red) / (nir + red + 0.5)) * 1.5}
  if (whatIndex == 'msavi') {index <- (2*nir + 1 - sqrt((2*nir + 1)^2 - 8*(nir - red))) / 2}
  if (whatIndex == 'rcc') {index <- red / (red + green + blue)}
  if (whatIndex == 'gcc') {index <- green / (red + green + blue)}
  if (whatIndex == 'grvi') {index <- (green-red) / (red + green)}
  if (whatIndex == 'ndsi') {index <- (green-swirOne) / (green+swirOne)}
  if (whatIndex == 'npci') {index <- (red - blue) / (red + blue)}     #Hatfield et al. 2010
  if (whatIndex == 'CIgreen') {index <- (nir/green) - 1}     #Hatfield et al. 2010
  if (whatIndex == 'psri') {index <- (red-green) / (nir)}    #Hatfield et al. 2010
  
  if (whatIndex == 'brightness') {index <- 0.3029*blue + 0.2786*green + 0.4733*red + 0.5599*nir + 0.508*swirOne + 0.1872*swirTwo}
  if (whatIndex == 'greenness') {index <- -0.2941*blue + -0.243*green + -0.5424*red + 0.7276*nir + 0.0713*swirOne + -0.1608*swirTwo}
  if (whatIndex == 'wetness') {index <- 0.1511*blue + 0.1973*green + 0.3283*red + 0.3407*nir + -0.7117*swirOne + -0.4559*swirTwo}
  
  #Red edge indices
  if (whatIndex == 'ndvi_re1') {index <- (nir-edge1) / (nir + edge1)}
  if (whatIndex == 'nredi1') {index <- (edge2-edge1) / (edge2 + edge1)}
  if (whatIndex == 'nredi2') {index <- (edge3-edge1) / (edge3 + edge1)}
  if (whatIndex == 'nredi3') {index <- (edge3-edge2) / (edge3 + edge2)}
  if (whatIndex == 'psri1') {index <- (red-green) / edge1} 
  
  #Just bands
  if (whatIndex == 'blue') {index <- blue}
  if (whatIndex == 'green') {index <- green}
  if (whatIndex == 'red') {index <- red}
  if (whatIndex == 'nir') {index <- nir}
  if (whatIndex == 'swirOne') {index <- swirOne} 
  if (whatIndex == 'swirTwo') {index <- swirTwo} 

  return(index)
}




#---------------------------------------------------------------------
#Find dates with > 1 image
#Average the images if they both have observations
#Take the good image if one is masked by clouds, snow, etc
#Also adjust the snow mask, sensor vector, and date vector to reflect the removal of duplicates
#Douglas Bolton
#---------------------------------------------------------------------
averageDuplicates <-function(bands, snowPix, sensor, dates) {
  
    dups <- duplicated(dates)
    dupDates <- dates[dups]
    
    if (length(dupDates) > 0) {
  
      numBands <- length(bands)
      
      #Average each band and index
      bandOut <- vector("list",numBands)
      for (b in 1:length(bands)) {
        band <- bands[[b]]
        for (d in 1:length(dupDates)) {
          dd <- dates == dupDates[d]
          band[,dd] = rowMeans(band[,dd],na.rm = T)}
        bandOut[[b]] <- band[,dups==0]                    #Remove duplicated dates
        }
  
      #Change sensor to Avg if there are two observations on the same day
      sensorOut <- sensor
      for (d in 1:length(dupDates)) {sensorOut[dates == dupDates[d]] <- 'Avg'} 
      sensorOut <- sensorOut[dups==0]  #Remove duplicated dates
      
      
      #Label as snow only if both images are labelled as snow
      noSnowPix <- snowPix == 0    #Change so 0 = snow, 1 = no snow
      for (d in 1:length(dupDates)) {
        dd <- dates == dupDates[d]
        noSnowPix[,dd] <- rowSums(noSnowPix[,dd],na.rm = T)}  #Check if any observation from the day was NO snow (sum will be greater than zero)
      snowPix <- noSnowPix == 0   #Label as snow if ALL images on the particular day are labelled as snow (noSnowPix == 0)
      snowOut <- snowPix[,dups==0] #Remove duplicated dates
      
      #Remove duplicated dates from date vector
      datesOut <- dates[dups==0]
      
      out <- list(bandOut,snowOut,sensorOut,datesOut)
      
    } else  {out <- list(bands, snowPix, sensor, dates)}
    
    return(out)
}




#---------------------------------------------------------------------
#Despike time-series
#Uses the Maja approach to get positive time-series spikes associated with clouds
#See Section 2.1.4 in this document: https://www.theia-land.fr/sites/default/files/imce/produits/atbd_maja_071217.pdf
#The maja approach relies on both the blue and red bands (If change is red is 1.5x the change in blue, it is likely a land condition change)
#In addition, a 3 point method is used to look for negative spikes in EVI2. This step is added because the Maja approach won't capture
#spikes associated with cloud shadow
#Vectorized to avoid looping through dates
#Douglas Bolton
#---------------------------------------------------------------------
CheckSpike_MultiBand <- function(blue,red,evi2, dates, pheno_pars){
  
  blue_og <- blue # preserve original vector
  red_og <- red # preserve original vector
  evi2_og <- evi2 # preserve original vector
  dates_og <- as.numeric(dates)
  
  good <- !is.na(blue_og) & !is.na(red_og) & !is.na(evi2_og)
  
  x_outs <- matrix(F, length(blue_og)) # create the outlier output vector
  count <- 0 
  while (count < pheno_pars$maxDespikeIterations) {
    count <- count + 1
    
    bS <- blue_og[good] # subset to non missing values
    rS <- red_og[good]
    eS <- evi2_og[good] 
    dS <- dates_og[good] #subset date vector 
    
    ind1 <- 1:(length(dS)-2)  #Get indices for first, second, and third images
    ind2 <- 2:(length(dS)-1)
    ind3 <- 3:length(dS)
    
    dDiff1 <- dS[ind2] - dS[ind1]
    bDiff1 <- bS[ind2] - bS[ind1]
    rDiff1 <- rS[ind2] - rS[ind1]
    bTest1 <- bDiff1 > (0.03 * (1 + dDiff1/30))
    rTest1 <- rDiff1 < (1.5 * bDiff1)
    
    dDiff2 <- dS[ind3] - dS[ind2]
    bDiff2 <- bS[ind2] - bS[ind3]    #2 minus 3 because we are investigating 2 as a peak
    rDiff2 <- rS[ind2] - rS[ind3]
    bTest2 <- bDiff2 > (0.03 * (1 + dDiff2/30))
    rTest2 <- rDiff2 < (1.5 * bDiff2)
    
    majaTest <- bTest1 & rTest1 & bTest2 & rTest2
    
    dayFrac <- (dS[ind2]-dS[ind1]) / (dS[ind3]-dS[ind1])   #Calculate time fraction of date 1 to 2 compared to date 1 to 3
    fitVal <- eS[ind1] + (eS[ind3] - eS[ind1]) * dayFrac   #Calculate value at point 2 if a straight line is drawn from point 1 to 3.  
    dev1 <- eS[ind2] - eS[ind1]
    dev2 <- eS[ind2] - eS[ind3]
    dev <- fitVal - eS[ind2]
    devRatio <- dev / (eS[ind3] - eS[ind1])
    dDiff <- dS[ind3] - dS[ind1]
    
    #look for negative spikes in evi2
    eTest <- (dev > pheno_pars$minResid) & (abs(devRatio) > pheno_pars$spikeThresh) & (dDiff < pheno_pars$maxDistance)   
    
    check <- majaTest | eTest
    check <- c(FALSE,check,FALSE)
    check[is.na(check) | is.infinite(check)] <- FALSE
    if (sum(check) == 0) {break}    #Break if no observations have been despiked
    x_outs[good] <- check   #expand to size of original x, accounting for missing values
    good[x_outs] <- FALSE              #remove the despiked values from the pixels of interest and try again 
  }
  return(x_outs)
}







#---------------------------------------------------------------------
#Check the snow flags to see if the dormant value would be a good fit
#Specifically, if filling with the dormant value creates a negative spike in the time-series, we will remove the snow flag
#Outputs an updated vector of snow flags
#Douglas Bolton
#---------------------------------------------------------------------
Screen_SnowFills <- function(x, dormVal, snowPix, dates){
  #Reduce vectors to valid obervations and snow pixels
  sORg<- !is.na(x) | snowPix==1
  snowSub <- snowPix[sORg]
  xSub <- x[sORg]
  d <- as.numeric(dates[sORg])
  
  padded <- c(0,rollsum(snowSub,3),0)   #Get rolling sum of snow detections
  padded2 <- c(3,rollmax(padded,3),3)   #Get rolling max of the sums. Setting first and last to 3 keeps these values if they are snow
  
  #If the rolling max is less than 3, means center value isn't part of three consecutive snow values
  #So we will only check values less than 3 to determine if they are spikes
  toCheck <- which(padded2 < 3 & snowSub == 1)  
  noNA <- which(!is.na(xSub))
  rmIt <- matrix(FALSE,length(snowSub))
  for (ind2 in toCheck) {
    inds1 <- noNA[noNA < ind2]
    inds3 <- noNA[noNA > ind2]
    if (length(inds1)==0 | length(inds3)==0) {rmIt[ind2] <- TRUE;next}   #If there are no observations 
    ind1 <- max(inds1)
    ind3 <- min(inds3)
    
    dayFrac <- (d[ind2]-d[ind1]) / (d[ind3]-d[ind1])   #Calculate time fraction of date 1 to 2 compared to date 1 to 3
    fitVal <- xSub[ind1] + (xSub[ind3] - xSub[ind1]) * dayFrac   #Calculate value at point 2 if a straight line is drawn from point 1 to 3.  
    dev <- fitVal - dormVal
    devRatio <- dev / (xSub[ind3] - xSub[ind1])
    dDiff <- d[ind3] - d[ind1]
    
    rmIt[ind2] <- (dev > pheno_pars$minResid) & (abs(devRatio) > pheno_pars$spikeThresh) 
  }
  
  rmFull <- matrix(FALSE,length(snowPix))
  rmFull[sORg] <- rmIt
  snowPix[rmFull] <- FALSE
  return(snowPix)
}













#----------------------------------------------------------
# Fit a cubic spline to the HLS time-series
# Written by Josh Gray and Douglas Bolton
#----------------------------------------------------------
Just_Smooth <- function(x, dates, pred_dates, weights, dormant_value, pheno_pars){
    #Get index of pixels with good values
    ind <- !is.na(x)  
    # smooth with a spline to get continuous daily series
    spl <- smooth.spline(dates[ind], x[ind], spar=pheno_pars$splineSpar, w=weights[ind])
    # weighted version
    xSmooth <- predict(spl, as.numeric(pred_dates))$y
    # screen and fill values less than the the dormant value
    xSmooth[xSmooth < dormant_value] <- dormant_value
    return(xSmooth)
}





#----------------------------------------------------------
#New from Josh - 2018-10-31
#Finds time-series peaks
#Josh Gray
#----------------------------------------------------------
FindPeaks <- function(x, mag_order=T){
  # Function to identify peaks in time series x (or troughs if x=-x), supports "flat top" peaks
  # if mag_order is TRUE, peaks are returned in order of increasing magnitude (of x)
  d <- diff(x)
  d_code <- (d > 0) + (2 * (d < 0)) # 0=no change, 1=inc, 2=dec
  peaks <- unlist(gregexpr("12", paste(d_code, collapse=""))) # no match is -1
  if(peaks[1] == -1) peaks <- NULL
  flat_peaks <- unlist(gregexpr("10+2", paste(d_code, collapse=""))) # no match is -1
  if(flat_peaks[1] == -1) flat_peaks <- NULL
  d_code_rle <- rle(d_code)
  flat_peaks <- flat_peaks + round(d_code_rle$l[match(flat_peaks, cumsum(d_code_rle$l)) + 1] / 2)
  # all_peaks <- c(ifelse(peaks[1] == -1, NULL, peaks + 1), ifelse(flat_peaks[1] == -1, NULL, flat_peaks + 1))
  peaks <- sort(c(peaks + 1, flat_peaks + 1))
  if(mag_order) return(peaks[order(x[peaks])])
  return(peaks)
}


#----------------------------------------------------------
#New from Josh - 2018-10-31
#Determines valid segments 
#Josh Gray
#----------------------------------------------------------
GetSegs <- function(peaks, x, pars, peak=NA){
  # identifies valid increasing-decreasing segments in x subject to the parameters in pars
  # returns a list of segments: c(start, peak, end). DON'T call directly w/ peak!=NA
  # NOTE: returned segments will not necessarily be in order, and may not completely partition x
  
  # ensure that peaks are in increasing order of x's magnitude
  tmp_peaks <- peaks[order(x[peaks])] # so we only have to sort once if they're in the wrong order
  if(!identical(tmp_peaks, peaks)) peaks <- tmp_peaks
  
  # if no peak is specified, we start at the beginnning
  if(is.na(peak)) peak <- peaks[1]
  
  # get the next largest peak; will be NA if this peak is the highest (last one to do)
  next_highest_peak <- peaks[which(peaks == peak) + 1]
  
  # check if we're doing C5-style relative amplitude and peak identification
  # if(!is.na(pars$rel_amp_frac) & !is.na(pars$rel_peak_frac)){
  #   global_max <- max(x, na.rm=T)
  #   seg_thresh <- (global_max - min(x, na.rm=T)) * pars$rel_amp_frac
  #   peak_thresh <- global_max * pars$rel_peak_frac
  # }else{
  #   seg_thresh <- pars$min_seg_amplitude
  #   peak_thresh <- 0
  # }
  
  # we could have any combinaton of rel_amp_frac, rel_peak_frac, and min_seg_amplitude specified
  # initialize seg_thresh and peak_thresh to zero
  # determine the "global max/min", if peak_frac is specified, set it, if amp_frac is specified, set it
  # if min_seg_amplitude is set, choose the max of that and amp_frac
  seg_thresh <- peak_thresh <- 0
  global_max <- max(x, na.rm=T)
  global_min <- min(x, na.rm=T)
  if(!is.na(pars$rel_amp_frac)) seg_thresh <- (global_max - global_min) * pars$rel_amp_frac
  #if(!is.na(pars$rel_peak_frac)) peak_thresh <- global_max * pars$rel_peak_frac
  if(!is.na(pars$min_seg_amplitude)) seg_thresh <- max(pars$min_seg_amplitude, seg_thresh)
  
  # checks if the period preceding the peak covers enough amplitude
  # search before the peak up to the maximum of: previous peak, the head of x, or the peak - max_increase_length
  previous_peaks <- peaks[peaks - peak < 0]
  previous_peak <- NA
  if(length(previous_peaks) > 0) previous_peak <- max(previous_peaks)
  search_start <- max(1, peak - pars$max_increase_length, previous_peak, na.rm=T)
  search_end <- peak
  # get the index of the closest minimum value within the search window
  # NOTE: should maybe retrieve the troughs here with FindPeaks(-x) instead
  # in the event of repeated minimum values, we take the closest one here
  inc_min_ind <- max(which(x[search_start:search_end] == min(x[search_start:search_end], na.rm=T)) + search_start - 1, na.rm=T)
  seg_amp <- x[peak] - x[inc_min_ind] # get the increasing segment amplitude
  # if(seg_amp > pars$min_seg_amplitude){
  if((seg_amp >= seg_thresh) & (x[peak] >= peak_thresh)){
    # check for a valid decreasing segment
    next_peaks <- peaks[peaks - peak > 0]
    next_peak <- NA
    if(length(next_peaks) > 0) next_peak <- min(next_peaks)
    # search after the peak up to the minimum of: next peak, the tail of x, or the max_decrease_length
    search_start <- peak
    search_end <- min(length(x), peak + pars$max_decrease_length, next_peak, na.rm=T)
    # get the index of the closest minimum value within the search window
    # NOTE: see above note about finding troughs instead
    dec_min_ind <- min(which(x[search_start:search_end] == min(x[search_start:search_end], na.rm=T)) + search_start - 1, na.rm=T)
    seg_amp <- x[peak] - x[dec_min_ind] # get the decreasing segment amplitude
    # if(seg_amp > pars$min_seg_amplitude){
    if(seg_amp >= seg_thresh){
      # we found a valid segment, store it as a list with a single vector: c(start, peak, end)
      tmp_seg <- list(c(inc_min_ind, peak, dec_min_ind))
      # if this isn't the last peak, then call CheckSegRec again w/ next highest peak
      if(!is.na(next_highest_peak)){
        return(c(tmp_seg, GetSegs(peaks, x, pars, peak=next_highest_peak)))
      }else{
        # that was the last peak, and it was valid
        return(tmp_seg) # covers the case where there's only one valid peak
      }
    }else{
      # increase was valid, but decrease was not
      peaks <- peaks[-which(peaks == peak)] # remove peak from peaks list
      # if this isn't the last peak, then call CheckSegRec again w/ next highest peak
      if(!is.na(next_highest_peak)){
        return(GetSegs(peaks, x, pars, peak=next_highest_peak))
      }else{
        # that was the last peak, and it was invalid
        return(NULL)
      }
    }
  }else{
    # increase segment not valid
    peaks <- peaks[-which(peaks == peak)] # remove peak from peaks list
    # if this isn't the last peak, then call CheckSegRec again w/ next highest peak
    if(!is.na(next_highest_peak)){
      return(GetSegs(peaks, x, pars, peak=next_highest_peak))
    }else{
      # that was the last peak, and it was invalid
      return(NULL)
    }
  }
}



#----------------------------------------------------------
#Get phenology dates from segments
#Josh Gray
#----------------------------------------------------------
GetPhenoDates <- function(segs, x, dates, pheno_pars){
  pheno_dates <- list()
  for(gup_thresh in pheno_pars$gup_threshes){
    tmp_dates <- list(dates[unlist(lapply(segs, GetSegThresh, x, gup_thresh, gup=T), use.names=F)])
    if(all(is.na(unlist(tmp_dates, use.names=F)))){
      tmp_dates <- NA
    }
    pheno_dates <- c(pheno_dates, tmp_dates)
  }
  
  for(gdown_thresh in pheno_pars$gdown_threshes){
    tmp_dates <- list(dates[unlist(lapply(segs, GetSegThresh, x, gdown_thresh, gup=F), use.names=F)])
    if(all(is.na(unlist(tmp_dates, use.names=F)))){
      tmp_dates <- NA
    }
    pheno_dates <- c(pheno_dates, tmp_dates)
  }
  return(pheno_dates)
}


#----------------------------------------------------------
#Josh Gray
#----------------------------------------------------------
GetThresh <- function(thresh_value, x, first_greater=T, gup=T){
  # returns the index of the first/last value  of x that is greater/less than the value of thresh.
  # If gup is False (greendown) then it returns the first/last value of x that is less/greater than
  # the value of thresh. first/last and greater/less determined by first_greater
  # NOTE: if thresh is 1 or 0, rounding error can be a problem. Now we round the threshold and each
  # of the evi values to 6 decimal places to compensate
  
  if(gup){
    if(first_greater){
      return(min(which(round(x, 6) >= round(thresh_value, 6))))
    }else{
      return(max(which(round(x, 6) <= round(thresh_value, 6))))
    }
  }else{
    if(first_greater){
      return(min(which(round(x, 6) <= round(thresh_value, 6))))
    }else{
      return(max(which(round(x, 6) >= round(thresh_value, 6))))
    }
  }
}

#----------------------------------------------------------
#Josh Gray
#----------------------------------------------------------
GetSegThresh <- function(seg, x, thresh, gup=T){
  if(gup){
    # check for valid greenup segment
    if(!is.na(seg[1]) & !is.na(seg[2])){
      gup_thresh <- x[seg[1]] + ((x[seg[2]] - x[seg[1]]) * thresh)
      gup_thresh_index <- GetThresh(gup_thresh, x[seg[1]:seg[2]], first_greater=T, gup=T)
      return(gup_thresh_index + seg[1] - 1)
    }else{
      return(NA)
    }
  }else{
    # check for valid greendown segment
    if(!is.na(seg[2]) & !is.na(seg[3])){
      gdown_thresh <- x[seg[3]] + ((x[seg[2]] - x[seg[3]]) * thresh)
      gdown_thresh_index <- GetThresh(gdown_thresh, x[seg[2]:seg[3]], first_greater=T, gup=F)
      return(gdown_thresh_index + seg[2] - 1)
    }else{
      return(NA)
    }
  }
}



#----------------------------------------------------------
#Developed by Josh Gray, updated by Douglas Bolton for HLS
#----------------------------------------------------------
GetSegMetrics <- function(seg, x_smooth, x_raw, smooth_dates, raw_dates){
  if(any(is.na(seg))){return(NA)}
  # get the subset of the smoothed and original time series
  tmp_seg_smooth <- x_smooth[seg[1]:seg[3]]
  tmp_gup_smooth <- x_smooth[seg[1]:seg[2]]
  tmp_gdown_smooth <- x_smooth[seg[2]:seg[3]]
  
  # get the full segment minimum/maximum SVI
  seg_min <- min(tmp_seg_smooth, na.rm=T)
  seg_max <- max(tmp_seg_smooth, na.rm=T)
  
  
  # get the segment integrated SVI: the sum of values.
  #For MODIS C6, this is the sum of values above the minimum evi. 
  seg_int <- sum(tmp_seg_smooth)
  
  
  # organize greenup segment
  ######################################
  gup_raw_date_inds <- which(raw_dates >= smooth_dates[seg[1]] & raw_dates <= smooth_dates[seg[2]]) # indices in raw data of gup segment
  gup_smooth_date_inds <- match(raw_dates[gup_raw_date_inds], smooth_dates) # indices of raw dates in smooth dates
  
  raw_dates_gup <- raw_dates[gup_raw_date_inds]
  gup_raw_data <- x_raw[gup_raw_date_inds] # get the raw data associated with the gup segment (this is the pre-filled, despiked version)
  gup_smooth_data <- x_smooth[gup_smooth_date_inds] # get the smoothed values associated with each raw data value
  
  gup_numObs <- sum(!is.na(gup_raw_data)) 
  
  
  # organize greendown segment
  ######################################
  gdown_raw_date_inds <- which(raw_dates >= smooth_dates[seg[2]] & raw_dates <= smooth_dates[seg[3]]) # indices in raw data of gdown segment
  gdown_smooth_date_inds <- match(raw_dates[gdown_raw_date_inds], smooth_dates) # indices of raw dates in smooth dates
  
  raw_dates_gdown <- raw_dates[gdown_raw_date_inds]
  gdown_raw_data <- x_raw[gdown_raw_date_inds] # get the raw data associated with the gdown segment (this is the pre-filled, despiked version)
  gdown_smooth_data <- x_smooth[gdown_smooth_date_inds] # get the smoothed values associated with each raw data value
  
  gdown_numObs <- sum(!is.na(gdown_raw_data)) 
  
  
  if (gup_numObs == 0 | gdown_numObs == 0) {return(rep(NA,9))}
  
  
  ###Get the observation density for each period
  #This approach counts snow filled values as good values, since snow images are valuable for pinning down dormant period
  ###ind the biggest gap between images
  gup_seg_rsquared <- 1 - (sum((gup_raw_data - gup_smooth_data)^2, na.rm=T) / sum((gup_raw_data - mean(gup_raw_data, na.rm=T))^2, na.rm=T))
  gup_seg_rsquared[is.infinite(gup_seg_rsquared)] <- NA
  gup_maxgap <- max(diff(c(smooth_dates[seg[1]],raw_dates_gup[!is.na(gup_raw_data)],smooth_dates[seg[2]])))
  gup_maxgap_frac <- gup_maxgap / (seg[2] - seg[1]) 
    
  gdown_seg_rsquared <- 1 - (sum((gdown_raw_data - gdown_smooth_data)^2, na.rm=T) / sum((gdown_raw_data - mean(gdown_raw_data, na.rm=T))^2, na.rm=T))
  gdown_seg_rsquared[is.infinite(gdown_seg_rsquared)] <- NA
  gdown_maxgap <- max(diff(c(smooth_dates[seg[2]],raw_dates_gdown[!is.na(gdown_raw_data)],smooth_dates[seg[3]])))
  gdown_maxgap_frac <- gdown_maxgap / (seg[3] - seg[2]) 
  
  return(c(seg_min, seg_max, seg_int, 
           gup_seg_rsquared, gup_numObs, gup_maxgap_frac, 
           gdown_seg_rsquared, gdown_numObs, gdown_maxgap_frac))
}





#----------------------------------------------------------
#Developed by Josh Gray, updated by Douglas Bolton for HLS
#----------------------------------------------------------
GetSegMetricsLight <- function(seg, smooth_dates, raw_dates){
  if(any(is.na(seg))){return(NA,NA)}
  
  # organize greenup segment
  ######################################
  gup_raw_date_inds <- which(raw_dates >= smooth_dates[seg[1]] & raw_dates <= smooth_dates[seg[2]]) # indices in raw data of gup segment
  raw_dates_gup <- raw_dates[gup_raw_date_inds]

  # organize greendown segment
  ######################################
  gdown_raw_date_inds <- which(raw_dates >= smooth_dates[seg[2]] & raw_dates <= smooth_dates[seg[3]]) # indices in raw data of gdown segment
  raw_dates_gdown <- raw_dates[gdown_raw_date_inds]
  
  ###Get the observation density for each period
  #This approach counts snow filled values as good values, since snow images are valuable for pinning down dormant period
  ###ind the biggest gap between images
  gup_maxgap <- max(diff(c(smooth_dates[seg[1]],raw_dates_gup,smooth_dates[seg[2]])))
  gup_maxgap_frac <- gup_maxgap / (seg[2] - seg[1]) 
  
  gdown_maxgap <- max(diff(c(smooth_dates[seg[2]],raw_dates_gdown,smooth_dates[seg[3]])))
  gdown_maxgap_frac <- gdown_maxgap / (seg[3] - seg[2]) 
  
  return(c(gup_maxgap_frac, gdown_maxgap_frac))
}



#----------------------------------------------------------
#When a cycle is not detected, return a subset of metrics for the calendar year
#For now, returning evi minimum and evi maximum
#Douglas Bolton
#----------------------------------------------------------
annualMetrics <- function(smooth_evi2,pred_dates,yr,pheno_pars) {
  
  out <- matrix(NA,pheno_pars$numLyrs)
  try({
  ind <- as.numeric(format(pred_dates,'%Y')) == yr
  evi2_inyear  <- smooth_evi2[ind]
  out[1] <- 0    #Zero = No cycles detected
  out[44] <- min(evi2_inyear,na.rm=T)  * 10000    
  out[45] <- max(evi2_inyear,na.rm=T)  * 10000
  },silent=T)
  
  return(out)}




#----------------------------------------------------------
#Make image composites for specific dates from a smoothed time-series
#Does all bands at the same time
#Douglas Bolton
#----------------------------------------------------------
MakeComposites <- function(compositeDates,pred_dates,smooth_b1, smooth_b2, smooth_b3, smooth_b4, smooth_b5, smooth_b6) {

  days <- which(pred_dates %in% compositeDates)
  
  #Ensure that all days were found. If not, don't return any, as this means something is incorrect
  if (length(days) != length(compositeDates)) {return(NA)}
  
  out <- rbind(smooth_b1[days],smooth_b2[days],smooth_b3[days],
               smooth_b4[days],smooth_b5[days],smooth_b6[days])
  
  return(out)
  
}










#---------------------------------------------------------------------
#Calculate pheno metrics for each pixel
#This version using alternate years to gap fill
#Code adapted by Douglas Bolton
#Based on MODIS C6 algorithm developed by Josh Gray
#---------------------------------------------------------------------
DoPhenologyHLS_Climatology <- function(b1, b2, b3, b4, b5, b6, evi2, snowPix, dates, yrs, splineStart, numDaysFit, pheno_pars){
  
  
  numYrs <- length(yrs)
  daysVec <- 1:numDaysFit
  vecLength <- numDaysFit*numYrs
  
  
  #Despike, calculate dormant value, fill snow with dormant value, despike poorly fit snow values
  log <- try({
    #Despike
    spikes <- CheckSpike_MultiBand(b1/10000,b3/10000,evi2, dates, pheno_pars)
    
    #Set despiked values to NA
    evi2[spikes] <- NA
    b1[spikes] <- NA; b2[spikes] <- NA; b3[spikes] <- NA
    b4[spikes] <- NA; b5[spikes] <- NA; b6[spikes] <- NA
    
    dormIms <- dates >= pheno_pars$dormStart & dates <= pheno_pars$dormEnd
    
    evi2_dorm <- quantile(evi2[dormIms],probs=pheno_pars$dormantQuantile,na.rm=T)   #Calc evi2 dormant value
    snowPix <- Screen_SnowFills(evi2,evi2_dorm,snowPix,dates)              #Screen poorly filled snow values
    evi2[snowPix] <- evi2_dorm                                             #Fill remaining snow values with dormant value  
    
    #now calculate dormancy values and fill each band
    b1_dorm <- quantile(b1[dormIms],probs=pheno_pars$dormantQuantile,na.rm=T); b1[snowPix] <- b1_dorm
    b2_dorm <- quantile(b2[dormIms],probs=pheno_pars$dormantQuantile,na.rm=T); b2[snowPix] <- b2_dorm
    b3_dorm <- quantile(b3[dormIms],probs=pheno_pars$dormantQuantile,na.rm=T); b3[snowPix] <- b3_dorm
    b4_dorm <- quantile(b4[dormIms],probs=pheno_pars$dormantQuantile,na.rm=T); b4[snowPix] <- b4_dorm
    b5_dorm <- quantile(b5[dormIms],probs=pheno_pars$dormantQuantile,na.rm=T); b5[snowPix] <- b5_dorm
    b6_dorm <- quantile(b6[dormIms],probs=pheno_pars$dormantQuantile,na.rm=T); b6[snowPix] <- b6_dorm
    
    #First, we will fit splines to each year invidually
    #To line up observations from each year, we will create a matrix for evi2 and each band (numDaysFit x numYears)
    
    #Assign base Weights
    baseW <- matrix(1,length(evi2))
    baseW[snowPix == 1] <- pheno_pars$snowWeight
    
    #Option to not include observations from buffer period
    #climInds <- matrix(TRUE,numDaysFit)
    #climInds[pheno_pars$splineBuffer:(numDaysFit-pheno_pars$splineBuffer)] <- FALSE
    
    splineEnd <- splineStart+(numDaysFit-1)
    
    #Single spline through all years
    fullDates <- seq(splineStart[1], splineEnd[numYrs], by="day")
    smoothed <- Just_Smooth(evi2, dates, fullDates, baseW, evi2_dorm, pheno_pars)
    
    gDates <- dates[!is.na(evi2)]
    
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
    b1Mat <- smoothMat; b2Mat <- smoothMat; b3Mat <- smoothMat
    b4Mat <- smoothMat; b5Mat <- smoothMat; b6Mat <- smoothMat
    baseWeights <- smoothMat
    for (y in 1:numYrs) {
      dateRange <- dates >= splineStart[y] & dates <= splineEnd[y] & !is.na(evi2)   
      fillDs <- as.numeric(dates[dateRange] - splineStart[y]) + 1
      
      smoothMat[,y] <- smoothed[fullDates >= splineStart[y] & fullDates <= splineEnd[y]]
      baseWeights[fillDs,y] <- baseW[dateRange]
      fillMat[fillDs,y] <- evi2[dateRange]
      b1Mat[fillDs,y] <- b1[dateRange]; b2Mat[fillDs,y] <- b2[dateRange]; b3Mat[fillDs,y] <- b3[dateRange]
      b4Mat[fillDs,y] <- b4[dateRange]; b5Mat[fillDs,y] <- b5[dateRange]; b6Mat[fillDs,y] <- b6[dateRange] }
    
    
    xs <- rep(daysVec,numYrs)
    ys <- matrix(fillMat,vecLength)
    ysGood <- !is.na(ys)
    baseW <- matrix(baseWeights,vecLength)   #Base Weights are 1=clear observation, 0.5=snow-filled
    
    b1Mat <- matrix(b1Mat,vecLength); b2Mat <- matrix(b2Mat,vecLength) 
    b3Mat <- matrix(b3Mat,vecLength); b4Mat <- matrix(b4Mat,vecLength) 
    b5Mat <- matrix(b5Mat,vecLength); b6Mat <- matrix(b6Mat,vecLength) 
    
    
  },silent=TRUE)
  #If there is an error despiking, return NAs
  if(inherits(log, "try-error")){return(matrix(NA,pheno_pars$numLyrs*numYrs))}   
  
  #Loop through years, compare spline to other years, weight each year based on similiarity, fit spline, calculate phenology
  outAll=c()
  for (y in 1:numYrs) {
    log <- try({
      
      pred_dates <- seq(splineStart[y], splineEnd[y], by="day")  
      #Only compare on dates that have splined data in target year
      ind <- !is.na(smoothMat[,y])
      sub_evi2 <- smoothMat[ind,]
      numGoodDays <- colSums(!is.na(sub_evi2))  #how many days actually have splined data?
      
      #What approach to use for weighting
      ######
      #Calculate euclidean distance
      eucl <- colSums((sub_evi2 - sub_evi2[,y])^2,na.rm=T)^0.5 
      
      #Now calculate euclidean distance assuming the average through the year
      #Scale euculidean distances between this value and a perfect fit (0 to 1)
      theAvg <- matrix(mean(sub_evi2[,y],na.rm=T),length(sub_evi2[,y]),numYrs,byrow=T)
      theAvg[is.na(sub_evi2)] <- NA   #only calculate for days that have data
      
      max_eucl <- colSums((theAvg - sub_evi2[,y])^2,na.rm=T)^0.5    #calculate eucidean distance for this case
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
      smooth_evi2 <- Just_Smooth(ys[theInds], xs_sub, daysVec, w_sub, evi2_dorm, pheno_pars)  #Fit spline
      
      #Fit phenology
      peaks <- FindPeaks(smooth_evi2)
      if (all(is.na(peaks))) {outAll <- c(outAll,annualMetrics(smooth_evi2,pred_dates,yrs[y],pheno_pars));next}   #If no peaks, calc annual metrics, and move to next year
      
      #Find full segments
      full_segs <- GetSegs(peaks, smooth_evi2, pheno_pars)
      if (is.null(full_segs)) {outAll <- c(outAll,annualMetrics(smooth_evi2,pred_dates,yrs[y],pheno_pars));next}  #If no valid segments, calc annual metrics, and move to next year
      
      #Only keep segments with peaks within year *****
      kept_peaks <- unlist(full_segs, use.names=F)  
      kept_peaks  <- kept_peaks[seq(2, length(kept_peaks), by=3)]    #extract peaks. If there are no peaks, an error will thrown, and will move to next year
      full_segs <- full_segs[as.numeric(format(pred_dates[kept_peaks],'%Y')) == yrs[y]]  #check if peaks are in the year
      
      #Get PhenoDates
      pheno_dates <- GetPhenoDates(full_segs, smooth_evi2, pred_dates, pheno_pars)
      phen <- unlist(pheno_dates, use.names=F)
      if (all(is.na(phen))) {outAll <- c(outAll,annualMetrics(smooth_evi2,pred_dates,yrs[y],pheno_pars));next} #If no dates detected, calc annual metrics, and move to next year
      
      #Get the segment metrics
      filled_evi2 <- fillMat[,y]
      filled_evi2[baseWeights[,y] < 1] <- NA    #If weight is less than 1, implies it is a snow-fill, and we don't want to count snow-filled as a valid observation. So set to NA.
      seg_metrics <- lapply(full_segs, GetSegMetrics, smooth_evi2, filled_evi2[!is.na(filled_evi2)], pred_dates, pred_dates[!is.na(filled_evi2)])
      
      seg_metricsFill <- lapply(full_segs, GetSegMetricsLight, daysVec, sort(xs[ysGood & weights > 0 & baseW == 1]))
      
      #Unlist and scale the seg metrics
      un <- unlist(seg_metrics, use.names=F)
      ln <- length(un)
      seg_min <- un[seq(1, ln, by=9)] * 10000
      seg_max <- un[seq(2, ln, by=9)] * 10000
      seg_int <- un[seq(3, ln, by=9)] * 100
      gup_rsq <- un[seq(4, ln, by=9)] * 10000
      gup_numObs <- un[seq(5, ln, by=9)]
      gup_maxgap_frac <- un[seq(6, ln, by=9)] * 100
      gdown_rsq <- un[seq(7, ln, by=9)] * 10000
      gdown_numObs <- un[seq(8, ln, by=9)]
      gdown_maxgap_frac <- un[seq(9, ln, by=9)] * 100
      
      un <- unlist(seg_metricsFill, use.names=F)
      ln <- length(un)      
      gup_maxgap_frac_filled <- un[seq(1, ln, by=2)] * 100
      gdown_maxgap_frac_filled <- un[seq(2, ln, by=2)] * 100
      
      numRecords <- length(seg_min)  #how many cycles were recorded
      
      naCheck <- is.na(seg_min)
      numCyc <- sum(naCheck == 0)  #how many cycles have good data (seg metrics has valid observations)
      
      #If no cycles have good data, record NA output and move to next
      if (numCyc == 0) {outAll <- c(outAll,annualMetrics(smooth_evi2,pred_dates,yrs[y],pheno_pars));next}
      
      
      # Fit spline to individual bands. Waiting until now so that we don't fit for pixels with no cycle
      smooth_b1 <- Just_Smooth(b1Mat[theInds], xs_sub, daysVec, w_sub, b1_dorm, pheno_pars)
      smooth_b2 <- Just_Smooth(b2Mat[theInds], xs_sub, daysVec, w_sub, b2_dorm, pheno_pars)
      smooth_b3 <- Just_Smooth(b3Mat[theInds], xs_sub, daysVec, w_sub, b3_dorm, pheno_pars)
      smooth_b4 <- Just_Smooth(b4Mat[theInds], xs_sub, daysVec, w_sub, b4_dorm, pheno_pars)
      smooth_b5 <- Just_Smooth(b5Mat[theInds], xs_sub, daysVec, w_sub, b5_dorm, pheno_pars)
      smooth_b6 <- Just_Smooth(b6Mat[theInds], xs_sub, daysVec, w_sub, b6_dorm, pheno_pars)

      
      if (numRecords == 1) {
        comp <- as.vector(MakeComposites(phen,pred_dates,smooth_b1, smooth_b2, smooth_b3, smooth_b4, smooth_b5, smooth_b6))
        
        #If only one cycle was recorded, report it, and fill NA values for second cycle
        out <- c(1, phen, comp, seg_min, seg_max, seg_int,
                 gup_rsq, gup_numObs, gup_maxgap_frac, gup_maxgap_frac_filled, 
                 gdown_rsq, gdown_numObs, gdown_maxgap_frac, gdown_maxgap_frac_filled, 
                 matrix(NA,pheno_pars$halfLyrs))
      } else {
        #If there are multiple cycles, sort by amplitude and report two highest amplitudes (highest amplitude first)
        seg_amp <- seg_max - seg_min
        theOrd <- order(seg_amp,decreasing=T)   
        
        phen1 <- phen[seq(theOrd[1], length(phen), by = numRecords)]
        phen2 <- phen[seq(theOrd[2], length(phen), by = numRecords)]
        
        comp1 <- as.vector(MakeComposites(phen1,pred_dates,smooth_b1, smooth_b2, smooth_b3, smooth_b4, smooth_b5, smooth_b6))
        comp2 <- as.vector(MakeComposites(phen2,pred_dates,smooth_b1, smooth_b2, smooth_b3, smooth_b4, smooth_b5, smooth_b6))
        
        if (naCheck[theOrd[2]]) {
          #If the second cycle did not have enough observations (seg_metrics = NA), only report the first cycle
          out <- c(numCyc, 
                   phen1, comp1, seg_min[theOrd[1]], seg_max[theOrd[1]], seg_int[theOrd[1]],
                   gup_rsq[theOrd[1]], gup_numObs[theOrd[1]], gup_maxgap_frac[theOrd[1]], gup_maxgap_frac_filled[theOrd[1]], 
                   gdown_rsq[theOrd[1]],gdown_numObs[theOrd[1]], gdown_maxgap_frac[theOrd[1]], gdown_maxgap_frac_filled[theOrd[1]],
                   matrix(NA,pheno_pars$halfFill))
          
        } else {
          #If the second cycle had enough observations, report both
          out <- c(numCyc, 
                   phen1, comp1, seg_min[theOrd[1]], seg_max[theOrd[1]], seg_int[theOrd[1]],
                   gup_rsq[theOrd[1]], gup_numObs[theOrd[1]], gup_maxgap_frac[theOrd[1]], gup_maxgap_frac_filled[theOrd[1]], 
                   gdown_rsq[theOrd[1]],gdown_numObs[theOrd[1]], gdown_maxgap_frac[theOrd[1]], gdown_maxgap_frac_filled[theOrd[1]],
                   phen2, comp2, seg_min[theOrd[2]], seg_max[theOrd[2]], seg_int[theOrd[2]],
                   gup_rsq[theOrd[2]], gup_numObs[theOrd[2]], gup_maxgap_frac[theOrd[2]], gup_maxgap_frac_filled[theOrd[2]], 
                   gdown_rsq[theOrd[2]], gdown_numObs[theOrd[2]], gdown_maxgap_frac[theOrd[2]], gdown_maxgap_frac_filled[theOrd[2]])
        }
      }
    },silent=TRUE)  #End of the try block
    if(inherits(log, "try-error")){outAll <- c(outAll,matrix(NA,pheno_pars$numLyrs))
    } else {outAll <- c(outAll,out);remove(out)}
  }
  return(outAll)
}






CreateQA <- function(tile,yr,phenoFold,waterMask,baseImage,pheno_pars) {
  
  dType <- 'INT1U';naVal <- 255
  
  NumCycles <- readGDAL(paste0(phenoFold,'y',yr,'/LSP_',tile,'_',yr,'_NumCycles.tif'),silent=T)$band1
  
  qaNames <- c('gupQA','gupQA_2','gdownQA','gdownQA_2')
  rNames <- c('gupRsq','gupRsq_2','gdownRsq','gdownRsq_2')
  mgNames <- c('gupMaxGap','gupMaxGap_2','gdownMaxGap','gdownMaxGap_2')
  mgfNames <- c('gupMaxGapFilled','gupMaxGapFilled_2','gdownMaxGapFilled','gdownMaxGapFilled_2')
  
  for (q in 1:length(qaNames)) {
    #GUP QA for cycle 1
    Rsq <- readGDAL(paste0(phenoFold,'y',yr,'/LSP_',tile,'_',yr,'_',rNames[q],'.tif'),silent=T)$band1 / 10000
    maxGap <- readGDAL(paste0(phenoFold,'y',yr,'/LSP_',tile,'_',yr,'_',mgNames[q],'.tif'),silent=T)$band1
    maxGapFill <- readGDAL(paste0(phenoFold,'y',yr,'/LSP_',tile,'_',yr,'_',mgfNames[q],'.tif'),silent=T)$band1

    Qual <- matrix(0,length(NumCycles))
    Qual[Rsq > pheno_pars$r2_qa & maxGap <= pheno_pars$maxGap_qa[1]]  <- 1     #High quality pixel
    Qual[Rsq > pheno_pars$r2_qa & (maxGap > pheno_pars$maxGap_qa[1] & maxGap <= pheno_pars$maxGap_qa[2])]  <- 2      #Moderate quality pixel
    Qual[Rsq > pheno_pars$r2_qa & (maxGap > pheno_pars$maxGap_qa[2] & maxGapFill <= pheno_pars$maxGap_qa[2])]  <- 3  #Poor quality, but a successful fill with other years
    Qual[Rsq <= pheno_pars$r2_qa | (maxGap > pheno_pars$maxGap_qa[2] & maxGapFill > pheno_pars$maxGap_qa[2])]  <- 4   #Poor quality, and unsuccessful fill with other years, or low R2 fit (< 0.75)
    Qual[NumCycles == 0 | is.na(Rsq)] <- 5                                                                          #No cycles detected for this pixel
    Qual[waterMask == 1] <- 6                                                                                        #Water pixel, phenology code wasn't run
    outName <- paste0(phenoFold,'y',yr,'/LSP_',tile,'_',yr,'_',qaNames[q],'.tif')
    rast <- setValues(baseImage,Qual)  #Convert to raster using base image
    writeRaster(rast,filename=outName,format='GTiff',overwrite=TRUE,datatype=dType,NAflag=naVal) 
  }
}

