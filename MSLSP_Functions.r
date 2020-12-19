#Functions for MuSLI_LSP - Multisource Land Imaging of Land Surface Phenology

#Functions written by Josh Gray, Douglas Bolton, and Eli Melaas
###############



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
  mask <- matrix(FALSE,length(qaM),1)
  snow <- matrix(FALSE,length(qaM),1)

  #Out of the image value for S30 is NA. Set mask to 1
  mask[is.na(qaM)] <- TRUE
  qaM[is.na(qaM)] <- 0
  
  #Out of the image value for L30 is 255. Set mask to 1
  mask[qaM == 255] <- TRUE
  qaM[qaM == 255] <- 0
  
  #If bits 6-7 are 11, 01 or 10 (high/average/low aerosol), don't label, just subtract values  (we aren't confident enough in these values yet)
  qaM[qaM >= 192] = qaM[qaM >= 192] - 192 
  qaM[qaM >= 128] = qaM[qaM >= 128] - 128 
  qaM[qaM >= 64] = qaM[qaM >= 64] - 64
  
  #check if water, but don't label mask! Could be a false water detection, and we don't want to remove those. We'll use ancillary data for water
  qaM[qaM >= 32] <- qaM[qaM >= 32] - 32
  
  #check if snow/ice, label snow layer as 1
  snow[(qaM >= 16)] <- TRUE
  qaM[qaM >= 16] <- qaM[qaM >= 16] - 16
  
  #if cloud shadow, adjacent to cloud, cloud, or cirrus, label as 1
  mask[(qaM >= 1)] <- TRUE
  
  #mask <- mask == 1
  #snow <- snow == 1
  
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
ApplyMask_QA_and_Fmask <-function(imgName, tile, waterMask, chunkStart, chunkEnd, params, deleteInput=F) {
  
  pheno_pars <- params$phenology_parameters   #Pull phenology parameters
  
  #Get sensor, names sorted
  imgName_strip = tail(unlist(strsplit(imgName,'/')),n = 1)
  sensor = unlist(strsplit(imgName_strip,'.',fixed = T))[2]
  tileTxt = unlist(strsplit(imgName_strip,'.',fixed = T))[3]
  date = unlist(strsplit(imgName_strip,'.',fixed = T))[4]
  outBase = paste0('HLS_',sensor,'_',tileTxt,'_',date)
  
  #Mask for clouds and snow. Use QA flags for Landsat. Impliment Fmask 4.0 for Sentinel
  ################################################################################
  if (sensor == 'L30' | params$setup$runFmaskSentinel == FALSE) {
    qaName = paste0('HDF4_EOS:EOS_GRID:',imgName, ':Grid:QA')
    maskList <- getMasks(qaName)
    mask <- maskList[[1]]
    snow <- maskList[[2]]
    snow[waterMask] <- FALSE   #Remove snow flags that are over water
    remove(maskList)
  } else {
    if (sensor == 'S30') {
      fmaskName <- paste0(params$dirs$fmaskDir, gsub('hdf','',imgName_strip), "Fmask.tif")
      if (!file.exists(fmaskName)) {
        RunFmask(imgName, tile, params$dirs$imgDir, params$dirs$fmaskDir, params$dirs$fmaskDir, params)}
    } else {
      #If it's S10 data, resample fmask. A 30m version must already exist!
      imgName_strip30 <- gsub(".S10",".S30",imgName_strip)
      fmaskName <- paste0(params$dirs$fmask10m, gsub('hdf','',imgName_strip), "Fmask.tif")
      if (!file.exists(fmaskName)) {
        fmaskIn <- paste0(params$dirs$fmaskDir, gsub('hdf','',imgName_strip30), "Fmask.tif")
        if (file.exists(fmaskIn)) {
          run <- try({system2("gdalwarp",paste("-overwrite -r near -ts 10980 10980 -of GTiff",fmaskIn,fmaskName),stdout=T,stderr=T)},silent=T)
        }
      }
    }
    
    logIt <- try({
      fmask <- readGDAL(fmaskName,silent=T)$band1     #If Fmask doesn't exist (didn't run), function will fail and error will be logged
      mask <- (fmask == 4) | (fmask == 2) | (fmask == 255)  #If cloud, cloud shadow, or no data
      snow  <- fmask == 3
      snow[waterMask] <- FALSE  #Remove snow flags that are over water (shouldn't be any though)
      remove(fmask)
    },silent=TRUE)
    
    if (inherits(logIt, 'try-error')) {
      qaName = paste0('HDF4_EOS:EOS_GRID:',imgName, ':Grid:QA')    #If there was an error with fmask, just use HLS QA
      maskList <- getMasks(qaName)
      mask <- maskList[[1]]
      snow <- maskList[[2]]
      snow[waterMask] <- FALSE   #Remove snow flags that are over water
      remove(maskList)
      print(paste('Using backup mask for',imgName))}
  }
  
  #Open the image bands
  ##############################################################################
  theBase <- paste0('HDF4_EOS:EOS_GRID:',imgName,':Grid:')
  if (sensor == 'L30') {
      bNames <- c('band02','band03','band04','band05','band06','band07')
      bands <- matrix(as.integer(0),length(mask),length(bNames))
      for (i in 1:length(bNames)) {bands[,i] <- as.integer(readGDAL(paste0(theBase, bNames[i]),silent=T)$band1*10000)}

  } else if (sensor == 'S30') {
     bNames <- c('B02','B03','B04','B8A','B11','B12','B05','B06','B07')
     bands <- matrix(as.integer(0),length(mask),length(bNames))
     for (i in 1:length(bNames)) {bands[,i] <- as.integer(readGDAL(paste0(theBase, bNames[i]),silent=T)$band1*10000)}
     
  } else if (sensor == 'S10') {
    #Need to get all bands to 10m first. Using Broad band NIR (10m) instead of Narrow (20m)
    bNames <- c('B08','B11','B12','B05','B06','B07')
    tifBase <- paste0(params$dirs$tempDir,gsub('hdf','',imgName_strip))
    for (bName in bNames) {
      inName <- paste0(theBase,bName);outName <- paste0(tifBase,bName,'.tif')
      run <- try({system2("gdalwarp",paste("-overwrite -r near -ts 10980 10980 -of GTiff",inName,outName),stdout=T,stderr=T)},silent=T)
    }

    fullNames <- c(paste0(theBase, c('B02','B03','B04')), 
                   paste0(tifBase,bNames,'.tif'))
    
    bands <- matrix(as.integer(0),length(mask),length(fullNames))
    for (i in 1:length(fullNames)) {bands[,i] <- as.integer(readGDAL(fullNames[i],silent=T)$band1*10000)}
    
    for (bName in bNames) {file.remove(paste0(tifBase,bName,'.tif'))}

  }
  

  #If we aren't topo correcting, then we need to do the NDMI check now. If we are topo correcting, we will do this check after topo correction
  #If NDMI > 0.5 AND the pixel is within 5 km of detected snow, then mask
  if (!params$topocorrection_parameters$topoCorrect) {
    if (sum(snow) > 0) {
      additional_snow_mask <- runSnowScreen(snow,bands[,4],bands[,5],params)
      mask[additional_snow_mask] <- TRUE
      remove(additional_snow_mask)}
  }
      

  bands[bands < 0] <- NA    #Remove negative reflectance values
  bands[mask,]  <-  NA
  bands[waterMask,]  <-  NA
  
  #Mask pixels that have missing data in ANY band. These are likely problem pixels (and we need most bands for despiking, kmeans, snow detection, etc).
  check <- rowSums(is.na(bands)) > 0
  bands[check,] <- NA
  
  #We remove potentially snowy observations using NDMI, but we are only going to label pixels that were actually detected as snow as snow
  bands[snow,]  <-  as.integer(pheno_pars$snowFillVal)    

  for (n in 1:length(chunkStart)) {
    mat <- as.integer(bands[chunkStart[n]:chunkEnd[n],])
    if (all(is.na(mat))) {next} #if there is no good data in the chunk, move to next chunk
    saveRDS(mat,paste0(params$dirs$chunkDir,'c',n,'/',outBase,'.Rds'))
  }
  #If requested, delete the input file once it is processed
  if (deleteInput == TRUE) {file.remove(imgName)}
}




RunFmask <-function(imgName, tile, auxFold, fmaskFold, outDir, params) {
  
  #Different approach depending on if we are on SCC or AWS
  #If on SCC, need to set environmental variables 
  if (params$setup$AWS_or_SCC == 'SCC') {
    system(paste('MCRROOT="/share/pkg/matlab/2018b/install"',
                 'LD_LIBRARY_PATH=.:${MCRROOT}/runtime/glnxa64',
                 'LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64',
                 'LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64',
                 'LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/opengl/lib/glnxa64',
                 'export LD_LIBRARY_PATH',
                 paste('eval',params$setup$fmaskFunction,tile,imgName,auxFold,fmaskFold),sep='\n'))
  
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
    ulx <- as.numeric(unlist(strsplit(info[pmatch('  ULX',info)],'='))[2]) + params$setup$image_res/2    #Adjust for matlab (needs center of pixel)
    uly <- as.numeric(unlist(strsplit(info[pmatch('  ULY',info)],'='))[2]) - params$setup$image_res/2    #Adjust for matlab (needs center of pixel)
    sza <- as.numeric(unlist(strsplit(info[pmatch('  MEAN_SUN_ZENITH_ANGLE(B01)',info)],'='))[2])
    saa <- as.numeric(unlist(strsplit(info[pmatch('  MEAN_SUN_AZIMUTH_ANGLE(B01)',info)],'='))[2])
    vza <- as.numeric(unlist(strsplit(info[pmatch('  MEAN_VIEW_ZENITH_ANGLE(B01)',info)],'='))[2])
    vaa <- as.numeric(unlist(strsplit(info[pmatch('  MEAN_VIEW_AZIMUTH_ANGLE(B01)',info)],'='))[2])  
  
    zone <- unlist(strsplit(info[pmatch('  HORIZONTAL_CS_NAME',info)],' '))
    zone <- zone[length(zone)]
    zone <- as.numeric(substr(zone,1,2))
  
    system(paste('MCR_CACHE_ROOT=$TMPDIR ',params$setup$fmaskFunction,tile,fileBase,auxFold,fmaskFold,zone,ulx,uly,sza,saa,vza,vaa)) 
    
    #Delete the temporary images
    for (b in bandNames) {file.remove(paste0(fileBase,b,'.tif'))}
  }
}





#
#---------------------------------------------------------------------
#Run additional snow screening
#If NDMI > 0.5 AND the pixel is within 5 km of detected snow, then mask pixel
#Douglas Bolton
#---------------------------------------------------------------------
runSnowScreen <- function(snow, nir, swir1, params) {
  pheno_pars <- params$phenology_parameters
  
  ndmi <- calcIndex(nir=nir,swirOne=swir1,whatIndex='ndmi')
  ndmiCheck = ndmi > pheno_pars$ndmiSnowThresh   #Find pixels that meet threshold (potential snow pixels)
  remove(ndmi)
  
  #We will only keep snow pixels when 50% of 5x5 pixel window is also snow (to remove speckle)
  imBlur  <- boxblur(as.cimg(matrix(snow,length(snow)^.5,length(snow)^.5)),pheno_pars$snowWindow,neumann = FALSE)    #neumann = FALSE assumes snow values of zero outside image edges
  snowCheck <- imBlur > pheno_pars$snowFraction & snow
  remove(imBlur)
  
  #Calcuate the distance from each pixel to the nearest pixel detected as snow by fmask or lasrc
  sD <- matrix(as.matrix(distance_transform(snowCheck,1)),length(snow),1)
  sCheck <- (sD*params$setup$image_res) < pheno_pars$distanceToSnow  #Find pixels where distance to snow is less than threshold (5 km)
  remove(snowCheck)
  remove(sD)
  
  #Find pixels meeting ndmi rule that are within 5km of detected snow pixels (or pixels already detected as snow)
  toMask <- (ndmiCheck & sCheck) | snow
  remove(ndmiCheck)
  remove(sCheck)
  
  #Buffer snow pixels by 100m
  mD <- matrix(as.matrix(distance_transform(as.cimg(matrix(toMask,length(snow)^.5,length(snow)^.5)),1)),length(snow),1)

  additional_snow_mask <- (mD*params$setup$image_res) < pheno_pars$snowBuffer

  #Mask these pixels. BUT we won't label as snow, as our confidence is lower, and we don't want to falsely fill dormant values. 

  return(additional_snow_mask)
}


#---------------------------------------------------------------------
#Calculate annual quantiles for spectral indices
#We will use these quantiles to create kmeans classes for topographic correction
#Douglas Bolton
#---------------------------------------------------------------------
getIndexQuantile <- function(chunk, numPix, yr, errorLog, params) {
  
  topo_pars <- params$topocorrection_parameters  #isolate topo parameters
  
  #Make empty output
  indexOut <- matrix(NA,numPix,length(topo_pars$topoVIs)*length(topo_pars$viQuantiles))
  
  log <- try({
  
    #Get all images to process
    ######################
    chunkFold <- paste0(params$dirs$chunkDir,'c',chunk,'/') 
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
        b2 <- matrix(NA,numPix,numImgs); b3 <- matrix(NA,numPix,numImgs)
        b4 <- matrix(NA,numPix,numImgs); b5 <- matrix(NA,numPix,numImgs)
        b6 <- matrix(NA,numPix,numImgs); b7 <- matrix(NA,numPix,numImgs)
        for (i in 1:numImgs) {
          imgData <- try(matrix(readRDS(paste0(chunkFold,imgList[i])),nrow=numPix),silent = TRUE)
          if (inherits(imgData, 'try-error')) {cat(paste('getIndexQuantile: Error for chunk',chunk,imgList[i]), file=errorLog, append=T);next} 
          imgData[imgData == params$phenology_parameters$snowFillVal] <- NA
          imgData <- imgData / 10000
          b2[,i] <- imgData[,1]; b3[,i] <- imgData[,2]; b4[,i] <- imgData[,3]
          b5[,i] <- imgData[,4]; b6[,i] <- imgData[,5]; b7[,i] <- imgData[,6]
        }
        
        #Loop through requested indices in topo_pars$topoVIs
        #For each index, loop through requested quantiles in topo_pars$viQuantiles
        column <- 0
        for (i in 1:length(topo_pars$topoVIs)) { 
          ind <- calcIndex(blue=b2,green=b3,red=b4,nir=b5,swirOne=b6,swirTwo=b7,whatIndex=topo_pars$topoVIs[i])
          ind[is.infinite(ind)] <- NA
          for (j in 1:length(topo_pars$viQuantiles)) {
            column <- column + 1
            indexOut[,column] <- rowQuantileC(ind,topo_pars$viQuantiles[j])
          }
        }
      }
    }
  },silent=T)
  
  return(indexOut)
}



#---------------------------------------------------------------------
#Run topographic correction
#Overwrite existing image chunks with new, corrected data
#Douglas Bolton
#---------------------------------------------------------------------
runTopoCorrection <- function(imgName, groups, slopeVals, aspectVals, chunkStart,chunkEnd, errorLog, params, writeImages=FALSE, slope=NULL) {
  
  pheno_pars <- params$phenology_parameters
  topo_pars <- params$topocorrection_parameters
  
  #Get sensor, names sorted
  imgName_strip <- tail(unlist(strsplit(imgName,'/')),n = 1)
  sensor <- unlist(strsplit(imgName_strip,'.',fixed = T))[2]
  tileTxt <- unlist(strsplit(imgName_strip,'.',fixed = T))[3]
  date <- unlist(strsplit(imgName_strip,'.',fixed = T))[4]
  outBase <- paste0('HLS_',sensor,'_',tileTxt,'_',date)

  #Get sun zenith and azimuth.
  #Some L8 scenes will have two values (two L8 images as inputs). We will take the average of these two values
  info <- gdalinfo(imgName)
  line <- info[pmatch('  MEAN_SUN_ZENITH_ANGLE',info)]
  line <- unlist(strsplit(line,'='))[2]
  zMean <- mean(as.numeric(unlist(strsplit(line,','))))
  sunzenith <- (pi/180) * zMean
  
  line = info[pmatch('  MEAN_SUN_AZIMUTH_ANGLE',info)]
  line <- unlist(strsplit(line,'='))[2]
  aMean <- mean(as.numeric(unlist(strsplit(line,','))))
  sunazimuth <- (pi/180) * aMean
  
  numChunks <- length(chunkEnd)
  
  if (sensor == 'L30') {bands <- matrix(NA,chunkEnd[numChunks],6)}
  if (sensor == 'S30' | sensor == 'S10') {bands <- matrix(NA,chunkEnd[numChunks],9)}
  
  #Reconstruct each image from the chunked data
  for (n in 1:numChunks) {
    fName <- paste0(params$dirs$chunkDir,'c',n,'/',outBase,'.Rds')
    if (file.exists(fName)) {
      numPix <- length(chunkStart[n]:chunkEnd[n])
      imgData <- try(matrix(readRDS(fName),nrow=numPix),silent = TRUE)
      if (inherits(imgData, 'try-error')) {cat(paste('runTopoCorrection: Error for chunk',n,outBase), file=errorLog, append=T);next} 
      bands[chunkStart[n]:chunkEnd[n],] <- imgData
      remove(imgData)
    }
  }
  
  
  if (!all(is.na(bands))) {                   #Only continue if there is actually good data
    
    snow <- bands[,1] == pheno_pars$snowFillVal       #Determine snow pixels
    snow[is.na(snow)] <- FALSE                        #Set NA snow values to FALSE (no snow)
    bands[bands == pheno_pars$snowFillVal] <-  NA     #Mask snow values
    bands  <- bands / 10000                         
    
    #If image outputs are requested, create a base name for plots
    plotBaseName <- NULL
    if (writeImages) {plotBaseName <- paste0(params$dirs$tempDir,'IL_plots/',outBase)}   #Base name for outputing scatterplots

    #Run topographic correction function by kmeans group
    corr <- topocorr_rotational_by_group(x=bands,groups=groups, 
                                           slope=slopeVals,aspect=aspectVals, 
                                           sunzenith = sunzenith, sunazimuth=sunazimuth,topo_pars,plotBaseName)

    #NDMI check
    #Now that we have topo corrected the imagery, let's check for potentail snow....
    #If NDMI > 0.5 AND the pixel is within 5 km of detected snow, then mask pixel
    if (sum(snow) > 0) {
      additional_snow_mask <- runSnowScreen(snow,corr[,4],corr[,5],params)
      corr[additional_snow_mask,]  <- NA}
      

    corr = round(corr * 10000)                  #Convert back to integer
    corr[corr < 0]  <-  NA                      #Remove negative reflectance values
    
    #Mask pixels that have missing data in ANY band. These are likely problem pixels (and we need most bands for despiking, kmeans, snow detection, etc).
    check <- rowSums(is.na(corr)) > 0
    corr[check,] <- NA
    
    #We removed potentially snowy observations using NDMI, but we are only going to label pixels that were actually detected as snow (this impacts dormant filling)
    corr[snow,]  <-  pheno_pars$snowFillVal   
    corr <- round(corr)
    
    #Overwrite existing chunks with new corrected data
    for (n in 1:numChunks) {
      fName <- paste0(params$dirs$chunkDir,'c',n,'/',outBase,'.Rds')
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
      afterName <- paste0(params$dirs$tempDir,'/sampleImages/',outBase,'_corrected.tif')
      writeRaster(after,filename=afterName,format='GTiff',overwrite=TRUE,datatype=dType)
      
      ##Write out an uncorrected version
      bands = round(bands * 10000)
      bands <- round(bands)
      before <- stack(setValues(slope,bands[,5]),setValues(slope,bands[,4]),setValues(slope,bands[,3])) #SWIR/NIR/RED composite
      beforeName <- paste0(params$dirs$tempDir,'/sampleImages/',outBase,'.tif')
      writeRaster(before,filename=beforeName,format='GTiff',overwrite=TRUE,datatype=dType) 
    }
  }   
}



#---------------------------------------------------------------------
#Run Phenology code for an image chunk
#Write phenology results for each chunk to disk
#Douglas Bolton
#---------------------------------------------------------------------
runPhenoChunk <- function(chunk, numPix, imgYrs, phenYrs, errorLog, params) {

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
    
    
    #Loop through each pixel and estimate phenometrics
    pheno_mat <- matrix(NA,numPix,pheno_pars$numLyrs*length(phenYrs))
    for (i in 1:numPix) {pheno_mat[i,] <- DoPhenologyHLS(b2[i,],  b3[i,],  b4[i,],  b5[i,],  b6[i,], b7[i,],  vi[i,],
                                                         snowPix[i,],dates, imgYrs, phenYrs, splineStart, splineEnd, numDaysFit, pheno_pars)}

    pheno_mat <- round(pheno_mat)
    
    #Write results to disk (.Rds files)
    for (y in 1:length(phenYrs)) {
      yr <- phenYrs[y]
      ind <- ((y-1)*pheno_pars$numLyrs+1):(y*pheno_pars$numLyrs)
      subTab <- pheno_mat[,ind]
      for (i in 1:pheno_pars$numLyrs) {saveRDS(as.integer(subTab[,i]),paste0(outFold,'y',yr,'/lyr',i,'/c',chunk,'.Rds'))}}
}



#---------------------------------------------------------------------
#Correct reflectance for topography
#Developed by Eli Melaas
#Rewritten by Douglas Bolton to only run rotational model (Tan et al. 2010) and to have a stratified sample of pixels based on IL
#---------------------------------------------------------------------
topocorr_rotational_by_group <-function(x, groups, slope, aspect, sunzenith, sunazimuth, topo_pars, plotBaseName, IL.epsilon=0.000001) {
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
  numPer <- round(topo_pars$numSamples/topo_pars$numILclass)
  
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
      breaks <- seq(quantile(ilSub,prob=0.02),quantile(ilSub,prob=0.98),length.out=(topo_pars$numILclass+1))    #Using 2 and 98% percentiles to reduce influece of outliers
      breaks[1] <- min(ilSub); breaks[length(breaks)] <- max(ilSub)
      ilGroup <- try({as.numeric(cut(ilSub,breaks=breaks,labels=1:topo_pars$numILclass,include.lowest=T))},silent=T)  #Group pixels
      if (inherits(ilGroup, 'try-error')) {next}    #If pixels can't be grouped, move on
      pixID <- 1:length(ilGroup)
      pixToSample <- matrix(0,length(pixID))
      for (i in 1:topo_pars$numILclass) {
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
CheckSpike_MultiBand <- function(blue,red,vi, dates, pheno_pars){
  
  blue_og <- blue # preserve original vector
  red_og <- red # preserve original vector
  vi_og <- vi # preserve original vector
  dates_og <- as.numeric(dates)
  
  good <- !is.na(blue_og) & !is.na(red_og) & !is.na(vi_og)
  
  x_outs <- matrix(F, length(blue_og)) # create the outlier output vector
  count <- 0 
  while (count < pheno_pars$maxDespikeIterations) {
    count <- count + 1
    
    bS <- blue_og[good] # subset to non missing values
    rS <- red_og[good]
    eS <- vi_og[good] 
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
    
    #look for negative spikes in vi
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
Screen_SnowFills <- function(x, dormVal, snowPix, dates, pheno_pars){
  #Reduce vectors to valid obervations and snow pixels
  sORg<- !is.na(x) | snowPix==1
  snowSub <- snowPix[sORg]
  xSub <- x[sORg]
  d <- as.numeric(dates[sORg])
  
  padded <- c(0,roll_sum(snowSub,3),0)   #Get rolling sum of snow detections
  padded2 <- c(3,roll_max(padded,3),3)   #Get rolling max of the sums. Setting first and last to 3 keeps these values if they are snow
  
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
Smooth_VI <- function(x, dates, pred_dates, weights, pheno_pars, dormant_value) {
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


Smooth_Bands <- function(x, dates, pred_dates, weights, pheno_pars){
  #Get index of pixels with good values
  ind <- !is.na(x)  
  # smooth with a spline to get continuous daily series
  spl <- smooth.spline(dates[ind], x[ind], spar=pheno_pars$splineSpar, w=weights[ind])
  # weighted version
  xSmooth <- predict(spl, as.numeric(pred_dates))$y
  
  # determine upper and lower bound of good data (set in json file)
  lowBound <- quantile(x, probs=pheno_pars$bandLimits[1],na.rm=T) 
  upperBound <- quantile(x, probs=pheno_pars$bandLimits[2],na.rm=T) 
  
  xSmooth[xSmooth < lowBound] <- lowBound
  xSmooth[xSmooth > upperBound] <- upperBound
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
#Get phenology dates from segments. Also pull the peak date
#Josh Gray. Updated by Douglas Bolton to include peak date and cleaned
#----------------------------------------------------------
GetPhenoDates <- function(segs, x, dates, pheno_pars){
  pheno_dates <- list()
  
  #Pull greenup dates
  for(gup_thresh in pheno_pars$gup_threshes){
    pheno_dates <- c(pheno_dates, list(dates[unlist(lapply(segs, GetSegThresh, x, gup_thresh, gup=T), use.names=F)]))
  }
  
  #Pull peak dates
  pheno_dates <- c(pheno_dates, list(dates[sapply(segs, "[[", 2)]))
  
  #Pull greendown dates
  for(gdown_thresh in pheno_pars$gdown_threshes){
    pheno_dates <- c(pheno_dates, list(dates[unlist(lapply(segs, GetSegThresh, x, gdown_thresh, gup=F), use.names=F)]))
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
      gdown_thresh_index <- GetThresh(gdown_thresh, x[seg[2]:seg[3]], first_greater=F, gup=F)
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
  seg_amp <- seg_max - seg_min
  
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
  
  return(c(seg_amp, seg_max, seg_int, 
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
#For now, returning evi amplitude and evi maximum
#Douglas Bolton
#----------------------------------------------------------
annualMetrics <- function(smoothed_vi, pred_dates, filled_vi, baseWeight, yr, pheno_pars) {
  
  out <- matrix(NA,pheno_pars$numLyrs)
  try({
    inyear <- as.numeric(format(pred_dates,'%Y')) == yr
    vi_inyear  <- smoothed_vi[inyear]
    seg_min <- min(vi_inyear,na.rm=T) 
    seg_max <- max(vi_inyear,na.rm=T) 

    out[pheno_pars$loc_numCycles] <- 0                        #Zero cycles detected
    out[pheno_pars$loc_max] <- seg_max  * 10000               #Cycle maximum
    out[pheno_pars$loc_amp] <- (seg_max - seg_min) * 10000    #Cycle amplitude
  

    #And get calendar year metrics with snow counted as good
    ind <- !is.na(filled_vi) & inyear
    out[pheno_pars$loc_numObs_count_snow]<- sum(ind)
    out[pheno_pars$loc_maxGap_count_snow] <- as.numeric(max(diff(c( min(pred_dates[inyear]), pred_dates[ind], max(pred_dates[inyear])))))
  
    #Now get the full segment metrics, not counting snow and not counting gap filled
    ind <- ind & baseWeight == 1   #If weight is less than 1, implies it is a snow-fill, and we don't want to count snow-filled as a valid observation. 
    out[pheno_pars$loc_numObs] <- sum(ind)
    out[pheno_pars$loc_maxGap] <- as.numeric(max(diff(c( min(pred_dates[inyear]), pred_dates[ind], max(pred_dates[inyear])))))
    },silent=T)
  return(out)}



#----------------------------------------------------------
#Make image composites for specific dates from a smoothed time-series
#Does all bands at the same time
#Douglas Bolton
#----------------------------------------------------------
MakeComposites <- function(compositeDates,pred_dates,smooth_b2, smooth_b3, smooth_b4, smooth_b5, smooth_b6, smooth_b7) {

  days <- match(compositeDates,pred_dates)
  
  #Ensure that all days were found. If not, don't return any, as this means something is incorrect
  if (length(days) != length(compositeDates)) {return(NA)}
  
  out <- rbind(smooth_b2[days],smooth_b3[days],smooth_b4[days],
               smooth_b5[days],smooth_b6[days],smooth_b7[days])
  
  return(out)
  
}





calculateWeights <- function(smoothMat_Masked, numDaysFit, numYrs, pheno_pars) {
  outWeights <- array(0, dim = c(numDaysFit, numYrs, numYrs))
  for (y in 1:numYrs) {
    #Only compare on dates that have splined data in target year
    ind <- !is.na(smoothMat_Masked[,y])
    sub_vi <- smoothMat_Masked[ind,]
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
    weight <- pheno_pars$maxWeight * scaled_eucl
    weight[numGoodDays < pheno_pars$minDaysForSplineComparison]  <- 0
    weight[is.na(weight)] <- 0
    weight[is.infinite(weight)] <- 0
    outWeights[,,y] <- matrix(weight,numDaysFit,numYrs,byrow=T)
  }
  return(outWeights)
}










#---------------------------------------------------------------------
#Calculate pheno metrics for each pixel
#This version using alternate years to gap fill
#Code adapted by Douglas Bolton
#Based on MODIS C6 algorithm developed by Josh Gray
#---------------------------------------------------------------------
DoPhenologyHLS <- function(b2, b3, b4, b5, b6, b7, vi, snowPix, dates, imgYrs, phenYrs, splineStart, splineEnd, numDaysFit, pheno_pars){
  
  
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
    if(inherits(log, "try-error")){return(matrix(NA,pheno_pars$numLyrs*length(phenYrs)))}   
  
  outAll=c()
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
      
      
      #Fit phenology
      peaks <- FindPeaks(smoothed_vi)
      if (all(is.na(peaks))) {outAll <- c(outAll,annualMetrics(smoothed_vi,pred_dates,fillMat[,y], baseWeights[,y], yrs[y],pheno_pars));next}   #If no peaks, calc annual metrics, and move to next year
      
      #Find full segments
      full_segs <- GetSegs(peaks, smoothed_vi, pheno_pars)
      if (is.null(full_segs)) {outAll <- c(outAll,annualMetrics(smoothed_vi,pred_dates,fillMat[,y], baseWeights[,y], yrs[y],pheno_pars));next}  #If no valid segments, calc annual metrics, and move to next year
      
      #Only keep segments with peaks within year *****
      full_segs <- full_segs[inYear[sapply(full_segs, "[[", 2)] ]  #check if peaks are in the year
      
      #Get PhenoDates
      pheno_dates <- GetPhenoDates(full_segs, smoothed_vi, pred_dates, pheno_pars)
      phen <- unlist(pheno_dates, use.names=F)
      if (all(is.na(phen))) {outAll <- c(outAll,annualMetrics(smoothed_vi,pred_dates,fillMat[,y], baseWeights[,y], yrs[y],pheno_pars));next} #If no dates detected, calc annual metrics, and move to next year
      
      
      #Get metrics that describe the segments and the year
      
      
      #First, get metrics counting gap filled observations as "good" observations
      seg_metricsFill <- lapply(full_segs, GetSegMetricsLight, daysVec, sort(xs_sub))
      un <- unlist(seg_metricsFill, use.names=F)
      ln <- length(un)      
      gup_maxgap_frac_filled <- un[seq(1, ln, by=2)] * 100
      gdown_maxgap_frac_filled <- un[seq(2, ln, by=2)] * 100
      
      
      #Second, get segment metrics with snow observations counted as "good" observations
      filled_vi <- fillMat[,y]
      seg_metricsFill <- lapply(full_segs, GetSegMetricsLight, daysVec, daysVec[!is.na(filled_vi)])
      un <- unlist(seg_metricsFill, use.names=F)
      ln <- length(un)      
      gup_maxgap_frac_count_snow <- un[seq(1, ln, by=2)] * 100
      gdown_maxgap_frac_count_snow <- un[seq(2, ln, by=2)] * 100
      
      #And get calendar year metrics with snow counted as good
      numObs_count_snow <- sum(!is.na(filled_vi) & inYear)
      maxGap_annual_count_snow <- max(diff(c( pheno_pars$splineBuffer+1, daysVec[!is.na(filled_vi) & inYear], 365+pheno_pars$splineBuffer))) 
      
      
      #Now get the full segment metrics, not counting snow and not counting gap filled
      filled_vi[baseWeights[,y] < 1] <- NA    #If weight is less than 1, implies it is a snow-fill, and we don't want to count snow-filled as a valid observation. So set to NA.
      numObs <- sum(!is.na(filled_vi) & inYear)   #Number of observations in year
      maxGap_annual <- max(diff(c( pheno_pars$splineBuffer+1, daysVec[!is.na(filled_vi) & inYear], 365+pheno_pars$splineBuffer)))  #Max gap (in days) during year
      seg_metrics <- lapply(full_segs, GetSegMetrics, smoothed_vi, filled_vi[!is.na(filled_vi)], pred_dates, pred_dates[!is.na(filled_vi)]) #full segment metrics
      
      
      #Unlist and scale the seg metrics
      un <- unlist(seg_metrics, use.names=F)
      ln <- length(un)
      seg_amp <- un[seq(1, ln, by=9)] * 10000
      seg_max <- un[seq(2, ln, by=9)] * 10000
      seg_int <- un[seq(3, ln, by=9)] * 100
      gup_rsq <- un[seq(4, ln, by=9)] * 10000
      gup_numObs <- un[seq(5, ln, by=9)]
      gup_maxgap_frac <- un[seq(6, ln, by=9)] * 100
      gdown_rsq <- un[seq(7, ln, by=9)] * 10000
      gdown_numObs <- un[seq(8, ln, by=9)]
      gdown_maxgap_frac <- un[seq(9, ln, by=9)] * 100
      
      
      numRecords <- length(seg_amp)  #how many cycles were recorded
      
      naCheck <- is.na(seg_amp)
      numCyc <- sum(naCheck == 0)  #how many cycles have good data (seg metrics has valid observations)
      
      #If no cycles have good data, record NA output and move to next
      if (numCyc == 0) {outAll <- c(outAll,annualMetrics(smoothed_vi,pred_dates,fillMat[,y], baseWeights[,y], yrs[y],pheno_pars));next}
      
      
      # Fit spline to individual bands. Waiting until now so that we don't fit for pixels with no cycle
      # Since this takes the most time, there is an option to skip this and just fill NAs
      if (pheno_pars$doComposites) {
        smooth_b2 <- Smooth_Bands(b2Mat[theInds], xs_sub, daysVec, w_sub, pheno_pars)
        smooth_b3 <- Smooth_Bands(b3Mat[theInds], xs_sub, daysVec, w_sub, pheno_pars)
        smooth_b4 <- Smooth_Bands(b4Mat[theInds], xs_sub, daysVec, w_sub, pheno_pars)
        smooth_b5 <- Smooth_Bands(b5Mat[theInds], xs_sub, daysVec, w_sub, pheno_pars)
        smooth_b6 <- Smooth_Bands(b6Mat[theInds], xs_sub, daysVec, w_sub, pheno_pars)
        smooth_b7 <- Smooth_Bands(b7Mat[theInds], xs_sub, daysVec, w_sub, pheno_pars)
      } else {
        smooth_b2 <- matrix(NA,daysVec);smooth_b3 <- matrix(NA,daysVec);smooth_b4 <- matrix(NA,daysVec)
        smooth_b5 <- matrix(NA,daysVec);smooth_b6 <- matrix(NA,daysVec);smooth_b7 <- matrix(NA,daysVec)}
      
      
      if (numRecords == 1) {
        comp <- as.vector(MakeComposites(phen,pred_dates,smooth_b2, smooth_b3, smooth_b4, smooth_b5, smooth_b6, smooth_b7))
        
        #If only one cycle was recorded, report it, and fill NA values for second cycle
        out <- c(1, 
                 phen, seg_max, seg_amp, seg_int, comp,
                 gup_numObs,   gup_maxgap_frac,   gup_maxgap_frac_count_snow,   gup_maxgap_frac_filled,   gup_rsq,   
                 gdown_numObs, gdown_maxgap_frac, gdown_maxgap_frac_count_snow, gdown_maxgap_frac_filled, gdown_rsq, 
                 matrix(NA,pheno_pars$numLyrsCycle2),
                 numObs, numObs_count_snow, maxGap_annual, maxGap_annual_count_snow)
        
      } else {
        #If there are multiple cycles, sort by amplitude and report two highest amplitudes (highest amplitude first)
        theOrd <- order(seg_amp,decreasing=T)   
        
        phen1 <- phen[seq(theOrd[1], length(phen), by = numRecords)]
        phen2 <- phen[seq(theOrd[2], length(phen), by = numRecords)]
        
        comp1 <- as.vector(MakeComposites(phen1,pred_dates,smooth_b2, smooth_b3, smooth_b4, smooth_b5, smooth_b6, smooth_b7))
        comp2 <- as.vector(MakeComposites(phen2,pred_dates,smooth_b2, smooth_b3, smooth_b4, smooth_b5, smooth_b6, smooth_b7))
        
        if (naCheck[theOrd[2]]) {
          #If the second cycle did not have enough observations (seg_metrics = NA), only report the first cycle
          out <- c(numCyc, 
                   phen1, seg_max[theOrd[1]], seg_amp[theOrd[1]], seg_int[theOrd[1]], comp1, 
                   gup_numObs[theOrd[1]],   gup_maxgap_frac[theOrd[1]],   gup_maxgap_frac_count_snow[theOrd[1]],  gup_maxgap_frac_filled[theOrd[1]],    gup_rsq[theOrd[1]],   
                   gdown_numObs[theOrd[1]], gdown_maxgap_frac[theOrd[1]], gdown_maxgap_frac_count_snow[theOrd[1]], gdown_maxgap_frac_filled[theOrd[1]], gdown_rsq[theOrd[1]], 
                   matrix(NA,pheno_pars$numLyrsCycle2),
                   numObs, numObs_count_snow, maxGap_annual, maxGap_annual_count_snow)
          
        } else {
          #If the second cycle had enough observations, report both
          out <- c(numCyc, 
                   phen1, seg_max[theOrd[1]], seg_amp[theOrd[1]], seg_int[theOrd[1]], comp1, 
                   gup_numObs[theOrd[1]],   gup_maxgap_frac[theOrd[1]],   gup_maxgap_frac_count_snow[theOrd[1]],   gup_maxgap_frac_filled[theOrd[1]],   gup_rsq[theOrd[1]],   
                   gdown_numObs[theOrd[1]], gdown_maxgap_frac[theOrd[1]], gdown_maxgap_frac_count_snow[theOrd[1]], gdown_maxgap_frac_filled[theOrd[1]],  gdown_rsq[theOrd[1]], 
                   phen2, seg_max[theOrd[2]], seg_amp[theOrd[2]], seg_int[theOrd[2]], comp2, 
                   gup_numObs[theOrd[2]],   gup_maxgap_frac[theOrd[2]],   gup_maxgap_frac_count_snow[theOrd[2]],   gup_maxgap_frac_filled[theOrd[2]],   gup_rsq[theOrd[2]],   
                   gdown_numObs[theOrd[2]], gdown_maxgap_frac[theOrd[2]], gdown_maxgap_frac_count_snow[theOrd[2]], gdown_maxgap_frac_filled[theOrd[2]],  gdown_rsq[theOrd[2]], 
                   numObs, numObs_count_snow, maxGap_annual, maxGap_annual_count_snow)
        }
      }
    },silent=TRUE)  #End of the try block
    if(inherits(log, "try-error")){outAll <- c(outAll,matrix(NA,pheno_pars$numLyrs))
    } else {outAll <- c(outAll,out);remove(out)}
  }
  return(outAll)
}







#---------------------------------------------------------------------
#Function to create QA layer for each file
#Created using layers from "Extended_QA" netCDF file (must be created before this is run)
#Written by Douglas Bolton
#---------------------------------------------------------------------
CreateQA <- function(qaFile, waterMask, yr, params) {
  qa_pars <- params$qa_parameters
  ncFile <- nc_open(qaFile)
  
  
  rNames <- c('gupRsq','gdownRsq',  'gupRsq_2','gdownRsq_2')
  gNames <- c('gupMaxGap','gdownMaxGap', 'gupMaxGap_2','gdownMaxGap_2')
  gsNames <- c('gupMaxGapCountSnow','gdownMaxGapCountSnow',  'gupMaxGapCountSnow','gdownMaxGapCountSnow')
  gfNames <- c('gupMaxGapFilled','gdownMaxGapFilled',  'gupMaxGapFilled_2','gdownMaxGapFilled_2')
  
  qaSegs <- matrix(NA,numPix,length(rNames))
  for (q in 1:length(rNames)) {
    
    Rsq <- ncvar_get(ncFile, rNames[q]) / 10000
    maxGap <- ncvar_get(ncFile, gNames[q])
    maxGapCountSnow <- ncvar_get(ncFile, gsNames[q])
    maxGapFill <- ncvar_get(ncFile, gfNames[q])
 
    nrows <- dim(Rsq)[1]
    ncols <- dim(Rsq)[2] 
    
    #Fill QA classes backwards
    qual <- matrix(6,nrows,ncols)                                                           #Start with all pixels as "No cycles detected"
    qual[(Rsq <= qa_pars$min_r2_mod_quality | maxGapFill > qa_pars$maxGap_mod_quality)]  <- 5           #Phenology detected, but poor quality    
    qual[Rsq > qa_pars$min_r2_mod_quality & maxGapFill <= qa_pars$maxGap_mod_quality]  <- 4             #Moderate quality with snow filled values and fills from alternate years
    qual[Rsq > qa_pars$min_r2_mod_quality & maxGapCountSnow <= qa_pars$maxGap_mod_quality]  <- 3        #Moderate quality with snow filled values
    qual[Rsq > qa_pars$min_r2_mod_quality & maxGap <= qa_pars$maxGap_mod_quality]  <- 2                 #Moderate quality pixel
    qual[Rsq > qa_pars$min_r2_high_quality & maxGap <= qa_pars$maxGap_high_quality]  <- 1               #High quality pixel
  
    #For 2016 - Label first/last two rows/cols as class 6 (Unrealiable HLS data)
    if (yr == 2016) {
      qual[1:2,] <- 9
      qual[,1:2] <- 9
      qual[(nrows-1):nrows,] <- 9
      qual[,(ncols-1):ncols] <- 9
    }
    qual <- matrix(qual,nrows*ncols)
    qual[waterMask == 1] <- 10       #Water pixel, phenology code wasn't run
    
    qaSegs[,q] <- qual
}

  
  numObs <- ncvar_get(ncFile, "numObs")
  numObsCountSnow <- ncvar_get(ncFile, "numObsCountSnow")
  maxGap <- ncvar_get(ncFile, "maxGap")
  maxGapCountSnow <- ncvar_get(ncFile, "maxGapCountSnow")
  
  #quality based on annual layers (8 == poor, 7 = moderate, 6 = high)
  cycle_1 <- matrix(8,dim(numObs)[1],dim(numObs)[2])                         #start all as poor quality
  cycle_1[maxGapCountSnow <= 30 & numObsCountSnow > 23] <- 7                 #moderate quality if average of 2 images per month and no gap > a month (counting snow)
  cycle_1[maxGap <= 30 & numObs > 23] <- 6                                   #high quality if average of 2 images per month and no gap > a month (not counting snow)
  cycle_1 <- matrix(cycle_1, dim(numObs)[1]*dim(numObs)[2])    
  cycle_2 <- cycle_1
  
  cycle_1[qaSegs[,1] == 10] <- 10; cycle_1[qaSegs[,1] == 9] <- 9
  cycle_2[qaSegs[,1] == 10] <- 10; cycle_2[qaSegs[,1] == 9] <- 9
  
  for (i in seq(5,1)) {
    cycle_1[qaSegs[,1] <= i & qaSegs[,2] <= i] <- i
    cycle_2[qaSegs[,3] <= i & qaSegs[,4] <= i] <- i
  }
  
  qaOut <- cbind(cycle_1, cycle_2, qaSegs)
  
  
  return(qaOut)
}











#NetCDF functions
########################################
#Data are to be delivered to the LP-DAAC as netCDF files
#The following functions convert data into netCDF files




#---------------------------------------------------------------------
#Reconstruct output images from chunks (.Rds) in order to make netCDF
#Written by Douglas Bolton
#---------------------------------------------------------------------
readLyrChunks <- function(lyr, yr, numChunks, numPix, tempDir) {
  mat <- matrix(NA,numPix,1)
  for (n in 1:numChunks) {
    fileName <- paste0(tempDir,'outputs/y',yr,'/lyr',lyr,'/c',n,'.Rds')
    matSub <- try(readRDS(fileName),silent=T)
    if (inherits(matSub, 'try-error')) {next} 
    mat[chunkStart[n]:chunkEnd[n]] <- matSub
  }
  return(mat)
}


#---------------------------------------------------------------------
#Get spatial info from HLS file
#Will need in order to write spatial attribute info to netcdf file
#Written by Douglas Bolton
#---------------------------------------------------------------------
getNetCDF_projection_info <- function(baseImage) {

  #Get extent, and then define pixel centers in the x and y direction
  ext = extent(baseImage)
  res = res(baseImage)[1]
  if (res == 60) {res <-  10}              #Error in metadata for S10, where resolutions is written as 60m. Change to 10m.  

  x = seq(ext[1]+res/2,ext[1]+ncol(baseImage)*res, res)
  y = seq(ext[3]+res/2,ext[3]+nrow(baseImage)*res, res)
  
  #Define dimensions for netCDF file
  dimx = ncdim_def(name = 'x', longname = 'x coordinate', units='m', vals = as.double(x))
  dimy = ncdim_def(name = 'y', longname = 'y coordinate', units='m', vals = rev(as.double(y)))
  
  #Get projection in wkt format
  wkt <- showWKT(projection(baseImage))  
  
  #Need to pull the central meridian from the wkt 
  spt <- unlist(strsplit(gsub(']','',wkt),','))
  central_meridian <- as.numeric(spt[which(spt == "PARAMETER[\"central_meridian\"")+1])
  
  #Need "geoTransform" data for netcdf (top left corner coordinates and resolution)
  geoTransform <- paste(ext[1],res,0,ext[4],0,-res)
  
  return(list(dimx=dimx,dimy=dimy, wkt=wkt, central_meridian=central_meridian, geoTransform=geoTransform))
}




#---------------------------------------------------------------------
#Add attributes to each layer in netCDF file
#Attributes and layer names are defined in MuSLI_LSP_V1_Layers.csv file
#Written by Douglas Bolton
#---------------------------------------------------------------------
put_layer_attributes <- function(ncFile, lyr) {
  ncatt_put(ncFile,lyr$short_name,"scale",lyr$scale)
  ncatt_put(ncFile,lyr$short_name,"offset",lyr$offset)
  ncatt_put(ncFile,lyr$short_name,"data_type",lyr$data_type)
  ncatt_put(ncFile,lyr$short_name,"valid_min",lyr$valid_min)
  ncatt_put(ncFile,lyr$short_name,"valid_max",lyr$valid_max)
  ncatt_put(ncFile,lyr$short_name,"grid_mapping","transverse_mercator")  
}



#---------------------------------------------------------------------
#Add spatial info to netCDF file
#Most spatial info for UTM defined in .json file
#However, logitude_of_central_meridian, wkt, and geoTransoform vary per tile/utm
#Written by Douglas Bolton
#---------------------------------------------------------------------
put_prj_attributes <- function(ncFile, prj_info, prj_pars) {
  #First, define variables that do change across tile/utm
  prj_pars$longitude_of_central_meridian <- prj_info$central_meridian
  prj_pars$spatial_ref <- prj_info$wkt
  prj_pars$GeoTransform <- prj_info$geoTransform
  
  #Now loop through all variables in prj_pars and add to "transverse_mercator" variable
  for (i in 1:length(prj_pars)) {ncatt_put(ncFile,"transverse_mercator",names(prj_pars[i]),prj_pars[[i]])}
}


#---------------------------------------------------------------------
#Define global attributes for netCDF file
#All global attributes are defined in .json file
#Will add same attributes to both product layer and extended qa
#Written by Douglas Bolton
#---------------------------------------------------------------------
put_global_attributes <- function(ncFile, global_pars, extendedQA=FALSE) {
  #if it's the extended qa layer, add "- Extended QA" to the name
  if (extendedQA) {global_pars$title <- paste(global_pars$title,'- Extended QA')}  
  #Loop through all global attributes
  for (i in 1:length(global_pars)) {ncatt_put(ncFile,0,names(global_pars[i]),global_pars[[i]])}
}
  
  

#---------------------------------------------------------------------
#Function to create an extended QA layer for each cycle
#These are the layers that are used to build the QA classes and may be useful for refining qa
#Written by Douglas Bolton
#---------------------------------------------------------------------
CreateExtendedQA <- function(yr, qaFile, productTable, baseImage, params) {
  
  numChunks <- params$setup$numChunks                       #Get number of image chunks
  numPix <- dim(baseImage)[1] * dim(baseImage)[2]           #get total number of pixels 
  
  qaLyrs <- productTable[!is.na(productTable$qa_lyr),]      #determine which layers are to be included in extended qa (defined in MuSLI_LSP_V1_Layers.csv)
  qaLyrs <- qaLyrs[order(qaLyrs$qa_lyr),]                   #order the layers
  
  
  prj_info <- getNetCDF_projection_info(baseImage)          #Get required spatial info from an HLS image
  
  results<-vector("list", dim(qaLyrs)[1]+1)                                 #create list where netCDF variables will be defined (spatial variable + qa layers)
  results[[1]] <- ncvar_def("transverse_mercator","",list(),prec="char")    #Define a spatial variable for the file (1st variable)
  
  #Loop through the layers, define a variable for each layer
  for (i in 1:dim(qaLyrs)[1]) {
    lyr <- qaLyrs[i,]          #Pull info for layer
    
    #Create the variable, add to the list. Define the short_name, units, fill_value, long_name, and precision from the lyrs table
    results[[i+1]] <- ncvar_def(lyr$short_name, lyr$units, list(prj_info$dimx,prj_info$dimy), lyr$fill_value, lyr$long_name, prec='short', compression=2)    
  }
  
  
  # Now create the netCDF file with the defined variables
  if (file.exists(qaFile)) {file.remove(qaFile)};
  ncFile <- nc_create(qaFile,results,force_v4=T)
  
  #Loop through layers again, write the image data to the file
  for (i in 1:dim(qaLyrs)[1]) {
    lyr <- qaLyrs[i,]
    
    data <- readLyrChunks(lyr$calc_lyr, yr, numChunks, numPix, params$dirs$tempDir)     #Read all image chunks to recreate full image
    data <- matrix(data, dim(baseImage)[1], dim(baseImage)[2])           #Need to flip image for netcdf. The next lines do this
    
    data[data < -32767 | data > 32767] <- 32767     #Ensure no values are outside the range
    ncvar_put(ncFile,results[[i+1]], data)      #Now put the image into the file
    
    put_layer_attributes(ncFile, lyr)           #Put in attributes for the layer
      
  }
  
  put_prj_attributes(ncFile,prj_info, params$netcdf$transverse_mercator)   #Add spatial info to file
  put_global_attributes(ncFile, params$netcdf$global, extendedQA=TRUE)     #Add global attributes to file
  
  nc_close(ncFile)    #Close the file
  
}
  






#---------------------------------------------------------------------
#Function to generate netCDF file for the product
#The layers to be included in the product are defined in MuSLI_LSP_V1_Layers.csv
#Written by Douglas Bolton
#---------------------------------------------------------------------
CreateProduct <- function(yr,productFile, qaFile, productTable, baseImage, waterMask, params) {
  
  numChunks <- params$setup$numChunks                       #Get number of image chunks
  numPix <- dim(baseImage)[1] * dim(baseImage)[2]           #Get number of image pixels
  startDate <- as.numeric(as.Date(paste0(yr-1,'-12-31')))   #Determine date to subtract from each raster to get DOY
  
  qaLyrs <- CreateQA(qaFile, waterMask, yr, params)         #Create the QA layers (QA classes created from extended QA layer)
    
  lyrs <- productTable[!is.na(productTable$product_lyr),]   #Determine layers to include in product
  lyrs <- lyrs[order(lyrs$product_lyr),]                    #Determine order of layers in product
  
  prj_info <- getNetCDF_projection_info(baseImage)          #Get required spatial info from an HLS image
  
  results<-vector("list", dim(lyrs)[1]+1)                             #Create list where netCDF variables will be defined (spatial variable + qa layers)
  results[[1]] <- ncvar_def("transverse_mercator","",list(),prec="char")    #Define a spatial variable for the file (1st variable)
  
  #Loop through the layers, define a variable for each layer
  for (i in 1:dim(lyrs)[1]) {
    lyr <- lyrs[i,]   #Pull info for layer
    
    #Create the variable, add to the list. Define the short_name, units, fill_value, long_name, and precision from the lyrs table
    results[[i+1]] <- ncvar_def(lyr$short_name, lyr$units, list(prj_info$dimx,prj_info$dimy), lyr$fill_value, lyr$long_name, prec='short', compression=2)   
  }
  
  
  # Now create the netCDF file with the defined variables
  ncFile <- nc_create(productFile,results,force_v4=T)
  
  
  # Assign NA for the EVImax, EVIamp, and EVIare layers where their values exceed 10000, 10000, and 32766, respectively, 
  # and give 6 (i.e., "No cycles detected") for the pixels in all layers
  ngEVI <- matrix(0, dim(baseImage)[1], dim(baseImage)[2])           # A layer to screen all pixels having bad EVI values
  
  for (i in 1:length(lyrs$short_name)) {
    lyr <- lyrs[i,]
    
    if (lyr$short_name == 'EVImax' | lyr$short_name == 'EVIamp'| lyr$short_name == 'EVImax_2' | lyr$short_name == 'EVIamp_2') {
      data <- readLyrChunks(lyr$calc_lyr, yr, numChunks, numPix, params$dirs$tempDir)
      data <- matrix(data, dim(baseImage)[1], dim(baseImage)[2])           #Need to flip image for netcdf. The next lines do this
      
      ngEVI[data > 10000] <- 1      #Set as 1 for the pixels where EVImax or EVIamp exceed 10000
      
    }else if(lyr$short_name == 'EVIarea' | lyr$short_name == 'EVIarea_2'){
      data <- readLyrChunks(lyr$calc_lyr, yr, numChunks, numPix, params$dirs$tempDir)
      data <- matrix(data, dim(baseImage)[1], dim(baseImage)[2])           #Need to flip image for netcdf. The next lines do this
      
      ngEVI[data > 32766] <- 1      #Set as 1 for the pixels where EVIare exceed 32766
    }
  }
  
  #Loop through layers again, write the image data to the file
  for (i in 1:length(lyrs$short_name)) {
    lyr <- lyrs[i,]
    
    
    if (lyr$short_name == 'overallQA') {data <- qaLyrs[,1]
    } else if  (lyr$short_name == 'overallQA_2') {data <- qaLyrs[,2]
    } else if  (lyr$short_name == 'gupQA') {data <- qaLyrs[,3]
    } else if  (lyr$short_name == 'gdownQA') {data <- qaLyrs[,4]
    } else if  (lyr$short_name == 'gupQA_2') {data <- qaLyrs[,5]
    } else if  (lyr$short_name == 'gdownQA_2') {data <- qaLyrs[,6]
    } else {   data <- readLyrChunks(lyr$calc_lyr, yr, numChunks, numPix, params$dirs$tempDir)} #Read all image chunks to recreate full image
               
    data <- matrix(data, dim(baseImage)[1], dim(baseImage)[2])           #Need to flip image for netcdf. The next lines do this

    if (lyr$units == 'Day of year') {data <- data - startDate}    #If units are Day of year, convert from date to day of year
    
    
    if (lyr$short_name == 'overallQA' |lyr$short_name == 'overallQA_2' |lyr$short_name == 'gupQA' |lyr$short_name == 'gdownQA' |lyr$short_name == 'gupQA_2' |lyr$short_name == 'gdownQA_2'){
      data[ngEVI == 1] <- 6 #Give 6 (i.e., "No cycle detected") for the pixels having bad EVI values 
    }else{
      data[ngEVI == 1] <- 32767  #Give NA for the pixels having bad EVI values 
    }
    
    
    data[data < -32767 | data > 32767] <- 32767      #Ensure no values are outside the range
    ncvar_put(ncFile,results[[i+1]], data)           #Now put the image into the file
    
    put_layer_attributes(ncFile, lyr)                #Put in attributes for the layer
    
  }
  
  put_prj_attributes(ncFile,prj_info, params$netcdf$transverse_mercator)   #Put spatial info 
  put_global_attributes(ncFile, params$netcdf$global)                      #Put global attributes
  
  nc_close(ncFile)    #Close the file
  
}












