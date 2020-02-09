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




#Inputs
############

#Tile 
tile <- '10SGF'
codeVersion <- 'V1'      #run V0 or V1?

numCores <- 8

imgYrs <- 2016:2019
phenYrs <- 2016:2019

#What data folder
chunkBase <- '/projectnb/modislc/projects/landsat_sentinel/MSLSP_HLS30/'
imgBase <-   '/projectnb/modislc/projects/landsat_sentinel/v1_4/HLS30/'


#Define and create output directory
outFolder <- '/projectnb/modislc/projects/landsat_sentinel/MSLSP_assessment/plots/'


#Define path to function
functions_V0 <- '/usr2/postdoc/dbolt/CodeFolder/git_repos/MuSLI_LSP/AWS/r/MuSLI_LSP_Functions.r'
functions_V1 <- '/usr2/postdoc/dbolt/CodeFolder/git_repos/MSLSP/MSLSP_Functions.r'
functions_diagnostics <- "/usr2/postdoc/dbolt/CodeFolder/git_repos/MSLSP/supporting/MSLSP_Diagnostic_Functions_V1_0.r"

#Name of shapefile. Must be in same projection as the tile. Must have "id" column
shpName <- paste0('/projectnb/modislc/projects/landsat_sentinel/MSLSP_assessment/shps/',tile,'_pts.shp')


#Json file where phenology parameters are defined. Will just use the phenology paramaters from this 
jsonFile <- "/projectnb/modislc/projects/landsat_sentinel/MSLSP_assessment/MSLSP_Parameters.json"


showObservations <- T
showSpline <- T
showPhenDates <- F
showFilledData <- F
showDespiked <- F
showSnow <- F
yrsToPlot <- 2018:2019

#Get default parameters
params <- fromJSON(file=jsonFile)


#Alter the jsonFile?
params$phenology_parameters$splineSpar=0.2      #Alter the smoothing parameter





registerDoMC(cores=numCores)

#Load functions 
if (codeVersion=='V0') {source(file=functions_V0)}
if (codeVersion=='V1') {source(file=functions_V1)}
source(file=functions_diagnostics)






#Where is the data?
imgDir <- paste0(imgBase,tile,'/images/')
chunkDir <- paste0(chunkBase,tile,'/imageChunks/')




theTable <- Extract_Timeseries(tile, imgDir, chunkDir, imgYrs, phenYrs, numCores, params, shpName=shpName,codeVersion=codeVersion) 



#Create output directory if it doesn't exist
outDir <- paste0(outFolder,tile,'_',codeVersion,'/')
if (!dir.exists(outDir)) {dir.create(outDir)}



idNames <- names(theTable)

for (i in 1:length(theTable)) {
  
    pixTab <- theTable[[i]]
    idName <- idNames[i]
    
    for (y in 1:length(yrsToPlot)) {
      yr <- yrsToPlot[y]
      
      subTab <- pixTab[[paste0('y',yr)]]
      
  
      output_name = paste0(outDir,idName,'_',yr,'_',codeVersion,'.tif')
      png(output_name,res = 600,width = 6,height = 3.5,units = "in",pointsize = 8)
      par(fig=c(0.02,.98,.02,.98),mai=c(0.33,0.33,.3,.1),mgp=c(1.5,0.2,0))
          
      xlims <- as.Date(c(paste0(yrsToPlot[y]-1,'-11-1'),paste0(yrsToPlot[y]+1,'-02-1')))  #Always use the first of a month! (for labelling purposes)
      ylims <-  c(0.05,0.9)
      
      plot(subTab$filled_vi ~ as.Date(subTab$filled_dates),
               type='n',
               xlim=xlims,ylim=ylims,tck=0.025,
               xlab='Date',ylab='EVI2',xaxt='n',
               main=paste(idName,yr,codeVersion))
          
  
      if (showPhenDates) {abline(v=subTab$phen,lty='dashed',col='gray',lwd=2)}  
          
      if (showSpline) {lines(subTab$smoothed_vi ~  subTab$smoothed_dates,col='#31a354',lwd=2)}
          
      #Plot only the observations from other years. This will also plot snow from current year, but we'll cover those points up later
      if (showFilledData) {points(subTab$filled_vi[subTab$filled_weigth < 1] ~ subTab$filled_dates[subTab$filled_weigth < 1],pch=16,cex=0.6)}
        
      #Plot the original time-series before despiking  
      if (showDespiked) {points(pixTab$original_VI~pixTab$dates)}
       
      if (showObservations) {points(pixTab$filled_VI~pixTab$dates,pch=21,bg='#2c7fb8')}
      if (showSnow) {points(pixTab$filled_VI[pixTab$snowPix]~pixTab$dates[pixTab$snowPix],pch=21,bg='black')}
       

          
        axis.Date(1,at=seq(xlims[1], xlims[2],by="months"),
                    labels=format(seq(xlims[1], xlims[2],by="months"),'%b'),
                    tck=0.025)
          
        dev.off()
    }
    
    
    #Do a full time-series plot. Need to fit new spline
    output_name = paste0(outDir,idName,'_Full_Series_',codeVersion,'.tif')
    png(output_name,res = 600,width = 6,height = 3.5,units = "in",pointsize = 8)
    par(fig=c(0.02,.98,.02,.98),mai=c(0.33,0.33,.3,.1),mgp=c(1.5,0.2,0))
    
    xlims <- as.Date(c(paste0(min(yrsToPlot)-2,'-11-1'),paste0(max(yrsToPlot)+1,'-02-1')))  #Always use the first of a month! (for labelling purposes)
    ylims <-  c(0.05,0.9)
    
    plot(pixTab$original_VI~pixTab$dates,
         type='n',
         xlim=xlims,ylim=ylims,tck=0.025,
         xlab='Date',ylab='EVI2',xaxt='n',
         main=paste(idName,codeVersion))
    
    #Get index of pixels with good values
    ind <- !is.na(pixTab$original_VI)  
    # smooth with a spline to get continuous daily series
    spl <- smooth.spline(pixTab$dates[ind], pixTab$original_VI[ind], spar=0.1)
    # weighted version
    xSmooth <- predict(spl, as.numeric(pixTab$dates))$y
    
    if (showSpline) {lines(xSmooth~pixTab$dates,col='#31a354',lwd=2)}

    #Plot the original time-series before despiking  
    if (showDespiked) {points(pixTab$original_VI~pixTab$dates)}
    
    if (showObservations) {points(pixTab$filled_VI~pixTab$dates,pch=21,bg='#2c7fb8')}
    if (showSnow) {points(pixTab$filled_VI[pixTab$snowPix]~pixTab$dates[pixTab$snowPix],pch=21,bg='black')}
    
    
    axis.Date(1,at=seq(xlims[1], xlims[2],by="months"),
              labels=format(seq(xlims[1], xlims[2],by="months"),'%b'),
              tck=0.025)
    dev.off()
}
          
 




      
      
      
      

