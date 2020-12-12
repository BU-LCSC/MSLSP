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
tile <- '17SQD'
codeVersion <- 'V1'      #run V0 or V1?

numCores <- 8

imgYrs <- 2016:2019
phenYrs <- 2016:2019

#What data folder
chunkBase <- '/projectnb/modislc/projects/landsat_sentinel/MSLSP_HLS30/'
# chunkBase <- '/projectnb/modislc/users/mkmoon/MuSLI/V1_0/From_AWS/t2/'
imgBase <-   '/projectnb/modislc/projects/landsat_sentinel/v1_4/HLS30/'


#Define and create output directory
# outFolder <- '/projectnb/modislc/projects/landsat_sentinel/MSLSP_assessment/plots/'
outFolder <- '/projectnb/modislc/users/mkmoon/MuSLI/V1_0/figures/edge/'


#Define path to function
functions_V0 <- '/usr3/graduate/mkmoon/GitHub/musli/MuSLI_LSP_Functions.r'
functions_V1 <- '/usr3/graduate/mkmoon/GitHub/MSLSP/MSLSP_Functions.r'
functions_diagnostics <- "/usr3/graduate/mkmoon/GitHub/MSLSP/Development/MSLSP_Diagnostic_Functions_V1_0.r"

#Name of shapefile. Must be in same projection as the tile. Must have "id" column
# shpName <- paste0('/projectnb/modislc/projects/landsat_sentinel/MSLSP_assessment/shps/',tile,'_pts.shp')
shpName <- paste0('/projectnb/modislc/users/mkmoon/MuSLI/geo_data/shp/',tile,'_pts_edge.shp')
# shpName <- paste0('/projectnb/modislc/users/mkmoon/Planet/shp/rp_mf_1000.shp')

#Json file where phenology parameters are defined. Will just use the phenology paramaters from this 
jsonFile <- "/usr3/graduate/mkmoon/GitHub/MSLSP/MSLSP_Parameters.json"


showObservations <- T
showSpline <- T
showPhenDates <- T
showFilledData <- T
showDespiked <- T
showSnow <- T
yrsToPlot <- 2016:2019
PlanetS <- T

#Get default parameters
params <- fromJSON(file=jsonFile)


#Alter the jsonFile?
params$phenology_parameters$splineSpar=0.55      #Alter the smoothing parameter
# params$phenology_parameters$min_seg_amplitude=0.1    #Set minimum seg amplitude
#params$phenology_parameters$vegetation_index='ndvi_re1'  #Alter the vegetation index (Find list of vegetation indices in MSLSP_Functions.r file: CalcIndex function)




registerDoMC(cores=numCores)

#Load functions 
if (codeVersion=='V0') {source(file=functions_V0)}
if (codeVersion=='V1') {source(file=functions_V1)}
source(file=functions_diagnostics)






#Where is the data?
imgDir <- paste0(imgBase,tile,'/images/')
chunkDir <- paste0(chunkBase,tile,'/imageChunks/')




theTable <- Extract_Timeseries(tile, imgDir, chunkDir, imgYrs, phenYrs, numCores, params, shpName=shpName,codeVersion=codeVersion) 

# setwd('/projectnb/modislc/users/mkmoon/Planet/data/')
# save(theTable,file='rp_mf_hls_1000.rda')


# ################################################################
# #Tile 
# codeVersion <- 'V1'      #run V0 or V1?
# 
# numCores <- 8
# 
# imgYrs <- 2016:2018
# phenYrs <- 2016:2018
# 
# #What data folder
# chunkBase <- '/projectnb/modislc/projects/landsat_sentinel/MSLSP_HLS30/'
# # chunkBase <- '/projectnb/modislc/users/mkmoon/MuSLI/V1_0/From_AWS/t2/'
# imgBase <-   '/projectnb/modislc/projects/landsat_sentinel/v1_4/HLS30/'
# 
# 
# #Define and create output directory
# # outFolder <- '/projectnb/modislc/projects/landsat_sentinel/MSLSP_assessment/plots/'
# outFolder <- '/projectnb/modislc/users/mkmoon/MuSLI/V1_0/figures/edge/'
# 
# 
# #Define path to function
# functions_V0 <- '/usr3/graduate/mkmoon/GitHub/musli/MuSLI_LSP_Functions.r'
# functions_V1 <- '/usr3/graduate/mkmoon/GitHub/MSLSP/MSLSP_Functions.r'
# functions_diagnostics <- "/usr3/graduate/mkmoon/GitHub/MSLSP/Development/MSLSP_Diagnostic_Functions_V1_0.r"
# 
# #Name of shapefile. Must be in same projection as the tile. Must have "id" column
# # shpName <- paste0('/projectnb/modislc/projects/landsat_sentinel/MSLSP_assessment/shps/',tile,'_pts.shp')
# shpName <- paste0('/projectnb/modislc/users/mkmoon/MuSLI/geo_data/shp/',tile,'_pts_edge.shp')
# # shpName <- paste0('/projectnb/modislc/users/mkmoon/Planet/shp/rp_mf_1000.shp')
# 
# #Json file where phenology parameters are defined. Will just use the phenology paramaters from this 
# jsonFile <- "/usr3/graduate/mkmoon/GitHub/MSLSP/MSLSP_Parameters.json"
# 
# 
# showObservations <- T
# showSpline <- T
# showPhenDates <- T
# showFilledData <- T
# showDespiked <- T
# showSnow <- T
# yrsToPlot <- 2016:2018
# 
# #Get default parameters
# params <- fromJSON(file=jsonFile)
# 
# 
# #Alter the jsonFile?
# params$phenology_parameters$splineSpar=0.55      #Alter the smoothing parameter
# # params$phenology_parameters$min_seg_amplitude=0.1    #Set minimum seg amplitude
# #params$phenology_parameters$vegetation_index='ndvi_re1'  #Alter the vegetation index (Find list of vegetation indices in MSLSP_Functions.r file: CalcIndex function)
# 
# 
# 
# 
# registerDoMC(cores=numCores)
# 
# #Load functions 
# source(file=functions_V1)
# source(file=functions_diagnostics)
# 
# theTable1 <- Extract_Timeseries(tile, imgDir, chunkDir, imgYrs, phenYrs, numCores, params, shpName=shpName,codeVersion=codeVersion) 



#Create output directory if it doesn't exist
outDir <- paste0(outFolder,tile,'_V1_withPS/')
if (!dir.exists(outDir)) {dir.create(outDir)}



idNames <- names(theTable)
# idNames1 <- names(theTable1)

if (PlanetS) {
  source('/usr3/graduate/mkmoon/GitHub/PlanetLSP/PLSP_Functions.R')
  load('/projectnb/modislc/users/mkmoon/Planet/data/planet_eviStack.rda')
  shpPoints <- readOGR(shpName)
  
}


for (i in 1:length(theTable)) {
  
    pixTab <- theTable[[i]]
    idName <- idNames[i]
    
    # pixTab1 <- theTable1[[i]]
    # idName1 <- idNames1[i]
    
    # for (y in 1:length(yrsToPlot)) {
    for (y in 4) {
      yr <- yrsToPlot[y]
      
      subTab <- pixTab[[paste0('y',yr)]]
      # subTab1 <- pixTab1[[paste0('y',yr)]]
      
  
      output_name = paste0(outDir,idName,'_',yr,'_',codeVersion,'.png')
      png(output_name,res = 600,width = 6,height = 3,units = "in",pointsize = 8)
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
       
      legend('topright',c('Obs.','Filled','Despiked','Snow'),pch=c(21,16,1,21),pt.bg=c('#2c7fb8',NA,NA,'black'),
             pt.cex=c(1.2,0.8,1.2,1.2),bty='n',cex=1.1)
          
        axis.Date(1,at=seq(xlims[1], xlims[2],by="months"),
                    labels=format(seq(xlims[1], xlims[2],by="months"),'%b'),
                    tck=0.025)
        
      if (PlanetS) {
        
        datPTS1 <- extractTS(shpPoints,leviStack[[1]],limgBase[[1]],id=as.numeric(substr(idName,7,9)),rad=15)
        d2019 <- c(min(which(substr(lDates[[1]],1,7)=='2018-11')):max(which(substr(lDates[[1]],1,7)=='2020-02')))
        
        for(i in 1:dim(datPTS1)[1]){
          lines(lDates[[1]][d2019],datPTS1[i,d2019],pch=19,cex=0.3,col=rgb(0,0,0,0.2))
          
        }
      }  
      # Add PlanetScope
        
        # #
        # par(fig=c(0.02,.98,.02,.5),mai=c(0.33,0.33,.3,.1),mgp=c(1.5,0.2,0),new=T)
        # 
        # xlims <- as.Date(c(paste0(yrsToPlot[y]-1,'-11-1'),paste0(yrsToPlot[y]+1,'-02-1')))  #Always use the first of a month! (for labelling purposes)
        # ylims <-  c(0.05,0.9)
        # 
        # plot(subTab1$filled_vi ~ as.Date(subTab1$filled_dates),
        #      type='n',
        #      xlim=xlims,ylim=ylims,tck=0.025,
        #      xlab='Date',ylab='EVI2',xaxt='n',
        #      main=paste(idName1,yr,'V1'))
        # 
        # 
        # if (showPhenDates) {abline(v=subTab1$phen,lty='dashed',col='gray',lwd=2)}  
        # 
        # if (showSpline) {lines(subTab1$smoothed_vi ~  subTab1$smoothed_dates,col='#31a354',lwd=2)}
        # 
        # #Plot only the observations from other years. This will also plot snow from current year, but we'll cover those points up later
        # if (showFilledData) {points(subTab1$filled_vi[subTab1$filled_weigth < 1] ~ subTab1$filled_dates[subTab1$filled_weigth < 1],pch=16,cex=0.6)}
        # 
        # #Plot the original time-series before despiking  
        # if (showDespiked) {points(pixTab$original_VI~pixTab$dates)}
        # 
        # if (showObservations) {points(pixTab$filled_VI~pixTab$dates,pch=21,bg='#2c7fb8')}
        # if (showSnow) {points(pixTab$filled_VI[pixTab$snowPix]~pixTab$dates[pixTab$snowPix],pch=21,bg='black')}
        # 
        # legend('topright',c('Obs.','Filled','Despiked','Snow'),pch=c(21,16,1,21),pt.bg=c('#2c7fb8',NA,NA,'black'),
        #        pt.cex=c(1.2,0.8,1.2,1.2),bty='n',cex=1.1) 
        # 
        # axis.Date(1,at=seq(xlims[1], xlims[2],by="months"),
        #           labels=format(seq(xlims[1], xlims[2],by="months"),'%b'),
        #           tck=0.025)
          
        dev.off()
    }
    
    
    # #Do a full time-series plot. Need to fit new spline
    # output_name = paste0(outDir,idName,'_Full_Series_',codeVersion,'.png')
    # png(output_name,res = 600,width = 6,height = 3.5,units = "in",pointsize = 8)
    # par(fig=c(0.02,.98,.02,.98),mai=c(0.33,0.33,.3,.1),mgp=c(1.5,0.2,0))
    # 
    # ylims <-  c(0.05,0.9)
    # 
    # plot(pixTab$original_VI~pixTab$dates,
    #      type='n',
    #      ylim=ylims,tck=0.025,
    #      xlab='Date',ylab='EVI2',xaxt='n',
    #      main=paste(idName,codeVersion))
    # 
    # #Get index of pixels with good values
    # ind <- !is.na(pixTab$original_VI)  
    # # smooth with a spline to get continuous daily series
    # spl <- smooth.spline(pixTab$dates[ind], pixTab$original_VI[ind], spar=0.1)
    # # weighted version
    # xSmooth <- predict(spl, as.numeric(pixTab$dates))$y
    # 
    # if (showSpline) {lines(xSmooth~pixTab$dates,col='#31a354',lwd=2)}
    # 
    # #Plot the original time-series before despiking  
    # if (showDespiked) {points(pixTab$original_VI~pixTab$dates)}
    # 
    # if (showObservations) {points(pixTab$filled_VI~pixTab$dates,pch=21,bg='#2c7fb8')}
    # if (showSnow) {points(pixTab$filled_VI[pixTab$snowPix]~pixTab$dates[pixTab$snowPix],pch=21,bg='black')}
    # 
    # #axis.Date(1,at=seq(xlims[1], xlims[2],by="months"),
    # #          labels=format(seq(xlims[1], xlims[2],by="months"),'%b'),
    # #          tck=0.025)
    # dev.off()
}
          
 




      
      
      
      

