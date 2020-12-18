#Load required libraries
library(raster)
library(rgdal)
library(gdalUtils)
library(rgeos)
library(ncdf4)
library(RColorBrewer)
library(viridis)
library(lattice)

args <- commandArgs()
print(args)

tt <- as.numeric(args[3])
# tt <- 1

##
# tiles <- c('17SQD','18TYN','19TEL','16TCK','13TEF','10TGP')
path <- '/projectnb/modislc/users/mkmoon/MuSLI/V1_0/From_AWS/t6'
# tiles <- substr(list.dirs(path=path,full.names=T)[2:25],53,57)
tiles <- substr(list.dirs(path=path,full.names=T)[2:4],56,60)

year <- 2016:2019

spec <- viridis(11)
mycolRamp = colorRampPalette(c('White',spec))

for(yy in 1:length(year)){
  ##
  # pathSCC <- paste0('/projectnb/modislc/projects/landsat_sentinel/MSLSP_HLS30/',tiles[tt],'/phenoMetrics')
  pathSCC <- paste0('/projectnb/modislc/users/mkmoon/MuSLI/V0_11/',tiles[tt])             # V0 code with SplinePar 0.55
  # pathAWS <- paste0('/projectnb/modislc/users/mkmoon/MuSLI/V1_0/From_AWS/t6/',tiles[tt])  # V1 code with SplinePar 0.55
  pathAWS1 <- paste0('/projectnb/modislc/users/mkmoon/MuSLI/V1_0/From_AWS/t2/',tiles[tt]) # V1 code with SplinePar 0.4
  pathAWS2 <- paste0('/projectnb/modislc/users/mkmoon/MuSLI/V1_0/From_AWS/t5/',tiles[tt]) # Same code with t2
  pathAWS <- paste0('/projectnb/modislc/projects/landsat_sentinel/MSLSP_HLS30/',tiles[tt],'/phenoMetrics') # Same code with t2
  
    
  #
  sstr <- paste0('*',tiles[tt],'*',year[yy],'.nc')
  
  fileSCC <- list.files(pathSCC,pattern=glob2rx(sstr),recursive=F,full.names=T)
  fileAWS <- list.files(pathAWS,pattern=glob2rx(sstr),recursive=F,full.names=T)
  fileAWS1 <- list.files(pathAWS1,pattern=glob2rx(sstr),recursive=F,full.names=T)
  fileAWS2 <- list.files(pathAWS2,pattern=glob2rx(sstr),recursive=F,full.names=T)
  
  nc <- nc_open(fileSCC)
  var <- names(nc[['var']])
  
  setwd('/projectnb/modislc/users/mkmoon/MuSLI/V1_0/figures/v0_v1_t7/')
  for(pp in 3:11){
    rastSCC <- raster(fileSCC,varname=var[pp])
    rastAWS <- raster(fileAWS,varname=var[pp])  
    rastAWS1 <- raster(fileAWS1,varname=var[pp])  
    rastAWS2 <- raster(fileAWS2,varname=var[pp])  
    
    # #
    # rastEVIamp <- raster(fileAWS,varname='EVIamp')
    # rastAWS[rastEVIamp<1000] <- NA 
    # #
    
    png(filename=paste('diff_MSLSP_',tiles[tt],'_',year[yy],'_',sprintf('%02d',pp),'.png',sep=''),
        width=12,height=7,unit='in',res=150)
    
    par(mfrow=c(2,3),oma=c(1,1,1,1),mar=c(4,4,1,1),mgp=c(2.5,1,0))
    plot(rastSCC,zlim=c(quantile(values(rastAWS),0,na.rm=T),quantile(values(rastAWS),1,na.rm=T)),
         main=paste(var[pp],'_V0_',
                            round(sum(!is.na(values(rastSCC)))/sum(!is.na(values(rastAWS)))*100,0),'%',
                            sep=''),
         cex.main=1.5)
    plot(rastAWS,zlim=c(quantile(values(rastAWS),0,na.rm=T),quantile(values(rastAWS),1,na.rm=T)),
         main=paste(var[pp],'_V1_',sum(!is.na(values(rastAWS))),sep=''),cex.main=1.5)
    # plot(rastSCC,rastAWS,
    #      xlab=paste(var[3],'_SCC',sep=''),
    #      ylab=paste(var[3],'_AWS',sep=''))
    if(pp<9){
      limx <- c(-100,500)
      limy <- c(-100,500)
    }else if(pp>8 & pp<11){
      limx <- c(0,12000)
      limy <- c(0,12000)
    }else{
      limx <- c(0,25000)
      limy <- c(0,25000)
    }
    smoothScatter(values(rastSCC),values(rastAWS),
                  xlim=limx,ylim=limy,
                  colramp=mycolRamp,nbin=1000,nrpoints=0,
                  transformation=function(x)x^.4,axe=F,
                  xlab=paste(var[pp],'_V0',sep=''),
                  ylab=paste(var[pp],'_V1',sep=''),
                  cex.lab=1.2)
    if(pp>8){
      axis(1,at=seq(0,30000,5000),cex.axis=1.2)
      axis(2,at=seq(0,30000,5000),cex.axis=1.2)  
    }else{
      axis(1,at=seq(-100,500,100),cex.axis=1.2)
      axis(2,at=seq(-100,500,100),cex.axis=1.2)  
    }
    
    abline(0,1,lty=5)
    if(pp==11){
      hist(values(rastAWS),breaks=seq(-20000,50000,10),xlim=c(0,15000),col=rgb(1,1,1),border=rgb(1,0,0,0.3))  
      hist(values(rastSCC),breaks=seq(-20000,50000,10),xlim=c(0,15000),col=rgb(1,1,1),border=rgb(0,0,1,0.3),add=T)   
      legend('topright',c('V0','V1'),pch=15,col=c('blue','red'),pt.cex=1.8,bty='n',cex=1.3)
      
      # htPP  <- hist((values(rastAWS)-values(rastSCC)),breaks=seq(-100000,100000,10),plot=F)
      # htPP$counts <- htPP$counts/(sum(!is.na(values(rastAWS)-values(rastSCC))))*100
      # plot(htPP,xlim=c(-1000,1000),col=rgb(1,1,1),border=rgb(0,0,0),
      #      main='',xlab='V1 - V0',
      #      ylab='Frequency (%)',
      #      cex.lab=1.5,cex.axis=1.5)
      
      htPP  <- hist(abs(values(rastAWS)-values(rastSCC)),breaks=seq(-100000,100000,10),plot=F)
      htPP$counts <- htPP$counts/(sum(!is.na(values(rastAWS)-values(rastSCC))))*100
      plot(htPP$breaks[2:length(htPP$breaks)],cumsum(htPP$counts),
           xlim=c(0,2000),ylim=c(0,100),
           type='o',
           xlab='V1 - V0  (Abs.)',
           ylab='Cumulative (%)',
           cex.lab=1.5,cex.axis=1.5,col='red',pch=19)
      
      htPP  <- hist(abs(values(rastAWS1)-values(rastSCC)),breaks=seq(-100000,100000,10),plot=F)
      htPP$counts <- htPP$counts/(sum(!is.na(values(rastAWS1)-values(rastSCC))))*100
      points(htPP$breaks[2:length(htPP$breaks)],cumsum(htPP$counts),
             xlim=c(0,2000),ylim=c(0,100),
             type='o',
             xlab='V1 - V0  (Abs.)',
             ylab='Cumulative (%)',
             cex.lab=1.5,cex.axis=1.5,col='blue',pch=19)
      
      htPP  <- hist(abs(values(rastAWS2)-values(rastAWS1)),breaks=seq(-100000,100000,10),plot=F)
      htPP$counts <- htPP$counts/(sum(!is.na(values(rastAWS2)-values(rastAWS1))))*100
      points(htPP$breaks[2:length(htPP$breaks)],cumsum(htPP$counts),
             xlim=c(0,2000),ylim=c(0,100),
             type='o',
             xlab='V1 - V0  (Abs.)',
             ylab='Cumulative (%)',
             cex.lab=1.5,cex.axis=1.5,col='darkgreen',pch=19)
      
      legend('bottomright',c('SplinePar: 0.55','SplinePar: 0.40','Same code'),pch=19,cex=1.5,pt.cex=1.5,col=c('red','blue','darkgreen'),bty='n')
      
      
      htPP  <- hist((abs((values(rastAWS)-values(rastSCC))/values(rastSCC)*100)),breaks=seq(-100000,100000,1),plot=F)
      htPP$counts <- htPP$counts/(sum(!is.na(values(rastAWS)-values(rastSCC))))*100
      plot(htPP$breaks[2:length(htPP$breaks)],cumsum(htPP$counts),
           xlim=c(0,40),ylim=c(0,100),
           type='o',
           xlab='(V1-V0)/V0x100  (%)',
           ylab='Cumulative (%)',
           cex.lab=1.5,cex.axis=1.5,col='red',pch=19)
      
      abline(v=10,lty=5)
      
      htPP  <- hist((abs((values(rastAWS1)-values(rastSCC))/values(rastSCC)*100)),breaks=seq(-100000,100000,1),plot=F)
      htPP$counts <- htPP$counts/(sum(!is.na(values(rastAWS1)-values(rastSCC))))*100
      points(htPP$breaks[2:length(htPP$breaks)],cumsum(htPP$counts),
             xlim=c(0,40),ylim=c(0,100),
             type='o',
             xlab='(V1-V0)/V0x100  (%)',
             ylab='Cumulative (%)',
             cex.lab=1.5,cex.axis=1.5,col='blue',pch=19)
      
      htPP  <- hist((abs((values(rastAWS2)-values(rastAWS1))/values(rastAWS1)*100)),breaks=seq(-100000,100000,1),plot=F)
      htPP$counts <- htPP$counts/(sum(!is.na(values(rastAWS2)-values(rastAWS1))))*100
      points(htPP$breaks[2:length(htPP$breaks)],cumsum(htPP$counts),
             xlim=c(0,40),ylim=c(0,100),
             type='o',
             xlab='(V1-V0)/V0x100  (%)',
             ylab='Cumulative (%)',
             cex.lab=1.5,cex.axis=1.5,col='darkgreen',pch=19)
      
      legend('bottomright',c('SplinePar: 0.55','SplinePar: 0.40','Same code'),pch=19,cex=1.5,pt.cex=1.5,col=c('red','blue','darkgreen'),bty='n')
      
    }else if(pp>8 & pp<11){
      hist(values(rastAWS),breaks=seq(-20000,50000,10),xlim=c(0,15000),col=rgb(1,1,1),border=rgb(1,0,0,0.3))  
      hist(values(rastSCC),breaks=seq(-20000,50000,10),xlim=c(0,15000),col=rgb(1,1,1),border=rgb(0,0,1,0.3),add=T)   
      legend('topright',c('V0','V1'),pch=15,col=c('blue','red'),pt.cex=1.8,bty='n',cex=1.3)
    
      # htPP  <- hist((values(rastAWS)-values(rastSCC)),breaks=seq(-100000,100000,10),plot=F)
      # htPP$counts <- htPP$counts/(sum(!is.na(values(rastAWS)-values(rastSCC))))*100
      # plot(htPP,xlim=c(-1000,1000),col=rgb(1,1,1),border=rgb(0,0,0),
      #      main='',xlab='V1 - V0',
      #      ylab='Frequency (%)',
      #      cex.lab=1.5,cex.axis=1.5)
      
      htPP  <- hist(abs(values(rastAWS)-values(rastSCC)),breaks=seq(-1000000,1000000,10),plot=F)
      htPP$counts <- htPP$counts/(sum(!is.na(values(rastAWS)-values(rastSCC))))*100
      plot(htPP$breaks[2:length(htPP$breaks)],cumsum(htPP$counts),
           xlim=c(0,1000),ylim=c(0,100),
           type='o',
           xlab='V1 - V0  (Abs.)',
           ylab='Cumulative (%)',
           cex.lab=1.5,cex.axis=1.5,col='red',pch=19)
      
      htPP  <- hist(abs(values(rastAWS1)-values(rastSCC)),breaks=seq(-1000000,1000000,10),plot=F)
      htPP$counts <- htPP$counts/(sum(!is.na(values(rastAWS1)-values(rastSCC))))*100
      points(htPP$breaks[2:length(htPP$breaks)],cumsum(htPP$counts),
           xlim=c(0,1000),ylim=c(0,100),
           type='o',
           xlab='V1 - V0  (Abs.)',
           ylab='Cumulative (%)',
           cex.lab=1.5,cex.axis=1.5,col='blue',pch=19)
      
      htPP  <- hist(abs(values(rastAWS2)-values(rastAWS1)),breaks=seq(-100000,100000,10),plot=F)
      htPP$counts <- htPP$counts/(sum(!is.na(values(rastAWS2)-values(rastAWS1))))*100
      points(htPP$breaks[2:length(htPP$breaks)],cumsum(htPP$counts),
             xlim=c(0,1000),ylim=c(0,100),
             type='o',
             xlab='V1 - V0  (Abs.)',
             ylab='Cumulative (%)',
             cex.lab=1.5,cex.axis=1.5,col='darkgreen',pch=19)
      
      legend('bottomright',c('SplinePar: 0.55','SplinePar: 0.40','Same code'),pch=19,cex=1.5,pt.cex=1.5,col=c('red','blue','darkgreen'),bty='n')
      
      htPP  <- hist((abs((values(rastAWS)-values(rastSCC))/values(rastSCC)*100)),breaks=seq(-100000,1000000,1),plot=F)
      htPP$counts <- htPP$counts/(sum(!is.na(values(rastAWS)-values(rastSCC))))*100
      plot(htPP$breaks[2:length(htPP$breaks)],cumsum(htPP$counts),
           xlim=c(0,40),ylim=c(0,100),
           type='o',
           xlab='(V1-V0)/V0x100  (%)',
           ylab='Cumulative (%)',
           cex.lab=1.5,cex.axis=1.5,col='red',pch=19)
      
      abline(v=10,lty=5)
      
      htPP  <- hist((abs((values(rastAWS1)-values(rastSCC))/values(rastSCC)*100)),breaks=seq(-1000000,1000000,1),plot=F)
      htPP$counts <- htPP$counts/(sum(!is.na(values(rastAWS1)-values(rastSCC))))*100
      points(htPP$breaks[2:length(htPP$breaks)],cumsum(htPP$counts),
             xlim=c(0,40),ylim=c(0,100),
             type='o',
             xlab='(V1-V0)/V0x100  (%)',
             ylab='Cumulative (%)',
             cex.lab=1.5,cex.axis=1.5,col='blue',pch=19)
      
      htPP  <- hist((abs((values(rastAWS2)-values(rastAWS1))/values(rastAWS1)*100)),breaks=seq(-1000000,1000000,1),plot=F)
      htPP$counts <- htPP$counts/(sum(!is.na(values(rastAWS2)-values(rastAWS1))))*100
      points(htPP$breaks[2:length(htPP$breaks)],cumsum(htPP$counts),
             xlim=c(0,40),ylim=c(0,100),
             type='o',
             xlab='(V1-V0)/V0x100  (%)',
             ylab='Cumulative (%)',
             cex.lab=1.5,cex.axis=1.5,col='darkgreen',pch=19)
      
      legend('bottomright',c('SplinePar: 0.55','SplinePar: 0.40','Same code'),pch=19,cex=1.5,pt.cex=1.5,col=c('red','blue','darkgreen'),bty='n')
      
    }else{
      hist(values(rastAWS),breaks=seq(-500,600,1),xlim=c(-100,500),col=rgb(1,1,1),border=rgb(1,0,0,0.3))  
      hist(values(rastSCC),breaks=seq(-500,600,1),xlim=c(-100,500),col=rgb(1,1,1),border=rgb(0,0,1,0.3),add=T)  
      legend('topright',c('V0','V1'),pch=15,col=c('blue','red'),pt.cex=1.8,bty='n',cex=1.3)
      
      htPP  <- hist((values(rastAWS)-values(rastSCC)),breaks=seq(-100000,100000,5),plot=F)
      htPP$counts <- htPP$counts/(sum(!is.na(values(rastAWS)-values(rastSCC))))*100
      plot(htPP,xlim=c(-30,30),col=rgb(1,1,1),border=rgb(0,0,0),
           main='',xlab='V1 - V0',
           ylab='Frequency (%)',
           cex.lab=1.5,cex.axis=1.5)
      
      htPP  <- hist(abs(values(rastAWS)-values(rastSCC)),breaks=seq(-100000,100000,1),plot=F)
      htPP$counts <- htPP$counts/(sum(!is.na(values(rastAWS)-values(rastSCC))))*100
      plot(htPP$breaks[2:length(htPP$breaks)],cumsum(htPP$counts),
           xlim=c(0,30),ylim=c(0,100),
           type='o',
           xlab='V1 - V0 (Abs.)',
           ylab='Cumulative (%)',
           cex.lab=1.5,cex.axis=1.5,col='red',pch=19)
      
      htPP  <- hist(abs(values(rastAWS1)-values(rastSCC)),breaks=seq(-100000,100000,1),plot=F)
      htPP$counts <- htPP$counts/(sum(!is.na(values(rastAWS1)-values(rastSCC))))*100
      points(htPP$breaks[2:length(htPP$breaks)],cumsum(htPP$counts),
           xlim=c(0,30),ylim=c(0,100),
           type='o',
           xlab='V1 - V0 (Abs.)',
           ylab='Cumulative (%)',
           cex.lab=1.5,cex.axis=1.5,col='blue',pch=19)
      
      htPP  <- hist(abs(values(rastAWS2)-values(rastAWS1)),breaks=seq(-100000,100000,1),plot=F)
      htPP$counts <- htPP$counts/(sum(!is.na(values(rastAWS2)-values(rastAWS1))))*100
      points(htPP$breaks[2:length(htPP$breaks)],cumsum(htPP$counts),
             xlim=c(0,30),ylim=c(0,100),
             type='o',
             xlab='V1 - V0  (Abs.)',
             ylab='Cumulative (%)',
             cex.lab=1.5,cex.axis=1.5,col='darkgreen',pch=19)
      
      legend('bottomright',c('SplinePar: 0.55','SplinePar: 0.40','Same code'),pch=19,cex=1.5,pt.cex=1.5,col=c('red','blue','darkgreen'),bty='n')
    }
    
    
    dev.off()
    
    print(paste(year[yy],'; ',tiles[tt],'; ',var[pp]))
  }  
}



