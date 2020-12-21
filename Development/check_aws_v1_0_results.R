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

tt <- args[3]
# tt <- '18TYL'
years <- 2016:2019

###############################################
# All layers as raster
files <- list.files(path=paste(path,'/',tt,sep=''),pattern=glob2rx('*.nc'),full.names=T)
erepa <- shapefile('/projectnb/modislc/users/mkmoon/NEphenology/ecoregions/na_cec_eco_l1/NA_CEC_Eco_Level1.shp')

for(i in 1:length(years)){
  nc <- nc_open(files[i])
  var <- names(nc[['var']])
  
  options(warn=-1)
  
  setwd('/projectnb/modislc/users/mkmoon/MuSLI/V1_0/product_qc/raster/')
  png(filename=paste(tt,'_',years[i],'.png',sep=''),width=16,height=8,units='in',res=150)
  par(mfrow=c(4,7),oma=c(0,0,0,2),mar=c(0,1,2,4))
  
  for(j in 2:26){
    rast <- raster(files[i],varname=var[j])
    nv <- sum(!is.na(values(rast)))
    plot(rast,axes=F,box=F,
         main=paste(years[i],'_',var[j],'_',nv,sep=''),
         colNA='grey45',cex.main=1.2)
    print(paste(i,';',j))
  }
  
  pr3 <- projectExtent(rast,crs(erepa))
  temp <- crop(erepa,pr3)
  
  plot(erepa,col='grey90',border='grey90')
  plot(temp,add=T,col='red',border='red')
  
  dev.off()
}


###############################################
# Comparison with V0
spec <- viridis(11)
mycolRamp = colorRampPalette(c('White',spec))

for(yy in 1){
  
  pathSCC <- paste0('/projectnb/modislc/users/mkmoon/MuSLI/V0_11/',tt)             # V0 code with SplinePar 0.55
  pathAWS <- paste0('/projectnb/modislc/users/mkmoon/MuSLI/V1_0/From_AWS/product/',tt)  # 'master' as of 12/19/20
  
  sstr <- paste0('*',tt,'*',years[yy],'.nc')
  fileSCC <- list.files(pathSCC,pattern=glob2rx(sstr),recursive=F,full.names=T)
  fileAWS <- list.files(pathAWS,pattern=glob2rx(sstr),recursive=F,full.names=T)
  
  nc <- nc_open(fileSCC)
  var <- names(nc[['var']])
  
  setwd('/projectnb/modislc/users/mkmoon/MuSLI/V1_0/product_qc/diff/')
  for(pp in 3:11){
    rastSCC <- raster(fileSCC,varname=var[pp])
    rastAWS <- raster(fileAWS,varname=var[pp])  
   
    png(filename=paste('diff_MSLSP_',tt,'_',years[yy],'_',sprintf('%02d',pp),'.png',sep=''),
        width=12,height=7,unit='in',res=150)
    
    par(mfrow=c(2,3),oma=c(1,1,1,1),mar=c(4,4,1,1),mgp=c(2.5,1,0))
    plot(rastSCC,zlim=c(quantile(values(rastAWS),0,na.rm=T),quantile(values(rastAWS),1,na.rm=T)),
         main=paste(var[pp],'_V0_',
                            round(sum(!is.na(values(rastSCC)))/sum(!is.na(values(rastAWS)))*100,0),'%',
                            sep=''),
         cex.main=1.5)
    plot(rastAWS,zlim=c(quantile(values(rastAWS),0,na.rm=T),quantile(values(rastAWS),1,na.rm=T)),
         main=paste(var[pp],'_V1_',sum(!is.na(values(rastAWS))),sep=''),cex.main=1.5)
    
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
      
      htPP  <- hist(abs(values(rastAWS)-values(rastSCC)),breaks=seq(-100000,100000,10),plot=F)
      htPP$counts <- htPP$counts/(sum(!is.na(values(rastAWS)-values(rastSCC))))*100
      plot(htPP$breaks[2:length(htPP$breaks)],cumsum(htPP$counts),
           xlim=c(0,2000),ylim=c(0,100),
           type='o',
           xlab='V1 - V0  (Abs.)',
           ylab='Cumulative (%)',
           cex.lab=1.5,cex.axis=1.5,col='red',pch=19)
      
      htPP  <- hist((abs((values(rastAWS)-values(rastSCC))/values(rastSCC)*100)),breaks=seq(-100000,100000,1),plot=F)
      htPP$counts <- htPP$counts/(sum(!is.na(values(rastAWS)-values(rastSCC))))*100
      plot(htPP$breaks[2:length(htPP$breaks)],cumsum(htPP$counts),
           xlim=c(0,40),ylim=c(0,100),
           type='o',
           xlab='(V1-V0)/V0x100  (%)',
           ylab='Cumulative (%)',
           cex.lab=1.5,cex.axis=1.5,col='red',pch=19)
      
      abline(v=10,lty=5)
      
    }else if(pp>8 & pp<11){
      hist(values(rastAWS),breaks=seq(-20000,50000,10),xlim=c(0,15000),col=rgb(1,1,1),border=rgb(1,0,0,0.3))  
      hist(values(rastSCC),breaks=seq(-20000,50000,10),xlim=c(0,15000),col=rgb(1,1,1),border=rgb(0,0,1,0.3),add=T)   
      legend('topright',c('V0','V1'),pch=15,col=c('blue','red'),pt.cex=1.8,bty='n',cex=1.3)
    
      htPP  <- hist(abs(values(rastAWS)-values(rastSCC)),breaks=seq(-1000000,1000000,10),plot=F)
      htPP$counts <- htPP$counts/(sum(!is.na(values(rastAWS)-values(rastSCC))))*100
      plot(htPP$breaks[2:length(htPP$breaks)],cumsum(htPP$counts),
           xlim=c(0,1000),ylim=c(0,100),
           type='o',
           xlab='V1 - V0  (Abs.)',
           ylab='Cumulative (%)',
           cex.lab=1.5,cex.axis=1.5,col='red',pch=19)
      
      
      htPP  <- hist((abs((values(rastAWS)-values(rastSCC))/values(rastSCC)*100)),breaks=seq(-100000,1000000,1),plot=F)
      htPP$counts <- htPP$counts/(sum(!is.na(values(rastAWS)-values(rastSCC))))*100
      plot(htPP$breaks[2:length(htPP$breaks)],cumsum(htPP$counts),
           xlim=c(0,40),ylim=c(0,100),
           type='o',
           xlab='(V1-V0)/V0x100  (%)',
           ylab='Cumulative (%)',
           cex.lab=1.5,cex.axis=1.5,col='red',pch=19)
      
      abline(v=10,lty=5)
      
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
    }
    
    dev.off()
    print(paste(years[yy],'; ',tt,'; ',var[pp]))
  }  
}



