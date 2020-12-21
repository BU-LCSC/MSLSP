library(sp)
library(raster)
library(ncdf4)

args <- commandArgs()
print(args)

tt <- args[3]
# tt <- '18TXS'

path <- '/projectnb/modislc/users/mkmoon/MuSLI/V1_0/From_AWS/product'
files <- list.files(path=paste(path,'/',tt,sep=''),pattern=glob2rx('*.nc'),full.names=T)

years <- 2016:2019
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


# ###########################
# aa <- raster('/projectnb/modislc/users/dbolt/AWS_results/tiles/2017/15WWS/LSP_15WWS_2017_50PCGD_2.tif')
# bb <- values(aa)
# sum(!is.na(bb))
# plot(aa,colNA='gray')
# 
# 
# aa <- raster(matrix(NA,1000,1000))
# aa[1] <- 1
# plot(aa,colNA='gray')


# files <- list.files('/projectnb/modislc/users/mkmoon/lst/data/merra',
#                     pattern=glob2rx('MERRA2_300*'),full.names=T)
# merra <- nc_open(files[1])
# varmerra <- names(merra[['var']])
# ramerra <- brick(files[1],varname=varmerra[1])
# 
# daymet <- nc_open('/projectnb/modislc/data/climate/daymet/daymet_v3_dayl_1981_na.nc4')
# vardaymet <- names(daymet[['var']])
# radaymet <- brick('/projectnb/modislc/data/climate/daymet/daymet_v3_dayl_1981_na.nc4',varname=vardaymet[6])
# plot(radaymet)
