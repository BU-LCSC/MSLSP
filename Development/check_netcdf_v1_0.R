library(sp)
library(raster)
library(ncdf4)

args <- commandArgs()
print(args)

tt <- as.numeric(args[3])

# tt <- 1
years <- 2016:2019

path <- '/projectnb/modislc/users/mkmoon/MuSLI/V1_0/From_AWS'
tiles <- substr(list.dirs(path=path,full.names=T)[2:25],53,57)
files <- list.files(path=paste(path,'/',tiles[tt],sep=''),pattern=glob2rx('*.nc'),full.names=T)

setwd('/projectnb/modislc/users/mkmoon/NEphenology/data/bay/figure/')
erepa <- shapefile('/projectnb/modislc/users/mkmoon/NEphenology/ecoregions/na_cec_eco_l1/NA_CEC_Eco_Level1.shp')

for(ff in 1:4){
  nc <- nc_open(files[ff*2])
  var <- names(nc[['var']])
  
  options(warn=-1)
  
  setwd(paste('/projectnb/modislc/users/mkmoon/MuSLI/V1_0/rasters/',years[ff],sep=''))
  # pdf(file=paste(tiles[tt],'_',years[ff],'.pdf',sep=''),width=15,height=8)
  png(filename=paste(tiles[tt],'_',years[ff],'.png',sep=''),width=15,height=8,units='in',res=150)
  par(mfrow=c(4,6),oma=c(0,0,0,2),mar=c(0,1,2,4))
  for(vv in c(2:5,7:12,55:59,61:66,109:110)){
    rast <- brick(files[ff*2],varname=var[vv])
    temp <- raster(files[ff*2],varname=var[vv])
    nv <- sum(!is.na(values(temp)))
    plot(rast,axes=F,box=F,
         main=paste(years[ff],'_',var[vv],'_',nv,sep=''),
         colNA='grey45',cex.main=1.2)
    print(vv)
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
