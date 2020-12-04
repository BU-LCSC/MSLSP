library(raster)
library(sp)
library(gdalUtils)

tt <- 3
# for(tt in 1:8){
tile_list <- as.matrix(read.table('/usr3/graduate/mkmoon/GitHub/MSLSP/SCC/tileLists/europe.txt'))

path <- paste('/projectnb/modislc/projects/landsat_sentinel/v1_4/HLS30/',tile_list[tt],'/images',sep='')
sstr <- paste('HLS*',tile_list[tt],'*v1.4.hdf',sep='')
file <- list.files(path=path,pattern=glob2rx(sstr),full.names=T)
sds <- get_subdatasets(file[1])
img <- raster(sds[1])

shp_hls <- shapefile('/projectnb/modislc/projects/landsat_sentinel/shapefiles/sentinel2_tiles_world/sentinel2_tiles_world/sentinel2_tiles_world.shp')
shp_hls <- shp_hls[as.character(shp_hls$Name) %in% tile_list[tt],]
shp_hls <- spTransform(shp_hls,crs(img))

# setwd('/projectnb/modislc/projects/landsat_sentinel/eu_dem/hls_tiles/shapefile/')
# shapefile(shp_hls,filename=paste('shp_',tile_list[tt],sep=''),overwrite=T)
# }

dem <- raster('/projectnb/modislc/projects/landsat_sentinel/eu_dem/eu_dem_v11_E30_40N30.TIF')
dem1 <- raster('/projectnb/modislc/projects/landsat_sentinel/eu_dem/eu_dem_v11_E30N30.TIF')
dem2 <- raster('/projectnb/modislc/projects/landsat_sentinel/eu_dem/eu_dem_v11_E40N30.TIF')
# wtr <- raster('/projectnb/modislc/projects/landsat_sentinel/eu_dem/Hansen_GFC2013_datamask_6070N_010020E.tif')
wtr <- raster('/projectnb//modislc/data/water/original_data/Hansen_GFC2013_datamask_60N_000E.tif/Stornext/scienceweb1/development/gtc/downloads/WaterMask2010_UMD/Hansen_GFC2013_datamask_60N_000E.tif')

pr3 <- projectExtent(img,crs(img))
dem1 <- projectRaster(dem1,pr3)
dem2 <- projectRaster(dem2,pr3)
dem <- merge(dem1,dem2)

pr2 <- projectExtent(img,crs(wtr))
wtr <- crop(wtr,pr2)
wtr <- projectRaster(wtr,pr3,method='ngb')

asp <- terrain(dem,opt='aspect')*10000
slp <- terrain(dem,opt='slope')*10000

# sdem <- raster('/projectnb/modislc/data/dem/usgs_ned/hls_tiles/dem/dem_18STE.tif')
# sasp <- raster('/projectnb/modislc/data/dem/usgs_ned/hls_tiles/aspect/aspect_18STE.tif')
# sslp <- raster('/projectnb/modislc/data/dem/usgs_ned/hls_tiles/slope/slope_18STE.tif')
# swtr <- raster('/projectnb/modislc/data/water/hls_tiles/water_33VVC.tif')
# 
# asp <- terrain(sdem,opt='aspect')*10000
# slp <- terrain(sdem,opt='slope')*10000


# Save files
setwd('/projectnb/modislc/projects/landsat_sentinel/eu_dem/hls_tiles/dem/')
writeRaster(dem,filename=paste('dem_',tile_list[tt],'.tif',sep=''),format="GTiff",overwrite=T)

setwd('/projectnb/modislc/projects/landsat_sentinel/eu_dem/hls_tiles/aspect/')
writeRaster(asp,filename=paste('aspect_',tile_list[tt],'.tif',sep=''),format="GTiff",overwrite=T)

setwd('/projectnb/modislc/projects/landsat_sentinel/eu_dem/hls_tiles/slope/')
writeRaster(slp,filename=paste('slope_',tile_list[tt],'.tif',sep=''),format="GTiff",overwrite=T)

setwd('/projectnb/modislc/projects/landsat_sentinel/eu_dem/hls_tiles/water/')
writeRaster(wtr,filename=paste('water_',tile_list[tt],'.tif',sep=''),format="GTiff",overwrite=T)


# #############
# par(mfrow=c(3,3))
# for(i in 1:8){
#   temp <- raster(paste('/projectnb/modislc/projects/landsat_sentinel/eu_dem/hls_tiles/water/water_',tile_list[i],'.tif',sep=''))
# plot(temp)
#   }
  