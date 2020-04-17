

#install.packages('rasterVis')


#Load required libraries
library(raster)
library(rgdal)
library(gdalUtils)
library(rgeos)
library(maptools)
library(rasterVis)

library(iterators)
library(foreach)
library(doMC)

library(ncdf4)



numCores=28

yrs <- 2016:2018

#Register the parallel backend
registerDoMC(cores=numCores)

inBase <- '/projectnb/modislc/users/dbolt/AWS_results/tiles/'

outBase <- '/projectnb/modislc/users/mkmoon/MuSLI/V0_1/' 


lyrs <- read.csv('~/MSLSP/Development/make_netCDF/MSLSP_Layers_V0.csv',header=T,stringsAsFactors = F)

lyrs <- lyrs[!is.na(lyrs$product_lyr),]
lyrs <- lyrs[order(lyrs$product_lyr),]


for (yr in yrs) {

  inDir <- paste0(inBase,yr,'/')
  tiles <- list.dirs(inDir,recursive=F,full.names=F)

  imgLog <- foreach(t=1:length(tiles),.combine=c) %dopar% {
      
        tile <- tiles[t]
        inRast_Base <- paste0(inDir,tile,'/LSP_',tile,'_',yr,'_')  
        
        outFold <- paste0(outBase, tile, '/')
        outFile <- paste0(outFold,'MSLSP_',tile,'_',yr,'.nc')
        if (!dir.exists(outFold)) {dir.create(outFold)}
        
        
        #Get a base imageto pull raster info from     
        baseImage <- raster(paste0(inRast_Base, 'EVImax.tif'))
        
        #Get extent, and then define pixel centers in the x and y direction
        ext = extent(baseImage)
        res = res(baseImage)[1]
        x = seq(ext[1]+res/2,ext[2]-res/2, res)
        y = seq(ext[3]+res/2,ext[4]-res/2, res)
        
        #Define dimensions for netCDF file
        dimx = ncdim_def(name = 'x', longname = 'x coordinate', units='m', vals = as.double(x))
        dimy = ncdim_def(name = 'y', longname = 'y coordinate', units='m', vals = rev(as.double(y)))

        
        #Define a projection variable for the file
        prj_def <- ncvar_def("transverse_mercator","",list(),prec="char")
        
        
        #Loop through all the layers, and create a variable for each 
        #Make the first variable the projection info from above
        results<-vector("list", dim(lyrs)[1]+1) 
        results[[1]] <- ncvar_def("transverse_mercator","",list(),prec="char")
        
        for (i in 1:dim(lyrs)[1]) {
          
          lyr <- lyrs[i,]   #Pull the info for this layer from the lyrs table
          
          if (lyr$data_type == 'Int16') {precision <- "short"}  #All are int16, so this isn't necessary
          
          #Create the variable, add to the list. Define the short_name, units, fill_value, long_name, and precision from the lyrs table
          results[[i+1]] <- ncvar_def(lyr$short_name, lyr$units, list(dimx,dimy), lyr$fill_value, lyr$long_name, prec=precision, compression=2)  
        }
        
        
        # Now create the netCDF file with the defined variables
        if (file.exists(outFile)) {file.remove(outFile)};
        ncout <- nc_create(outFile,results,force_v4=T)
        
        
        #Now loop through the layers again, this time actually 
        #writing the image data to the file
        for (i in 1:dim(lyrs)[1]) {
          
          lyr <- lyrs[i,]
          
          #Open the data. Current solution to "WRITE_BOTTOMUP" is 
          #to manually flip the image 
          mat <- matrix(as.numeric(readGDAL(paste0(inRast_Base,lyr$short_name,'.tif'),silent=T)$band1),length(x),length(y))
          ncvar_put(ncout,results[[i+1]], mat)      #Now put the image into the file
          
          #Fill in the attributes for the layer from the lyrs table
          ncatt_put(ncout,lyr$short_name,"scale",lyr$scale)
          ncatt_put(ncout,lyr$short_name,"offset",lyr$offset)
          ncatt_put(ncout,lyr$short_name,"data_type",lyr$data_type)
          ncatt_put(ncout,lyr$short_name,"valid_min",lyr$valid_min)
          ncatt_put(ncout,lyr$short_name,"valid_max",lyr$valid_max)
          ncatt_put(ncout,lyr$short_name,"grid_mapping","transverse_mercator")  
          
        }
        
        
        #Now we need to write the projection info the the transverse_mercator variable
        
        wkt <- showWKT(projection(baseImage), morphToESRI = FALSE)  #Get projection in wkt format
        wkt2 <-capture.output(cat(wkt))
        #Need to pull the central meridian from the wkt (Hopefull this will always work?)
        spt <- unlist(strsplit(gsub(']','',wkt),','))
        central_meridian <- as.numeric(spt[which(spt == "PARAMETER[\"central_meridian\"")+1])
        
        #Fill in the info. Manaully wrote values that don't change across UTM zones (Is that correct?)
        ncatt_put(ncout,"transverse_mercator","grid_mapping_name","transverse_mercator")
        ncatt_put(ncout,"transverse_mercator","longitude_of_central_meridian",central_meridian)
        ncatt_put(ncout,"transverse_mercator","false_easting",5e+05)
        ncatt_put(ncout,"transverse_mercator","false_northing",0)
        ncatt_put(ncout,"transverse_mercator","latitude_of_projection_origin",0)
        ncatt_put(ncout,"transverse_mercator","scale_factor_at_central_meridian",0.9996)
        ncatt_put(ncout,"transverse_mercator","long_name","CRS definition")
        ncatt_put(ncout,"transverse_mercator","longitude_of_prime_meridian",0)
        ncatt_put(ncout,"transverse_mercator","semi_major_axis",6378137)
        ncatt_put(ncout,"transverse_mercator","inverse_flattening",298.257223563)
        ncatt_put(ncout,"transverse_mercator","spatial_ref",gsub("\\", "", wkt, fixed=TRUE))
        ncatt_put(ncout,"transverse_mercator","GeoTransform",paste(ext[1],res,0,ext[4],0,-res))  #I think this should work?
        
        
        
        #Define global attributes
        ncatt_put(ncout,0,"title","Multisource Land Surface Phenology (MS-LSP)")
        ncatt_put(ncout,0,"product_version","v0.1")
        ncatt_put(ncout,0,"summary","A 30m Land Surface Phenology Product for North America derived from Landsat and Sentinel-2 imagery")
        
        #ncatt_put(ncout,0,"algorithm theoretical basis document","Under development")
        #ncatt_put(ncout,0,"user guide","Under development")
        ncatt_put(ncout,0,"software_repository","git@github.com:BU-LCSC/MSLSP.git")
        
        ncatt_put(ncout,0,"reference","Bolton, D.K., Gray, J.M., Melaas, E.K., Moon, M., Eklundh, L, Friedl, M.A., 2020. Continental-Scale Land Surface Phenology from Harmonized Landsat 8 and Sentinel-2 Imagery. Remote Sensing of Environment")
        
        ncatt_put(ncout,0,"program","The Multi-Source Land Imaging (MuSLI) project funded by NASA's Land Cover Land Use Change (LCLUC) program")
        
        ncatt_put(ncout,0,"creator_name","Land Cover & Surface Climate Group, Department of Earth & Environment, Boston University")
        ncatt_put(ncout,0,"creator_type","group")
        ncatt_put(ncout,0,"creator_email","friedl@bu.edu")
        ncatt_put(ncout,0,"creator_institution","Boston University")
        
        ncatt_put(ncout,0,"contributor_name", "Douglas K. Bolton, Mark A. Friedl, Josh M. Gray, Lars Eklundh, Eli M. Melaas, Minkyu M. Moon")
        ncatt_put(ncout,0,"contributor_role", "Developer, Principal Investigator, Co-Investigator, Collaborator, Contributor, Contributor")
        
        ncatt_put(ncout,0,"acknowledgement","Developed with funding from NASA LCLUC Grant #80NSSC18K0334. Data archiving and distribution supported by the NASA NSIDC Distributed Active Archive Center (DAAC).")
        
        #Put additional attributes on coordinates
        ncatt_put(ncout,"x","standard_name","projection_x_coordinate")
        ncatt_put(ncout,"y","standard_name","projection_y_coordinate")
        
        nc_close(ncout)
        
  
  }
}