#######################################################
# Copy from AWS to SCC

# tiles <- c('13TEF','16TCK','17SQD','18TYN','19TEL')
tiles <- c('18TXS','18TXT','18TYL','18TYM','18TYN','18TYP','18TYQ','18TYR','18TYS','18TYT') # 12/20/2020
tiles <- c('13SBA','13SBB','13SBC','13SBD','13SBR','13SBS','13SBT','13SBU','13SBV','13SCA', # 12/21/2020
           '15TTE','15TTF','15TTG','15TUE','15TUF') 
tiles <- c('16TBK','16TBL','16TBM','16TCK','16TCL') # 1/12/2021
            


# Recursive - phenometrics file only
for(i in 1:length(tiles)){
  for(yy in 2016:2019){
    system(paste('aws s3 cp s3://mslsp/v1.0/product/',tiles[i],'/MSLSP_',tiles[i],'_',yy,'.nc',' /projectnb/modislc/users/mkmoon/MuSLI/V1_0/From_AWS/product/',tiles[i],'/ ',sep=''))  
  }
}

# Recursive - all files
for(i in 1:length(tiles)){
  system(paste('aws s3 cp s3://mslsp/v1.0/product/',tiles[i],' /projectnb/modislc/users/mkmoon/MuSLI/V1_0/From_AWS/product/',tiles[i],'/ --recursive',sep=''))
}

# # V0 log info
# system('aws s3 cp s3://mslsp/logs/instanceInfo/ /projectnb/modislc/users/mkmoon/MuSLI/V1_0/product_qc/v0_loginfo/ --recursive' )



#######################################################
# Submit jobs for quality check

# Check AWS results - rasters
setwd('/projectnb/modislc/users/mkmoon/MuSLI/V1_0/product_qc/')
for(tt in 1:length(tiles)){
  system(paste('qsub -V -pe omp 2 -j y -l h_rt=06:00:00 /usr3/graduate/mkmoon/GitHub/MSLSP/Development/run_script.sh ',tiles[tt],sep=''))  
}


