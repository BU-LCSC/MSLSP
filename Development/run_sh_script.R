setwd('/projectnb/modislc/users/mkmoon/Planet/data/HLS_fusion/')
for(tt in 1:10){
  system(paste('qsub -V -pe omp 2 -l h_rt=12:00:00 /usr3/graduate/mkmoon/GitHub/MSLSP/Development/run_script.sh ',tt,sep=''))  
}

setwd('/projectnb/modislc/users/mkmoon/Planet/data/HLS_fusion/')
for(tt in 1:10){
  system(paste('qsub -V -pe omp 2 -l h_rt=12:00:00 /usr3/graduate/mkmoon/GitHub/MSLSP/Development/run_script_1.sh ',tt,sep=''))  
}


# Check AWS results
setwd('/projectnb/modislc/users/mkmoon/MuSLI/V1_0/figures/')
for(tt in 1:3){
  system(paste('qsub -V -pe omp 2 -l h_rt=06:00:00 /usr3/graduate/mkmoon/GitHub/MSLSP/Development/run_script.sh ',tt,sep=''))  
}
