setwd('/projectnb/modislc/users/mkmoon/Planet/data/HLS_fusion/')
for(tt in 1:10){
  system(paste('qsub -V -pe omp 2 -l h_rt=12:00:00 /usr3/graduate/mkmoon/GitHub/MSLSP/Development/run_script.sh ',tt,sep=''))  
}

setwd('/projectnb/modislc/users/mkmoon/Planet/data/HLS_fusion/')
for(tt in 1:10){
  system(paste('qsub -V -pe omp 2 -l h_rt=12:00:00 /usr3/graduate/mkmoon/GitHub/MSLSP/Development/run_script_1.sh ',tt,sep=''))  
}


