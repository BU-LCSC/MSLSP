
#module load R/3.3.2
#module load python/2.7.13
#module load awscli


#This will only work if you have your AWS key set up. Ask Dennis for help!


args <- commandArgs(trailingOnly=T)
print(args)
tile <- args[1]
outDir <- args[2]
imgStartYr <- as.numeric(args[3])
imgEndYr <- as.numeric(args[4])

yrs <- imgStartYr:imgEndYr

for (yr in yrs) {

  path <- paste0("s3://hlsanc/PRO/v1.4/S2/S10/",yr,'/',
                 substring(tile,1,2),'/',substring(tile,3,3),'/',
                 substring(tile,4,4),'/',substring(tile,5,5),'/')
  #line <- paste("module load python2/2.7.16;module load awscli;aws s3 ls",path,"--request-payer requester | awk '{print $4}'")
  line <- paste("aws s3 ls",path,"--request-payer requester | awk '{print $4}'")
  list <- system(line,intern=T)
  
  for (file in list) {
    outPath <- paste0(outDir,file)
    if (!file.exists(outPath)) {
      inPath <- paste0(path,file)
      line <- paste("aws s3 cp",inPath,outPath,'--request-payer requester')
      system(line)}
  }
}
      
    