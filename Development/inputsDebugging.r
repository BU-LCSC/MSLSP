tile <- '17SKV'
jsonFile  <- "/usr2/postdoc/dbolt/CodeFolder/git_repos/MSLSP/MSLSP_Parameters.json"
jsonFile <- "~/test.json"
errorLog <- '~/test.txt'
runLog <- '~/test2.txt'
timeStamp <- 'test'

registerDoMC(cores=8)

imgName <- imgList[80]



#Get default parameters
params <- fromJSON(file=jsonFile)
params$dirs$imgDir=paste0(params$SCC$dataDir,'S10/',tile,'/images/')
params$dirs$fmaskDir=paste0(params$SCC$dataDir,'HLS30/',tile,'/fmask/')
params$dirs$tempDir=paste0(params$SCC$workDir,tile,'/temp/')	
params$dirs$chunkDir=paste0(params$SCC$workDir,tile,'/imageChunks/')
params$dirs$phenDir=paste0(params$SCC$workDir,tile,'phenoMetrics/')



#Identify missing data in S10
S10 <- unique(list.files(path="~/dbolt/Runs_S10/18TYN/imageChunks/", pattern=glob2rx("*S10*.Rds"), full.names=F, recursive=T))
S30  <- unique(list.files(path="~/dbolt/run_it_29/18TYN/imageChunks/", pattern=glob2rx("*S30*.Rds"), full.names=F, recursive=T) )
S10_2 <- unique(unlist(strsplit(S10,'/'))[seq(2,length(S10)*2,by=2)])
S30_2 <- unique(unlist(strsplit(S30,'/'))[seq(2,length(S30)*2,by=2)])
S30_m <- gsub("S30","S10",S30_2)
notin <- S10_2[!(S10_2 %in% S30_m)]


tile <- '18TYP'
inDir <- "/projectnb/modislc/users/dbolt/V1_0/18TYP/images/"
outDir <- "/projectnb/modislc/users/dbolt/V1_0/18TYP/phenoMetrics/"

jsonFile  <- "/usr2/postdoc/dbolt/CodeFolder/git_repos/MuSLI_LSP/AWS/sh/V0_1/MSLSP_V1_Parameters.json"

errorLog <- '~/test.txt'



chunk=100
numPix <- numPixPerChunk[100]





#Compare time of V1 vs V0
############################

source('~/CodeFolder/git_repos/MSLSP/MSLSP_Functions.r')

start_time <- Sys.time()
#Loop through each pixel and estimate phenometrics
numPix <- 500
#Loop through each pixel and est  imate phenometrics
pheno_mat <- matrix(NA,numPix,pheno_pars$numLyrs*length(phenYrs))
for (i in 1:numPix) {pheno_mat[i,] <- DoPhenologyHLS(b2[i,],  b3[i,],  b4[i,],  b5[i,],  b6[i,], b7[i,],  vi[i,],
                                                     snowPix[i,],dates, imgYrs, phenYrs, splineStart, splineEnd, numDaysFit, pheno_pars)}
as.numeric(difftime(Sys.time(),start_time,units="secs"))



source('~/CodeFolder/git_repos/MuSLI_LSP/AWS/r/MuSLI_LSP_Functions.r')
pheno_pars2 <- DefaultPhenoParameters_HLS()

start_time <- Sys.time()
#Loop through each pixel and estimate phenometrics
numPix <- 500
#Loop through each pixel and est  imate phenometrics
pheno_mat <- matrix(NA,numPix,pheno_pars2$numLyrs*length(imgYrs))

for (i in 1:numPix) {pheno_mat[i,] <- DoPhenologyHLS_Climatology(b2[i,],  b3[i,],  b4[i,],  b5[i,],  b6[i,], b7[i,],  vi[i,],
                                                     snowPix[i,],dates, imgYrs, splineStart, numDaysFit, pheno_pars2)}
as.numeric(difftime(Sys.time(),start_time,units="secs"))


























source('~/CodeFolder/git_repos/MuSLI_LSP/AWS/r/V0.1/MuSLI_LSP_Functions_V0_1_Diagnostics.r')
i=1
out <- DoPhenologyHLS_Diagnostic(b2[i,],  b3[i,],  b4[i,],  b5[i,],  b6[i,], b7[i,],  vi[i,],
                                                     snowPix[i,],dates, imgYrs, phenYrs, splineStart, splineEnd, numDaysFit, pheno_pars)
as.numeric(difftime(Sys.time(),start_time,units="secs"))



source('~/CodeFolder/git_repos/MuSLI_LSP/AWS/r/V0.1/MuSLI_LSP_Functions_V0_2.r')

start_time <- Sys.time()
#Loop through each pixel and estimate phenometrics
numPix <- 500
#Loop through each pixel and est  imate phenometrics
pheno_mat <- matrix(NA,numPix,pheno_pars$numLyrs*length(phenYrs))
for (i in 1:numPix) {pheno_mat[i,] <- DoPhenologyHLS(b2[i,],  b3[i,],  b4[i,],  b5[i,],  b6[i,], b7[i,],  vi[i,],
                                                     snowPix[i,],dates, imgYrs, phenYrs, splineStart, splineEnd, numDaysFit, pheno_pars)}
as.numeric(difftime(Sys.time(),start_time,units="secs"))



source('~/CodeFolder/git_repos/MuSLI_LSP/AWS/r/V0.1/MuSLI_LSP_Functions_V0_1.r')

start_time <- Sys.time()
#Loop through each pixel and estimate phenometrics
numPix <- 500
#Loop through each pixel and est  imate phenometrics
pheno_mat <- matrix(NA,numPix,pheno_pars$numLyrs*length(phenYrs))
for (i in 1:numPix) {pheno_mat[i,] <- DoPhenologyHLS(b2[i,],  b3[i,],  b4[i,],  b5[i,],  b6[i,], b7[i,],  vi[i,],
                                                     snowPix[i,],dates, imgYrs, phenYrs, splineStart, splineEnd, numDaysFit, pheno_pars)}
as.numeric(difftime(Sys.time(),start_time,units="secs"))















i <- 50
b2 <- b2[i,]; b3 <- b3[i,]; b4 <- b4[i,]
b5 <- b5[i,]; b6 <- b6[i,]; b7 <- b7[i,]
snowPix <- snowPix[i,]
vi <- vi[i,]





source('~/CodeFolder/git_repos/MuSLI_LSP/AWS/r/V0.1/MuSLI_LSP_Functions_V0_2.r')

start_time <- Sys.time()
for (i in 1:500) {
  DoPhenologyHLS(b2,  b3,  b4,  b5,  b6, b7,  vi,
                 snowPix,dates, imgYrs, phenYrs, splineStart, splineEnd, numDaysFit, pheno_pars)}
as.numeric(difftime(Sys.time(),start_time,units="secs"))


source('~/CodeFolder/git_repos/MuSLI_LSP/AWS/r/V0.1/MuSLI_LSP_Functions_V0_1.r')

start_time <- Sys.time()
for (i in 1:500) {
  DoPhenologyHLS(b2,  b3,  b4,  b5,  b6, b7,  vi,
                 snowPix,dates, imgYrs, phenYrs, splineStart, splineEnd, numDaysFit, pheno_pars)}
as.numeric(difftime(Sys.time(),start_time,units="secs"))




































start_time <- Sys.time()
for (i in 1:500) {
DoPhenologyHLS(b2,  b3,  b4,  b5,  b6, b7,  vi,
               snowPix,dates, imgYrs, phenYrs, splineStart, splineEnd, numDaysFit, pheno_pars)}
as.numeric(difftime(Sys.time(),start_time,units="secs"))













i=4
DoPhenologyHLS(b2[i,],  b3[i,],  b4[i,],  b5[i,],  b6[i,], b7[i,],  vi[i,],
               snowPix[i,],dates, imgYrs, phenYrs, splineStart, splineEnd, numDaysFit, pheno_pars)










yrs=2016:2018
pheno_pars <- DefaultPhenoParameters_HLS() 
yrs <- seq(startYear,endYear)
splineStart <- as.Date(as.Date(paste0(yrs,'-01-01')) - pheno_pars$splineBuffer) 
numDaysFit  <-  365 + (pheno_pars$splineBuffer * 2)
pheno_pars$halfLyrs <- (pheno_pars$numLyrs-1)/2    #Number of layers that require filling if there is no second cycle. First cycle has 7 extra layers

start_time <- Sys.time()
#Loop through each pixel and estimate phenometrics
numPix <- 1000
#Loop through each pixel and est  imate phenometrics
pheno_mat <- matrix(NA,numPix,pheno_pars$numLyrs*length(yrs))
for (i in 1:numPix) {pheno_mat[i,] <- DoPhenologyHLS_Climatology(b2[i,],  b3[i,],  b4[i,],  b5[i,],  b6[i,], b7[i,],  vi[i,],
                                                     snowPix[i,],dates, yrs, splineStart, numDaysFit, pheno_pars)}
as.numeric(difftime(Sys.time(),start_time,units="secs"))










#Old functions
source('~/CodeFolder/git_repos/MuSLI_LSP/AWS/r/MuSLI_LSP_Functions.r')








