#Load required libraries
library(raster)
library(rgdal)
library(gdalUtils)
library(rgeos)
library(foreach)
library(iterators)
library(doMC)

library(matrixStats)
library(WGCNA)
library(zoo)

library(RcppRoll)

library(rjson)

library(ncdf4)



#Inputs
############


inFold <- '/projectnb/modislc/users/dbolt/Crop_Analysis/anomTiles/'


crops=c('Corn','Soybeans')
for (crop in crops) {

phen2017 <- read.csv(paste0(inFold,'Phen_Dates_',crop,'_2017.csv'),header=T)
phen2018 <- read.csv(paste0(inFold,'Phen_Dates_',crop,'_2018.csv'),header=T)
phen2019 <- read.csv(paste0(inFold,'Phen_Dates_',crop,'_2019.csv'),header=T)

smooth2017 <- read.csv(paste0(inFold,'Smoothed_timeseries_',crop,'_2017.csv'),header=F)
smooth2018 <- read.csv(paste0(inFold,'Smoothed_timeseries_',crop,'_2018.csv'),header=F)
smooth2019 <- read.csv(paste0(inFold,'Smoothed_timeseries_',crop,'_2019.csv'),header=F)


p_2017 <- apply(smooth2017[,-1],1,quantile,probs=c(0.25,0.5,0.75),na.rm=T)
p_2018 <- apply(smooth2018[,-1],1,quantile,probs=c(0.25,0.5,0.75),na.rm=T)
p_2019 <- apply(smooth2019[,-1],1,quantile,probs=c(0.25,0.5,0.75),na.rm=T)


pcgi_2017 <- as.Date(quantile(phen2017$PCGI,probs=c(0.25,0.5,0.75),na.rm=T))+365
pcgi_2018 <- as.Date(quantile(phen2018$PCGI,probs=c(0.25,0.5,0.75),na.rm=T))
pcgi_2019 <- as.Date(quantile(phen2019$PCGI,probs=c(0.25,0.5,0.75),na.rm=T))-365

pcgd_2017 <- as.Date(quantile(phen2017$PCGD,probs=c(0.25,0.5,0.75),na.rm=T))+365 
pcgd_2018 <- as.Date(quantile(phen2018$PCGD,probs=c(0.25,0.5,0.75),na.rm=T)) 
pcgd_2019 <- as.Date(quantile(phen2019$PCGD,probs=c(0.25,0.5,0.75),na.rm=T)) -365


dates <- as.Date(smooth2018[,1])

polyDates <- c(dates,rev(dates))
poly17 <- c(p_2017[1,],rev(p_2017[3,]))
poly18 <- c(p_2018[1,],rev(p_2018[3,]))
poly19 <- c(p_2019[1,],rev(p_2019[3,]))

polyCol <- c('#78c679','#43a2ca','#feb24c')
transCol <- paste0(polyCol,'B3')



figName <- paste0(inFold,crop,'_2018_2019.png')
png(figName,res=300,width = 6,height = 3,units = "in",pointsize = 8)


xlims <- as.Date(c('2018-02-01','2018-12-01'))  #Always use the first of a month! (for labelling purposes)
plot(p_2017[2,]~dates,xlim=xlims,type='n',xlab='Date',ylab='EVI2',ylim=c(0.1,0.8))
#polygon(polyDates,poly17,col=transCol[1],border=polyCol[1])
polygon(polyDates,poly18,col=transCol[2],border=polyCol[2])
polygon(polyDates,poly19,col=transCol[3],border=polyCol[3])

axis.Date(1,at=seq(xlims[1], xlims[2],by="months"),
          labels=format(seq(xlims[1], xlims[2],by="months"),'%b'),
          tck=0.025)


abline(v=pcgi_2018[2],col=polyCol[2],lty='dashed')
abline(v=pcgd_2018[2],col=polyCol[2],lty='dashed')
abline(v=pcgi_2019[2],col=polyCol[3],lty='dashed')
abline(v=pcgd_2019[2],col=polyCol[3],lty='dashed')

legend('topleft',legend=c('2018','2019','50% threshold dates'),
       border=c(polyCol[2:3],NA),fill=c(transCol[2:3],NA),lty=c(NA,NA,'dashed'),bty='n')

dev.off()

}

