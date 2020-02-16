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


inFold <- '/projectnb/modislc/users/dbolt/Crop_Analysis/kansasTiles/'

crop='Corn'


progress <- read.csv(paste0(inFold,'kansas_crop_progress_2018_2019.csv'),header=T,stringsAsFactors = F)

progress$Week.Ending = as.Date(progress$Week.Ending)

#Pull 50% date for all metrics for each year

uniq = unique(progress$Data.Item)

planted_2018 = min(progress$Week.Ending[progress$Year == 2018 & progress$Data.Item == "CORN - PROGRESS, MEASURED IN PCT PLANTED" & progress$Value > 50])
emerged_2018 = min(progress$Week.Ending[progress$Year == 2018 & progress$Data.Item == "CORN - PROGRESS, MEASURED IN PCT EMERGED" & progress$Value > 50])
silking_2018 = min(progress$Week.Ending[progress$Year == 2018 & progress$Data.Item == "CORN - PROGRESS, MEASURED IN PCT SILKING" & progress$Value > 50])
dough_2018 = min(progress$Week.Ending[progress$Year == 2018 & progress$Data.Item == "CORN - PROGRESS, MEASURED IN PCT DOUGH" & progress$Value > 50])
dent_2018 = min(progress$Week.Ending[progress$Year == 2018 & progress$Data.Item == "CORN - PROGRESS, MEASURED IN PCT DENTED" & progress$Value > 50])
mature_2018 = min(progress$Week.Ending[progress$Year == 2018 & progress$Data.Item == "CORN - PROGRESS, MEASURED IN PCT MATURE" & progress$Value > 50])
harvested_2018 = min(progress$Week.Ending[progress$Year == 2018 & progress$Data.Item == "CORN, GRAIN - PROGRESS, MEASURED IN PCT HARVESTED" & progress$Value > 50])


planted_2019 = min(progress$Week.Ending[progress$Year == 2019 & progress$Data.Item == "CORN - PROGRESS, MEASURED IN PCT PLANTED" & progress$Value > 50]) - 365
emerged_2019 = min(progress$Week.Ending[progress$Year == 2019 & progress$Data.Item == "CORN - PROGRESS, MEASURED IN PCT EMERGED" & progress$Value > 50])- 365
silking_2019 = min(progress$Week.Ending[progress$Year == 2019 & progress$Data.Item == "CORN - PROGRESS, MEASURED IN PCT SILKING" & progress$Value > 50])- 365
dough_2019 = min(progress$Week.Ending[progress$Year == 2019 & progress$Data.Item == "CORN - PROGRESS, MEASURED IN PCT DOUGH" & progress$Value > 50])- 365
dent_2019 = min(progress$Week.Ending[progress$Year == 2019 & progress$Data.Item == "CORN - PROGRESS, MEASURED IN PCT DENTED" & progress$Value > 50])- 365
mature_2019 = min(progress$Week.Ending[progress$Year == 2019 & progress$Data.Item == "CORN - PROGRESS, MEASURED IN PCT MATURE" & progress$Value > 50])- 365
harvested_2019 = min(progress$Week.Ending[progress$Year == 2019 & progress$Data.Item == "CORN, GRAIN - PROGRESS, MEASURED IN PCT HARVESTED" & progress$Value > 50])- 365



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



figName <- paste0(inFold,crop,'_2018_2019_with_crop_progress.png')
png(figName,res=300,width = 6,height = 4,units = "in",pointsize = 8)


xlims <- as.Date(c('2018-02-01','2018-12-01'))  #Always use the first of a month! (for labelling purposes)
plot(p_2017[2,]~dates,xlim=xlims,type='n',xlab='Date',ylab='EVI2',ylim=c(0.1,0.8))
#polygon(polyDates,poly17,col=transCol[1],border=polyCol[1])
polygon(polyDates,poly18,col=transCol[2],border=polyCol[2])
polygon(polyDates,poly19,col=transCol[3],border=polyCol[3])

axis.Date(1,at=seq(xlims[1], xlims[2],by="months"),
          labels=format(seq(xlims[1], xlims[2],by="months"),'%b'),
          tck=0.025)

usda_col = c('#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f','#e5c494')

abline(v=planted_2018,col=usda_col[1],lty='solid',lwd=2)
abline(v=emerged_2018,col=usda_col[2],lty='solid',lwd=2)
abline(v=silking_2018,col=usda_col[3],lty='solid',lwd=2)
abline(v=dough_2018,col=usda_col[4],lty='solid',lwd=2)
abline(v=dent_2018,col=usda_col[5],lty='solid',lwd=2)
abline(v=mature_2018,col=usda_col[6],lty='solid',lwd=2)
abline(v=harvested_2018,col=usda_col[7],lty='solid',lwd=2)

abline(v=planted_2019,col=usda_col[1],lty='dashed',lwd=2)
abline(v=emerged_2019,col=usda_col[2],lty='dashed',lwd=2)
abline(v=silking_2019,col=usda_col[3],lty='dashed',lwd=2)
abline(v=dough_2019,col=usda_col[4],lty='dashed',lwd=2)
abline(v=dent_2019,col=usda_col[5],lty='dashed',lwd=2)
abline(v=mature_2019,col=usda_col[6],lty='dashed',lwd=2)
abline(v=harvested_2019,col=usda_col[7],lty='dashed',lwd=2)

legend('topleft',legend=c('2018','2019'),
       border=c(polyCol[2:3]),fill=c(transCol[2:3]),bty='n')

#legend('topleft',legend=c('2018','2019','50% threshold dates'),
#       border=c(polyCol[2:3],NA),fill=c(transCol[2:3],NA),lty=c(NA,NA,'dashed'),bty='n')

dev.off()


