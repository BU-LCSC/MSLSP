library(signal)
i <- 1
for(i in 1:100){
  pixTab <- theTable1[[i]]
  plot(pixTab$filled_VI~pixTab$dates,pch=21,bg='#2c7fb8',main=i)
}

ind <- !is.na(pixTab$original_VI)

##############
# Spline
for_a_minute_1 <- function() {
  for(i in 1:100000){
    spl <- smooth.spline(pixTab$dates[ind], pixTab$original_VI[ind], spar=0.55)
    xSmooth <- predict(spl, as.numeric(pixTab$dates))$y  
  }
}

start_time <- Sys.time()
for_a_minute_1()
end_time <- Sys.time()
end_time - start_time

##############
# Linear interpolate (and filtering)
for_a_minute_2 <- function() {
  for(i in 1:100000){
    lint <- approx(pixTab$dates[ind],pixTab$original_VI[ind],xout=as.numeric(pixTab$dates),rule=2)$y
    xInt <- sgolayfilt(lint,p=3)
  }
}

start_time <- Sys.time()
for_a_minute_2()
end_time <- Sys.time()
end_time - start_time


setwd('/projectnb/modislc/users/mkmoon/MuSLI/V1_0/figures/')
png(filename='test_filtering_1.png',width=8,height=6.5,units='in',res=150)

d678 <- which(substr(pixTab$dates,1,4)=='2016'|substr(pixTab$dates,1,4)=='2017'|substr(pixTab$dates,1,4)=='2018')
par(mfrow=c(3,1),oma=c(1,1,1,2),mar=c(3,3,1,1))
plot(pixTab$filled_VI[d678]~pixTab$dates[d678],ylim=c(0.2,0.8),pch=19,col='red',cex=0,cex.axis=1.5)
points(pixTab$original_VI[ind]~pixTab$dates[ind],pch=21,bg='#2c7fb8',cex=1.5)
plot(pixTab$filled_VI[d678]~pixTab$dates[d678],ylim=c(0.2,0.8),pch=19,col='red',cex=0,cex.axis=1.5)
points(pixTab$dates[ind],lint[ind],bg='#2c7fb8',pch=21,cex=1.5)
points(pixTab$dates[ind],xSmooth[ind],bg='red',pch=21,cex=1.5)

plot(pixTab$filled_VI[d678]~pixTab$dates[d678],ylim=c(0.2,0.8),pch=19,col='red',cex=0,cex.axis=1.5)
# points(pixTab$original_VI~pixTab$dates,ylim=c(0.2,0.9),pch=19,cex=1.5)

d2016 <- (substr(pixTab$y2016$filled_dates,1,4)=='2016')
d2017 <- (substr(pixTab$y2017$filled_dates,1,4)=='2017')
d2018 <- (substr(pixTab$y2018$filled_dates,1,4)=='2018')
points(pixTab$y2016$filled_vi[d2016 & pixTab$y2016$filled_weigth < 1] ~ pixTab$y2016$filled_dates[d2016 & pixTab$y2016$filled_weigth < 1],pch=4,cex=pixTab$y2016$filled_weigth[d2016 & pixTab$y2016$filled_weigth < 1]*5)
points(pixTab$y2017$filled_vi[d2017 & pixTab$y2017$filled_weigth < 1] ~ pixTab$y2017$filled_dates[d2017 & pixTab$y2017$filled_weigth < 1],pch=4,cex=pixTab$y2017$filled_weigth[d2017 & pixTab$y2017$filled_weigth < 1]*5)
points(pixTab$y2018$filled_vi[d2018 & pixTab$y2018$filled_weigth < 1] ~ pixTab$y2018$filled_dates[d2018 & pixTab$y2018$filled_weigth < 1],pch=4,cex=pixTab$y2018$filled_weigth[d2018 & pixTab$y2018$filled_weigth < 1]*5)

d2016 <- (substr(pixTab$y2016$smoothed_dates,1,4)=='2016')
d2017 <- (substr(pixTab$y2017$smoothed_dates,1,4)=='2017')
d2018 <- (substr(pixTab$y2018$smoothed_dates,1,4)=='2018')
lines(pixTab$y2016$smoothed_vi[d2016]~pixTab$y2016$smoothed_dates[d2016],ylim=c(0.2,0.9),lwd=2,col='#31a354')
lines(pixTab$y2017$smoothed_vi[d2017]~pixTab$y2017$smoothed_dates[d2017],ylim=c(0.2,0.9),lwd=2,col='#31a354')
lines(pixTab$y2018$smoothed_vi[d2018]~pixTab$y2018$smoothed_dates[d2018],ylim=c(0.2,0.9),lwd=2,col='#31a354')

points(pixTab$filled_VI[pixTab$snowPix]~pixTab$dates[pixTab$snowPix],ylim=c(0.2,0.8),cex=2,pch=21,bg='black')
points(pixTab$original_VI~pixTab$dates,ylim=c(0.2,0.9),pch=21,bg='#2c7fb8',cex=1.5)

dev.off()


pixTab$y2016$filled_day
sum(diff(pixTab$dates[ind])>20)
plot(diff(pixTab$dates[d678]))

#Plot only the observations from other years. This will also plot snow from current year, but we'll cover those points up later
if (showFilledData) {points(subTab$filled_vi[subTab$filled_weigth < 1] ~ subTab$filled_dates[subTab$filled_weigth < 1],pch=16,cex=0.6)}

#Plot the original time-series before despiking  
if (showDespiked) {points(pixTab$original_VI~pixTab$dates)}

if (showObservations) {points(pixTab$filled_VI~pixTab$dates,pch=21,bg='#2c7fb8')}
if (showSnow) {points(pixTab$filled_VI[pixTab$snowPix]~pixTab$dates[pixTab$snowPix],pch=21,bg='black')}




plot(xSmooth,xInt)
abline(0,1)
summary(lm(xInt~xSmooth))


