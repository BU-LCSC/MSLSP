library(signal)
i <- 4
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
png(filename='test_filtering.png',width=7.5,height=8,units='in',res=150)
par(mfrow=c(3,1),oma=c(1,1,1,2),mar=c(3,3,1,1))
plot(pixTab$filled_VI~pixTab$dates,ylim=c(0.2,0.9))
points(pixTab$dates[ind],xSmooth[ind],col='blue',pch=19,cex=1.5)
plot(pixTab$filled_VI~pixTab$dates,ylim=c(0.2,0.9))
points(pixTab$dates[ind],xInt[ind],col='red',pch=19,cex=1.5)
plot(pixTab$filled_VI~pixTab$dates,ylim=c(0.2,0.9))
points(pixTab$dates[ind],lint[ind],col='forestgreen',pch=19,cex=1.5)
dev.off()

plot(xSmooth,xInt)
abline(0,1)
summary(lm(xInt~xSmooth))


