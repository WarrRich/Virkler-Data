setwd('C:/Users/rlw86/Box/Statistics-RLW86/Statistics Research/Degradation Reliability/Application')

rawData <- read.csv('VirklerData-Raw.csv')

data <- matrix(NA,ncol=3,nrow=68*9)
colnames(data) <- c('Item Number','Cycles x 1000','Length -- mm')

data[,1] <- rep(1:68,each=9)
data[,3] <- rep(c(9,11,13,17,20,26,33,39,49.8),68)
  
for (i in 1:68) {
  data[(i-1)*9+(1:9),2] <- t(3200*(rawData[i,2:10]-rawData[i,2])/(rawData[i,11]-rawData[i,2]))
}
data

plot(data[1:9,2]/10,data[1:9,3],xlim=c(0,320),ylim=c(0,50),ylab='Crack Length in mm',xlab='Thousands of Cycles')
par(new=T)
plot(data[1:9,2]/10,data[1:9,3],xlim=c(0,320),ylim=c(0,50),xlab='',ylab='',type='l',yaxt='n',xaxt='n')
for (i in 2:68) {
  par(new=T)
  plot(data[(i-1)*9+(1:9),2]/10,data[(i-1)*9+(1:9),3],xlim=c(0,320),ylim=c(0,50),xlab='',ylab='',yaxt='n',xaxt='n')
  par(new=T)
  plot(data[(i-1)*9+(1:9),2]/10,data[(i-1)*9+(1:9),3],xlim=c(0,320),ylim=c(0,50),xlab='',ylab='',type='l',yaxt='n',xaxt='n')
}


hist(data[(1:68)*9-7,2],br=30,xlim=c(0,3200))
hist(data[(1:68)*9-6,2],br=30,xlim=c(0,3200))
hist(data[(1:68)*9-5,2],br=30,xlim=c(0,3200))
hist(data[(1:68)*9-4,2],br=30,xlim=c(0,3200))
hist(data[(1:68)*9-3,2],br=30,xlim=c(0,3200))
hist(data[(1:68)*9-2,2],br=30,xlim=c(0,3200))
hist(data[(1:68)*9-1,2],br=30,xlim=c(0,3200))
hist(data[(1:68)*9-0,2],br=30,xlim=c(0,3200))

##############
# Splines

plot(data[1:9,2]/10,data[1:9,3],xlim=c(0,330),ylim=c(0,50),ylab='Crack Length in mm',xlab='Thousands of Cycles')
par(new=T)
tempFunct <- function(x) splinefun(x=data[1:9,2]/10,y=data[1:9,3])(x)
curve(tempFunct(x),xlim=c(0,330),ylim=c(0,50),xlab='',ylab='',type='l',yaxt='n',xaxt='n')
for (i in 2:68) {
  #  par(new=T)
  #  plot(data[(i-1)*9+(1:9),2]/10,data[(i-1)*9+(1:9),3],xlim=c(0,330),ylim=c(0,50),xlab='',ylab='',yaxt='n',xaxt='n')
  par(new=T)
  tempFunct <- function(x) splinefun(x=data[(i-1)*9+(1:9),2]/10,y=data[(i-1)*9+(1:9),3])(x)
  curve(tempFunct(x),xlim=c(0,330),ylim=c(0,50),xlab='',ylab='',type='l',yaxt='n',xaxt='n')
}

## Generating data with a random length

newdata <- matrix(NA,ncol=3,nrow=816)
for (i in 1:68) {
  newdata[(i-1)*12+(1:12),1] <- i
  newdata[(i-1)*12+(1:12),2] <- t(((1:12)*200))
  newdata[(i-1)*12+(1:12),3] <- t(splinefun(x=data[(i-1)*9+(1:9),2],y=data[(i-1)*9+(1:9),3])((1:12)*200))
}

# Remove some data to fit the mold of real data
#newdata <- newdata[which(newdata[,3]<49.8),]
smalldata <- NULL

for (i in 1:68) {
  indices <- which(newdata[,1]==i)
  limit <- TRUE
  for (j in indices) {
    if (limit) smalldata <- rbind(smalldata,newdata[j,])
    if (newdata[j,3] > 30) limit <- FALSE
  }
}

# Plot the data
tempFunct <- function(x) splinefun(x=data[1:9,2]/10,y=data[1:9,3])(x)
curve(tempFunct(x),xlim=c(0,330),ylim=c(0,50),ylab='Crack Length in mm',xlab='Thousands of Cycles',col='gray')
for (i in 2:68) {
  par(new=T)
  tempFunct <- function(x) splinefun(x=data[(i-1)*9+(1:9),2]/10,y=data[(i-1)*9+(1:9),3])(x)
  curve(tempFunct(x),xlim=c(0,330),ylim=c(0,50),xlab='',ylab='',type='l',yaxt='n',xaxt='n',col='gray')
}
par(new=T)
plot(smalldata[,2]/10,smalldata[,3],xlim=c(0,330),ylim=c(0,50),xlab='',ylab='',yaxt='n',xaxt='n')

abline(h=30)

######
# Saving the data

smalldata[,2] <- smalldata[,2]/10

write.csv(smalldata,file = "VirklerData.csv",row.names=F)





######## Extra

## Looking at the means
means <- matrix(NA,nrow=68,ncol=42)
for (i in 1:68) {
  means[i,] <- splinefun(x=data[(i-1)*9+(1:9),3],y=data[(i-1)*9+(1:9),2]/10)(9:50)
}
meanCycles <- apply(means,2,mean)
sdCycles <- apply(means,2,sd)
par(new=T);plot(meanCycles,9:50,type='l',col=2,xlim=c(0,330),ylim=c(0,50),xlab='',ylab='',yaxt='n',xaxt='n')
par(new=T);plot(meanCycles-1.96*sdCycles/sqrt(68),9:50,type='l',col=3,xlim=c(0,330),ylim=c(0,50),xlab='',ylab='',yaxt='n',xaxt='n')

