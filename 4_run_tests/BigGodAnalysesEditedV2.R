#library(plotrix)
#setwd("/Users/pesavage/Documents/Research/Oxford Seshat/Data/SCBigGodsOct2017")

rm(list = ls())
source("../project_support.R")

dir_init("./temp")

polities <- read.csv('./input/polities.csv', header=TRUE)

#New scripts for automated analysis of rates of change in social complexity pre/post moralising gods/doctrinal mode/writing

dat <- read.table("./input/PC1_traj_merged.csv", sep=",", header=TRUE) # from 3_run_pca
dat$NGA<-as.character(dat$NGA)
NGAs <- levels(polities$NGA)
NGAs <- NGAs[NGAs != "Crete"]    #### Remove new NGAs
NGAs <- NGAs[NGAs != "Galilee"]

#Overall rate (beginning to end of polity)
out <- matrix(NA, nrow=0, ncol=4)
for(i in 1:length(NGAs)){
  dt <- dat[dat$NGA == NGAs[i],]
  Earliest<-subset(dt,Time==min(Time))	
  Latest<-subset(dt,Time==max(Time))	
  MG<-subset(dt,MoralisingGods=="1") #Replace "MoralisingGods" with "DoctrinalMode" or "Writing" to do these analyses
  MGAppear<-subset(MG, Time==min(Time))
  rates<-cbind(MGAppear$NGA,(MGAppear$Mean-Earliest$Mean)/(MGAppear$Time-Earliest$Time),((Latest$Mean-MGAppear$Mean)/(Latest$Time-MGAppear$Time)),(MGAppear$End-MGAppear$Start))
  out <- rbind(out,rates)
}
colnames(out)<-c("NGA","PreRate","PostRate","MGUncertainty")
#mean(out$MGUncertainty) #Use this while replacing "MoralisingGods" as above to get uncertainty values for time-series

write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE) #Exporting/importing to force it to read as numeric (there is probably a more elegant way to do this)
out<-read.table("./temp/FullRates.csv", sep=",", header=TRUE)

out$Difference<-out[,2]-out[,3]
write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE)


#Full time windows (up to 10,000 years before and after moralizing gods)
NGAs <- levels(out$NGA)
out <- matrix(NA, nrow=0, ncol=5)
for(i in 1:length(NGAs)){
  dt <- dat[dat$NGA == NGAs[i],]
  MG<-subset(dt,MoralisingGods=="1") #Replace "MoralisingGods" with "DoctrinalMode" or "Writing" to do these analyses
  #library(dplyr)
  #MG<-as.data.frame(MG %>% group_by(PolID) %>% sample_n(size = 1)) #randomly samples so there is only one century per polity
  MGAppear<-subset(MG, Time==min(Time))
  for(j in 1: 100){
    Earliest<-subset(dt, Time==MGAppear$Time-j*100) 
    Latest<-subset(dt, Time==MGAppear$Time+j*100)
    rates<-cbind(MGAppear$NGA,ifelse(class(Earliest$Time)=="NULL","NA",(MGAppear$Mean-Earliest$Mean)/(MGAppear$Time-Earliest$Time)),ifelse(class(Latest$Time)=="NULL","NA",((Latest$Mean-MGAppear$Mean)/(Latest$Time-MGAppear$Time))),(MGAppear$End-MGAppear$Start),j*100)
    out <- rbind(out,rates)
  }
  out <- rbind(out,rates)
}
colnames(out)<-c("NGA","PreRate","PostRate","MGUncertainty","TimeWindow")

#mean(out$MGUncertainty) #Use this while replacing "MoralisingGods" as above to get uncertainty values for time-series

write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE) #Exporting/importing to force it to read as numeric (there is probably a more elegant way to do this)
out<-read.table("./temp/FullRates.csv", sep=",", header=TRUE)

out$Difference<-out[,3]-out[,2]

for(i in 2:length(out[,5])){
  out[i,7]<-out[i,5]-out[i-1,5]
}
out <-subset(out, out[,7]!=0) #getting rid of bug when the final row repeats in each NGA

write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE)

out <-subset(out, out[,5]<2050) #Change this to modify time-window restriction from 700 years pre/post moralizing gods (<750) or # out to use full time-window

write.csv(out, file="./temp/EqualRates.csv",  row.names=FALSE)

#bar chart paired
my.values<-mean(1000*out[,6],na.rm=TRUE)
err1<-1.96*std.error(1000*out[,6],na.rm=TRUE)
x <- barplot(my.values,names.arg="After - before moralizing gods",ylim = c(-1.3,1),ylab="Rate of increase in social complexity (SC/ky)")
arrows(x,my.values-err1 ,x,my.values+err1, code=3, angle=90, length=.1)
abline(h=0)

#significance test
print(t.test(out[,3], out[,2],paired=TRUE))

#histogram of differences
hist(1000*out[,6],col="gray",breaks="FD",xlim=c(-15,5),ylim=c(0,80),xlab="Rates of change in social complexity before vs. after moralizing gods (SC/ky)", ylab="Frequency")
abline(v = 0,col="black")

#Mean Pre:post ratio
print(mean(out[,2],na.rm=TRUE)/mean(out[,3],na.rm=TRUE))

###Calculate p-values etc.


data <- matrix(NA, nrow=0, ncol=5)

for(i in 1:length(NGAs)){
  dt <- out[out$NGA == NGAs[i],]
  ot <- dat[dat$NGA == NGAs[i],]
  MG<-subset(ot,MoralisingGods=="1") #Replace "MoralisingGods" with "DoctrinalMode" or "Writing" to do these analyses
  #library(dplyr)
  #MG<-as.data.frame(MG %>% group_by(PolID) %>% sample_n(size = 1)) #randomly samples so there is only one century per polity
  MGAppear<-subset(MG, Time==min(Time))
  DM<-subset(ot,DoctrinalMode=="1") #Replace "MoralisingGods" with "DoctrinalMode" or "Writing" to do these analyses
  DMAppear<-subset(DM, Time==min(Time))
  WR<-subset(ot,Writing=="1") #Replace "MoralisingGods" with "DoctrinalMode" or "Writing" to do these analyses
  WRAppear<-subset(WR, Time==min(Time))
  my.values<-c(NGAs[i],ifelse(class(mean(dt[,3],na.rm=TRUE))=="NULL","NA",mean(dt[,3],na.rm=TRUE)),ifelse(class(mean(dt[,2],na.rm=TRUE))=="NULL","NA",mean(dt[,2],na.rm=TRUE)),ifelse(class(std.error(dt[,3],na.rm=TRUE))=="NULL","NA",1.96*std.error(dt[,3],na.rm=TRUE)),ifelse(class(std.error(dt[,2],na.rm=TRUE))=="NULL","NA",1.96*std.error(dt[,2],na.rm=TRUE)))
  data <- rbind(data,my.values)
}
colnames(data)<-c("NGA","PostRate","PreRate","PostConfInt","PreConfInt")
write.csv(data, file="./temp/PrePostComparisonFull.csv",  row.names=FALSE)
data<-read.table("./temp/PrePostComparisonFull.csv", sep=",", header=TRUE)
data<-as.data.frame(data)
data$Difference<-data[,2]-data[,3]
data[,2:6]<-data[,2:6]*1000
write.csv(data, file="./temp/PrePostComparisonFull.csv",  row.names=FALSE)



#Full values for matching pre-/post-NGAs
out <-subset(out, out[,6]<1000) #Removing windows without matching pre-/post-MG rates
write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE)
out <- read.table("./temp/FullRates.csv", sep=",", header=TRUE)
NGAs <- levels(out$NGA)

data <- matrix(NA, nrow=0, ncol=15)

# error in this loop:
for(i in 1:length(NGAs)){
  if (sum(out$NGA == NGAs[i]) > 1) { # added b/c NGA = Cusco has only 1 observation, can't do t-test
    dt <- out[out$NGA == NGAs[i],]
    ot <- dat[dat$NGA == NGAs[i],]
    MG<-subset(ot,MoralisingGods=="1") #Replace "MoralisingGods" with "DoctrinalMode" or "Writing" to do these analyses
    #library(dplyr)
    #MG<-as.data.frame(MG %>% group_by(PolID) %>% sample_n(size = 1)) #randomly samples so there is only one century per polity
    MGAppear<-subset(MG, Time==min(Time))
    DM<-subset(ot,DoctrinalMode=="1") #Replace "MoralisingGods" with "DoctrinalMode" or "Writing" to do these analyses
    DMAppear<-subset(DM, Time==min(Time))
    WR<-subset(ot,Writing=="1") #Replace "MoralisingGods" with "DoctrinalMode" or "Writing" to do these analyses
    WRAppear<-subset(WR, Time==min(Time))
    my.values<-c(NGAs[i], mean(dt[,3]), mean(dt[,2]),
      1.96 * std.error(dt[,3]), 1.96 * std.error(dt[,2]),
      t.test(dt[,3],dt[,2])$p.value,
      t.test(dt[,3],dt[,2])$parameter, length(dt[,2]),
      t.test(dt[,3],dt[,2])$statistic, DMAppear$Time - MGAppear$Time,
      WRAppear$Time - MGAppear$Time,
      MGAppear$End - MGAppear$Start,
      DMAppear$End - DMAppear$Start,
      WRAppear$End - WRAppear$Start,
      MGAppear$Time
    )
    data <- rbind(data,my.values)
  }
}

colnames(data)<-c("NGA","PostRate","PreRate","PostConfInt","PreConfInt","p","df","n","t","PreMGRitual","PreMGWriting","RangeMGAppear","RangeDMAppear","RangeWRAppear","MGAppear")
write.csv(data, file="./temp/PrePostComparison.csv",  row.names=FALSE)
data<-read.table("./temp/PrePostComparison.csv", sep=",", header=TRUE)
data<-as.data.frame(data)
data$Difference<-data[,2]-data[,3]
data[,c(2:5,16)]<-data[,c(2:5,16)]*1000
write.csv(data, file="./temp/PrePostComparison.csv",  row.names=FALSE)

#Test signifiance of doctrinal ritual preceding moralizing gods
print(t.test(data$PreMGRitual))

######Normalize time-series centered around moralising god appearance
out <- matrix(NA, nrow=0, ncol=0)
for(i in 1:length(NGAs)){
  dt <- dat[dat$NGA == NGAs[i],]
  Earliest<-subset(dt,Time==min(Time))	
  Latest<-subset(dt,Time==max(Time))	
  MG<-subset(dt,MoralisingGods=="1") #Replace "MoralisingGods" with "DoctrinalMode" or "Writing" to do these analyses
  MGAppear<-subset(MG, Time==min(Time))
  dt$Time.norm<-dt$Time-MGAppear$Time
  out <- rbind(out,dt)
}
out$MGUncertainty<-out$End-out$Start
write.csv(out, file="./temp/TimeNorm.csv",  row.names=FALSE) 

#Merge Normalized times
dat <- read.table("./temp/TimeNorm.csv", sep=",", header=TRUE)
out<-unique(dat$Time.norm)

#library(plyr)
#library(dplyr) ##Bugs when I load dplyr and plyr (even when only loading dplyr after plyr)!

for(i in 1:length(NGAs)){
  dt <- dat[dat$NGA == NGAs[i],]
  out<-merge(out,dt[,c("Time.norm","Mean")],by.x="x",by.y="Time.norm",all=TRUE)
  out<-rename(out, c("Mean" = NGAs[i]))
}
MeanPCs<-out[,2:(1+length(NGAs))]
out$Mean <- apply(MeanPCs,1,mean,na.rm=TRUE)
out$Lower <- out$Mean - 1.96*apply(MeanPCs,1, std.error,na.rm=TRUE)
out$Upper <- out$Mean + 1.96*apply(MeanPCs,1, std.error,na.rm=TRUE)
write.csv(out, file="./temp/SCNorm.csv",  row.names=FALSE) 

#plot normalized times
FullImpDat<-read.csv('./temp/SCNorm.csv', header=TRUE)

####Plot 12 NGA time-series with pre- and post-moralizing god SC data
#First average normalized to moralizing gods time
out<-read.csv('./temp/SCNorm.csv', header=TRUE)
data<-read.csv('./temp/PrePostComparison.csv', header=TRUE)

#data <- FullImpDat
y <- out$Mean
x <- out$x
l <- out$Lower
u <- out$Upper
#earliest<-read.csv('EarliestAppearance.csv', header=TRUE)
#a<-earliest[,5]
#b<-earliest[,6]
#c<-earliest[,2]
#d<-earliest[,4]

#ylim <- c(min(na.rm=TRUE,y), max(na.rm=TRUE,y)) #This sets y-axis from the minimum to maximum values throughout the whole global dataset
ylim <- c(0,1) #This sets x-axis to min (0) to max (1) social complexity
#xlim <- c(min(na.rm=TRUE,x), max(na.rm=TRUE,x)) #Use this instead to set x-axis from earliest to latest time in dataset
xlim <- c(-2000,2000) #This sets x-axis to 1,000 years before and after the appearance of moralising gods
pch=19
cex=.5
xaxt='s'
yaxt='s'
xaxs="i"
yaxs="i"
ann=FALSE
#v1=min(subset(x,a==1))
#v2=min(subset(x,c==1))
linecol1<-"red"
linecol2<-"green"
linecol3<-"blue"
linecol4<-"orange"
lty1<-2
lty2<-4
lty3<-3
lty4<-2
lwd<-1.5
type="l"
h=.6

col1<-rgb(0,0,0,max=255,alpha=50)
#col2<-rgb(255,0,0,max=255,alpha=125)
#col3<-rgb(0,255,0,max=255,alpha=125)

plot(x, y, ylim=ylim, xlim=xlim, pch=pch, cex=cex, xaxt=xaxt, ann=ann, yaxt=yaxt,type=type,xaxs=xaxs,yaxs=yaxs)
panel.first = rect(c(mean(data$PreMGRitual)-(1.96*std.error(data$PreMGRitual)), 0), -1e6, c(mean(data$PreMGRitual)+ (1.96*std.error(data$PreMGRitual)), 0+mean(data$RangeMGAppear)), 1e6, col=col1, border=NA)
#panel.first = rect(c(mean(data$PreMGRitual)-(1.96*std.error(data$PreMGRitual)), mean(data$PreMGWriting)-(1.96*std.error(data$PreMGWriting)), 0), -1e6, c(mean(data$PreMGRitual)+ (1.96*std.error(data$PreMGRitual)), mean(data$PreMGWriting)+(1.96*std.error(data$PreMGWriting)), 0+mean(data$RangeMGAppear)), 1e6, col=c(col1,col2,col3), border=NA)
#This was earlier version that included writing
#polygon(c(x, rev(x)) , c(u, rev(l)) , col = 'grey' , border = NA) #out for now because of bugs
lines(x, y,type="l") 
lines(x, u,type="l",lty="dotted") 
lines(x, l,type="l",lty="dotted")
#abline(h=0.6,lty="dashed")

####Next plot each NGA individually without normalizing time
data <- read.table("./temp/TimeNorm.csv", sep=",", header=TRUE)
y <- data$Mean
x <- data$Time
l <- data$Lower
u <- data$Upper
earliest<-read.csv('./temp/PrePostComparison.csv', header=TRUE)
MG<-earliest[,15]
DM<-earliest[,10]+MG
WR<-earliest[,11]+MG
MGR<-earliest[,12]
DMR<-earliest[,13]
WRR<-earliest[,14]
#ylim <- c(min(na.rm=TRUE,y), max(na.rm=TRUE,y)) #This sets y-axis from the minimum to maximum values throughout the whole global dataset
ylim <- c(0,1) #This sets x-axis to min (0) to max (1) social complexity
#xlim <- c(min(na.rm=TRUE,x), max(na.rm=TRUE,x)) #Use this instead to set x-axis from earliest to latest time in dataset
xlim <- c(-10000,2000) #This sets x-axis to 1,000 years before and after the appearance of moralising gods
pch=19
cex=.5
xaxt='s'  
yaxt='s' 
xaxs="i"
yaxs="i"

ann=FALSE
#v1=min(subset(x,a==1))
#v2=min(subset(x,c==1))
linecol1<-"red"
linecol2<-"green"
linecol3<-"blue"
linecol4<-"orange"
lty1<-2
lty2<-4
lty3<-3
lty4<-2
lwd<-1.5
type="l"
h=0

col1<-rgb(0,0,255,max=255,alpha=125)
col2<-rgb(0,255,0,max=255,alpha=125)
col3<-rgb(255,0,0,max=255,alpha=125)

par(mfrow=c(2,6),mar=c(1.5,1.5,1.5,1.5),oma=c(0,0,0,0),mgp = c(.5, .5, 0),xpd = FALSE)

NGA<-"Upper Egypt"
plot(subset(x, data$NGA==NGA), subset(y, data$NGA==NGA),  ylim=ylim, xlim=xlim, pch=pch, cex=cex, xaxt=xaxt, ann=ann, yaxt=yaxt,type=type, xaxs=xaxs,yaxs=yaxs)
polygon(c(subset(x, data$NGA==NGA), rev(subset(x, data$NGA==NGA))) , c(subset(u, data$NGA==NGA) , rev(subset(l, data$NGA==NGA))) , col = 'grey' , border = NA)
lines(subset(x, data$NGA==NGA), subset(y, data$NGA==NGA),type="l")
abline(h=h,lty=lty1)
panel.first = rect(c(subset(DM, earliest$NGA==NGA), subset(MG, earliest$NGA==NGA)), -1e6, c(subset(DM, earliest$NGA==NGA) + subset(DMR, earliest$NGA==NGA), subset(MG, earliest$NGA==NGA) + subset(MGR, earliest$NGA==NGA)), 1e6, col=c(col1,col3), border=NA)

NGA<-"Susiana"
plot(subset(x, data$NGA==NGA), subset(y, data$NGA==NGA),  ylim=ylim, xlim=xlim, pch=pch, cex=cex, xaxt=xaxt, ann=ann, yaxt=yaxt,type=type, xaxs=xaxs,yaxs=yaxs)
polygon(c(subset(x, data$NGA==NGA), rev(subset(x, data$NGA==NGA))) , c(subset(u, data$NGA==NGA) , rev(subset(l, data$NGA==NGA))) , col = 'grey' , border = NA)
lines(subset(x, data$NGA==NGA), subset(y, data$NGA==NGA),type="l")
abline(h=h,lty=lty1)
panel.first = rect(c(subset(DM, earliest$NGA==NGA), subset(MG, earliest$NGA==NGA)), -1e6, c(subset(DM, earliest$NGA==NGA) + subset(DMR, earliest$NGA==NGA), subset(MG, earliest$NGA==NGA) + subset(MGR, earliest$NGA==NGA)), 1e6, col=c(col1,col3), border=NA)

NGA<-"Konya Plain"
plot(subset(x, data$NGA==NGA), subset(y, data$NGA==NGA),  ylim=ylim, xlim=xlim, pch=pch, cex=cex, xaxt=xaxt, ann=ann, yaxt=yaxt,type=type, xaxs=xaxs,yaxs=yaxs)
polygon(c(subset(x, data$NGA==NGA), rev(subset(x, data$NGA==NGA))) , c(subset(u, data$NGA==NGA) , rev(subset(l, data$NGA==NGA))) , col = 'grey' , border = NA)
lines(subset(x, data$NGA==NGA), subset(y, data$NGA==NGA),type="l")
abline(h=h,lty=lty1)
panel.first = rect(c(subset(DM, earliest$NGA==NGA), subset(MG, earliest$NGA==NGA)), -1e6, c(subset(DM, earliest$NGA==NGA) + subset(DMR, earliest$NGA==NGA), subset(MG, earliest$NGA==NGA) + subset(MGR, earliest$NGA==NGA)), 1e6, col=c(col1,col3), border=NA)

NGA<-"Middle Yellow River Valley"
plot(subset(x, data$NGA==NGA), subset(y, data$NGA==NGA),  ylim=ylim, xlim=xlim, pch=pch, cex=cex, xaxt=xaxt, ann=ann, yaxt=yaxt,type=type, xaxs=xaxs,yaxs=yaxs)
polygon(c(subset(x, data$NGA==NGA), rev(subset(x, data$NGA==NGA))) , c(subset(u, data$NGA==NGA) , rev(subset(l, data$NGA==NGA))) , col = 'grey' , border = NA)
lines(subset(x, data$NGA==NGA), subset(y, data$NGA==NGA),type="l")
abline(h=h,lty=lty1)
panel.first = rect(c(subset(DM, earliest$NGA==NGA), subset(MG, earliest$NGA==NGA)), -1e6, c(subset(DM, earliest$NGA==NGA) + subset(DMR, earliest$NGA==NGA), subset(MG, earliest$NGA==NGA) + subset(MGR, earliest$NGA==NGA)), 1e6, col=c(col1,col3), border=NA)

NGA<-"Kachi Plain"
plot(subset(x, data$NGA==NGA), subset(y, data$NGA==NGA),  ylim=ylim, xlim=xlim, pch=pch, cex=cex, xaxt=xaxt, ann=ann, yaxt=yaxt,type=type, xaxs=xaxs,yaxs=yaxs)
polygon(c(subset(x, data$NGA==NGA), rev(subset(x, data$NGA==NGA))) , c(subset(u, data$NGA==NGA) , rev(subset(l, data$NGA==NGA))) , col = 'grey' , border = NA)
lines(subset(x, data$NGA==NGA), subset(y, data$NGA==NGA),type="l")
abline(h=h,lty=lty1)
panel.first = rect(c(subset(DM, earliest$NGA==NGA), subset(MG, earliest$NGA==NGA)), -1e6, c(subset(DM, earliest$NGA==NGA) + subset(DMR, earliest$NGA==NGA), subset(MG, earliest$NGA==NGA) + subset(MGR, earliest$NGA==NGA)), 1e6, col=c(col1,col3), border=NA)

NGA<-"Sogdiana"
plot(subset(x, data$NGA==NGA), subset(y, data$NGA==NGA),  ylim=ylim, xlim=xlim, pch=pch, cex=cex, xaxt=xaxt, ann=ann, yaxt=yaxt,type=type, xaxs=xaxs,yaxs=yaxs)
polygon(c(subset(x, data$NGA==NGA), rev(subset(x, data$NGA==NGA))) , c(subset(u, data$NGA==NGA) , rev(subset(l, data$NGA==NGA))) , col = 'grey' , border = NA)
lines(subset(x, data$NGA==NGA), subset(y, data$NGA==NGA),type="l")
abline(h=h,lty=lty1)
panel.first = rect(c(subset(DM, earliest$NGA==NGA), subset(MG, earliest$NGA==NGA)), -1e6, c(subset(DM, earliest$NGA==NGA) + subset(DMR, earliest$NGA==NGA), subset(MG, earliest$NGA==NGA) + subset(MGR, earliest$NGA==NGA)), 1e6, col=c(col1,col3), border=NA)

NGA<-"Latium"
plot(subset(x, data$NGA==NGA), subset(y, data$NGA==NGA),  ylim=ylim, xlim=xlim, pch=pch, cex=cex, xaxt=xaxt, ann=ann, yaxt=yaxt,type=type, xaxs=xaxs,yaxs=yaxs)
polygon(c(subset(x, data$NGA==NGA), rev(subset(x, data$NGA==NGA))) , c(subset(u, data$NGA==NGA) , rev(subset(l, data$NGA==NGA))) , col = 'grey' , border = NA)
lines(subset(x, data$NGA==NGA), subset(y, data$NGA==NGA),type="l")
abline(h=h,lty=lty1)
panel.first = rect(c(subset(DM, earliest$NGA==NGA), subset(MG, earliest$NGA==NGA)), -1e6, c(subset(DM, earliest$NGA==NGA) + subset(DMR, earliest$NGA==NGA), subset(MG, earliest$NGA==NGA) + subset(MGR, earliest$NGA==NGA)), 1e6, col=c(col1,col3), border=NA)

NGA<-"Deccan"
plot(subset(x, data$NGA==NGA), subset(y, data$NGA==NGA),  ylim=ylim, xlim=xlim, pch=pch, cex=cex, xaxt=xaxt, ann=ann, yaxt=yaxt,type=type, xaxs=xaxs,yaxs=yaxs)
polygon(c(subset(x, data$NGA==NGA), rev(subset(x, data$NGA==NGA))) , c(subset(u, data$NGA==NGA) , rev(subset(l, data$NGA==NGA))) , col = 'grey' , border = NA)
lines(subset(x, data$NGA==NGA), subset(y, data$NGA==NGA),type="l")
abline(h=h,lty=lty1)
panel.first = rect(c(subset(DM, earliest$NGA==NGA), subset(MG, earliest$NGA==NGA)), -1e6, c(subset(DM, earliest$NGA==NGA) + subset(DMR, earliest$NGA==NGA), subset(MG, earliest$NGA==NGA) + subset(MGR, earliest$NGA==NGA)), 1e6, col=c(col1,col3), border=NA)

NGA<-"Paris Basin"
plot(subset(x, data$NGA==NGA), subset(y, data$NGA==NGA),  ylim=ylim, xlim=xlim, pch=pch, cex=cex, xaxt=xaxt, ann=ann, yaxt=yaxt,type=type, xaxs=xaxs,yaxs=yaxs)
polygon(c(subset(x, data$NGA==NGA), rev(subset(x, data$NGA==NGA))) , c(subset(u, data$NGA==NGA) , rev(subset(l, data$NGA==NGA))) , col = 'grey' , border = NA)
lines(subset(x, data$NGA==NGA), subset(y, data$NGA==NGA),type="l")
abline(h=h,lty=lty1)
panel.first = rect(c(subset(DM, earliest$NGA==NGA), subset(MG, earliest$NGA==NGA)), -1e6, c(subset(DM, earliest$NGA==NGA) + subset(DMR, earliest$NGA==NGA), subset(MG, earliest$NGA==NGA) + subset(MGR, earliest$NGA==NGA)), 1e6, col=c(col1,col3), border=NA)

NGA<-"Orkhon Valley"
plot(subset(x, data$NGA==NGA), subset(y, data$NGA==NGA),  ylim=ylim, xlim=xlim, pch=pch, cex=cex, xaxt=xaxt, ann=ann, yaxt=yaxt,type=type, xaxs=xaxs,yaxs=yaxs)
polygon(c(subset(x, data$NGA==NGA), rev(subset(x, data$NGA==NGA))) , c(subset(u, data$NGA==NGA) , rev(subset(l, data$NGA==NGA))) , col = 'grey' , border = NA)
lines(subset(x, data$NGA==NGA), subset(y, data$NGA==NGA),type="l")
abline(h=h,lty=lty1)
panel.first = rect(c(subset(DM, earliest$NGA==NGA), subset(MG, earliest$NGA==NGA)), -1e6, c(subset(DM, earliest$NGA==NGA) + subset(DMR, earliest$NGA==NGA), subset(MG, earliest$NGA==NGA) + subset(MGR, earliest$NGA==NGA)), 1e6, col=c(col1,col3), border=NA)

NGA<-"Kansai"
plot(subset(x, data$NGA==NGA), subset(y, data$NGA==NGA),  ylim=ylim, xlim=xlim, pch=pch, cex=cex, xaxt=xaxt, ann=ann, yaxt=yaxt,type=type, xaxs=xaxs,yaxs=yaxs)
polygon(c(subset(x, data$NGA==NGA), rev(subset(x, data$NGA==NGA))) , c(subset(u, data$NGA==NGA) , rev(subset(l, data$NGA==NGA))) , col = 'grey' , border = NA)
lines(subset(x, data$NGA==NGA), subset(y, data$NGA==NGA),type="l")
abline(h=h,lty=lty1)
panel.first = rect(c(subset(DM, earliest$NGA==NGA), subset(MG, earliest$NGA==NGA)), -1e6, c(subset(DM, earliest$NGA==NGA) + subset(DMR, earliest$NGA==NGA), subset(MG, earliest$NGA==NGA) + subset(MGR, earliest$NGA==NGA)), 1e6, col=c(col1,col3), border=NA)

NGA<-"Niger Inland Delta"
plot(subset(x, data$NGA==NGA), subset(y, data$NGA==NGA),  ylim=ylim, xlim=xlim, pch=pch, cex=cex, xaxt=xaxt, ann=ann, yaxt=yaxt,type=type, xaxs=xaxs,yaxs=yaxs)
polygon(c(subset(x, data$NGA==NGA), rev(subset(x, data$NGA==NGA))) , c(subset(u, data$NGA==NGA) , rev(subset(l, data$NGA==NGA))) , col = 'grey' , border = NA)
lines(subset(x, data$NGA==NGA), subset(y, data$NGA==NGA),type="l")
abline(h=h,lty=lty1)
panel.first = rect(c(subset(DM, earliest$NGA==NGA), subset(MG, earliest$NGA==NGA)), -1e6, c(subset(DM, earliest$NGA==NGA) + subset(DMR, earliest$NGA==NGA), subset(MG, earliest$NGA==NGA) + subset(MGR, earliest$NGA==NGA)), 1e6, col=c(col1,col3), border=NA)
