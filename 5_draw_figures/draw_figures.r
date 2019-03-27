
rm(list = ls())
source("../project_support.r")

dir_init("./temp")


# the only 12 NGAs with a pre-post comparison
NGAs <-  c("Deccan", "Kachi Plain", "Kansai", "Konya Plain", 
  "Latium", "Middle Yellow River Valley", "Niger Inland Delta",
  "Orkhon Valley", "Paris Basin", "Sogdiana", 
  "Susiana", "Upper Egypt")

dat <- read.table("./input/PC1_traj_merged.csv", sep=",", header=TRUE)
dat$NGA<-as.character(dat$NGA)

TimeNorm <- matrix(NA, nrow=0, ncol=0)
for(i in 1:length(NGAs)){
  dt <- dat[dat$NGA == NGAs[i],]
  Earliest<-subset(dt,Time==min(Time))  
  Latest<-subset(dt,Time==max(Time))  
  MG<-subset(dt,MoralisingGods=="1") #Replace "MoralisingGods" with "DoctrinalMode" or "Writing" to do these analyses
  MGAppear<-subset(MG, Time==min(Time))
  dt$Time.norm<-dt$Time-MGAppear$Time
  TimeNorm <- rbind(TimeNorm,dt)
}
TimeNorm$MGUncertainty<-TimeNorm$End-TimeNorm$Start
write.csv(TimeNorm, file="./temp/TimeNorm.csv",  row.names=FALSE) 

all(dim(TimeNorm) == c(709, 51))



####plot each NGA individually without normalizing time

png("./temp/twelveNGA_arrival_ritual.png", res = 300, height = 3, width = 10, units = "in")

data <- read.table("./temp/TimeNorm.csv", sep=",", header=TRUE)
y <- data$Mean
x <- data$Time
l <- data$Lower
u <- data$Upper
earliest<-read.csv('./input/PrePostComparison.csv', header=TRUE)
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

dev.off()









#Merge Normalized times

TimeNorm <- read.table("./temp/TimeNorm.csv", sep=",", header=TRUE)
SCNorm<-unique(TimeNorm$Time.norm)

for(i in 1:length(NGAs)){
  dt <- TimeNorm[TimeNorm$NGA == NGAs[i],]
  SCNorm<-merge(SCNorm,dt[,c("Time.norm","Mean")],by.x="x",by.y="Time.norm",all=TRUE)
  SCNorm<-plyr::rename(SCNorm, c("Mean" = NGAs[i]))
}
MeanPCs<-SCNorm[,2:(1+length(NGAs))]
SCNorm$Mean <- apply(MeanPCs,1,mean,na.rm=TRUE)
SCNorm$Lower <- SCNorm$Mean - 1.96*apply(MeanPCs,1, std.error,na.rm=TRUE)
SCNorm$Upper <- SCNorm$Mean + 1.96*apply(MeanPCs,1, std.error,na.rm=TRUE)
write.csv(SCNorm, file="./temp/SCNorm.csv",  row.names=FALSE) 

all(dim(SCNorm) == c(127, 16))


# plot main averaged figure

SCNorm<-read.csv('./temp/SCNorm.csv', header=TRUE)

data<-read.csv('./input/PrePostComparison.csv', header=TRUE)

png("./temp/twelveNGA_arrival_ritual_averaged.png", res = 300, height = 5, width = 5, units = "in")

col1<-rgb(0,0,0,max=255,alpha=50)

plot(SCNorm$x, SCNorm$Mean, type = "n", ylim = c(0, 1), xlim = c(-2000, 2000), ann = FALSE,
  xaxs = "i", yaxs = "i")

# Ritual period
rect(
  mean(data$PreMGRitual)- 1.96 * std.error(data$PreMGRitual), 
  0,
  mean(data$PreMGRitual) + 1.96 * std.error(data$PreMGRitual),
  1,
  border = NA, col = col1
)

# MG appear period
rect(0, 0, 0 + mean(data$RangeMGAppear), 1,
  border = NA, col = col1)
# why would you do that? why not center the range at 0?

lines(SCNorm$x, SCNorm$Mean, type="l") 
lines(SCNorm$x, SCNorm$Upper, type="l",lty="dotted") 
lines(SCNorm$x, SCNorm$Lower, type="l",lty="dotted")

dev.off()





png("./temp/twelveNGA_arrival_ritual_writing_averaged.png", res = 300, height = 5, width = 5, units = "in")

col1<-rgb(0,0,0,max=255,alpha=50)
col2<-rgb(255,0,0,max=255,alpha=125)
col3<-rgb(0,255,0,max=255,alpha=125)

plot(SCNorm$x, SCNorm$Mean, type = "n", ylim = c(0, 1), xlim = c(-2000, 2000), ann = FALSE,
  xaxs = "i", yaxs = "i")

# Ritual period
rect(
  mean(data$PreMGRitual)- 1.96 * std.error(data$PreMGRitual), 
  0,
  mean(data$PreMGRitual) + 1.96 * std.error(data$PreMGRitual),
  1,
  border = NA, col = col1
)

# MG appear period
rect(0, 0, 0 + mean(data$RangeMGAppear), 1,
  border = NA, col = col2)

# writing appear period
rect(
  mean(data$PreMGWriting) - 1.96 * std.error(data$PreMGWriting),
  0,
  mean(data$PreMGWriting) + 1.96 * std.error(data$PreMGWriting),
  1,
  col=col3, border=NA
)

lines(SCNorm$x, SCNorm$Mean, type="l") 
lines(SCNorm$x, SCNorm$Upper, type="l",lty="dotted") 
lines(SCNorm$x, SCNorm$Lower, type="l",lty="dotted")

dev.off()




png("./temp/twelveNGA_arrival_ritual_writing_averaged.png", res = 300, height = 5, width = 5, units = "in")

y <- SCNorm$Mean
x <- SCNorm$x
l <- SCNorm$Lower
u <- SCNorm$Upper

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
col2<-rgb(255,0,0,max=255,alpha=125)
col3<-rgb(0,255,0,max=255,alpha=125)

plot(x, y, ylim=ylim, xlim=xlim, pch=pch, cex=cex, xaxt=xaxt, ann=ann, yaxt=yaxt,type=type,xaxs=xaxs,yaxs=yaxs)
panel.first = rect(
  c(mean(data$PreMGRitual)-(1.96*std.error(data$PreMGRitual)), mean(data$PreMGWriting)-(1.96*std.error(data$PreMGWriting)), 0),
  -1e6,
  c(mean(data$PreMGRitual)+ (1.96*std.error(data$PreMGRitual)), mean(data$PreMGWriting)+(1.96*std.error(data$PreMGWriting)), 0+mean(data$RangeMGAppear)),
  1e6, col=c(col1,col2,col3), border=NA
)
#This was earlier version that included writing
lines(x, y,type="l") 
lines(x, u,type="l",lty="dotted") 
lines(x, l,type="l",lty="dotted")

dev.off()

print("publication figures created")


dir_init("./output")

file.copy("./temp/SCNorm.csv", "./output")
