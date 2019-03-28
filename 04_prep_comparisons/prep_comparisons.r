
rm(list = ls())
source("../project_support.r")

dir_init("./temp")


# Create a pre/post comparison table

polities <- read.csv('./input/polities.csv', header=TRUE)

dat <- read.table("./input/PC1_traj_merged.csv", sep=",", header=TRUE)

all(dim(dat) == c(864, 49))

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

write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE)
#Exporting/importing to force it to read as numeric (there is probably a more elegant way to do this)
out<-read.table("./temp/FullRates.csv", sep=",", header=TRUE)

out$Difference<-out[,2]-out[,3]
write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE)

expect_equal(dim(out), c(20, 5))

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

write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE)
#Exporting/importing to force it to read as numeric (there is probably a more elegant way to do this)
out<-read.table("./temp/FullRates.csv", sep=",", header=TRUE)

dim(out)

out$Difference<-out[,3]-out[,2]

for(i in 2:length(out[,5])){
  out[i,7]<-out[i,5]-out[i-1,5]
}
out <-subset(out, out[,7]!=0) #getting rid of bug when the final row repeats in each NGA

write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE)

all(dim(out) == c(1999, 7))

out <-subset(out, out[,5]<2050)
#Change this to modify time-window restriction from 700 years pre/post moralizing gods (<750) or # out to use full time-window
dim(out)

write.csv(out, file="./temp/EqualRates.csv",  row.names=FALSE)

expect_equal(dim(out), c(399, 7))

#bar chart paired
png("./temp/rate_of_increase.png", res = 300, units = "in", height = 5, width = 8)
my.values<-mean(1000*out[,6],na.rm=TRUE)
err1<-1.96*std.error(1000*out[,6],na.rm=TRUE)
x <- barplot(my.values,names.arg="After - before moralizing gods",ylim = c(-1.3,1),ylab="Rate of increase in social complexity (SC/ky)")
arrows(x,my.values-err1 ,x,my.values+err1, code=3, angle=90, length=.1)
abline(h=0)
dev.off()

#histogram of differences
png("./temp/rates_of_change.png", res = 300, units = "in", height = 5, width = 8)
hist(1000*out[,6],col="gray",breaks="FD",xlim=c(-15,5),ylim=c(0,80),xlab="Rates of change in social complexity before vs. after moralizing gods (SC/ky)", ylab="Frequency")
abline(v = 0,col="black")

dev.off()

#significance test
print(t.test(out[,3], out[,2],paired=TRUE))

#Mean Pre:post ratio
print(mean(out[,2],na.rm=TRUE)/mean(out[,3],na.rm=TRUE))

#Full values for matching pre-/post-NGAs
out <- read.csv("./temp/EqualRates.csv", stringsAsFactors = FALSE)
out <-subset(out, out[,6]<1000) #Removing windows without matching pre-/post-MG rates
write.csv(out, file="./temp/EqualRates.csv",  row.names=FALSE)



# Calculate PrePostComparison

out <- read.table("./temp/EqualRates.csv", sep=",", header=TRUE)
NGAs <- levels(out$NGA)

all(dim(out) == c(200, 7))

data <- matrix(NA, nrow=0, ncol=15)

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
  my.values<-c(NGAs[i],mean(dt[,3]),mean(dt[,2]),1.96*std.error(dt[,3]),1.96*std.error(dt[,2]),t.test(dt[,3],dt[,2])$p.value,t.test(dt[,3],dt[,2])$parameter,length(dt[,2]),t.test(dt[,3],dt[,2])$statistic, DMAppear$Time-MGAppear$Time,WRAppear$Time-MGAppear$Time,MGAppear$End-MGAppear$Start,DMAppear$End-DMAppear$Start,WRAppear$End-WRAppear$Start, MGAppear$Time)
  data <- rbind(data,my.values)
}
colnames(data)<-c("NGA","PostRate","PreRate","PostConfInt","PreConfInt","p","df","n","t","PreMGRitual","PreMGWriting","RangeMGAppear","RangeDMAppear","RangeWRAppear","MGAppear")
write.csv(data, file="./temp/PrePostComparison.csv",  row.names=FALSE)
data<-read.table("./temp/PrePostComparison.csv", sep=",", header=TRUE)
data<-as.data.frame(data)
data$Difference<-data[,2]-data[,3]
data[,c(2:5,16)]<-data[,c(2:5,16)]*1000
write.csv(data, file="./temp/PrePostComparison.csv", row.names=FALSE)

expect_equal(dim(data), c(12, 16))

NGA_column <- c(
  "Deccan",                    "Kachi Plain",              
  "Kansai",                    "Konya Plain",              
  "Latium",                    "Middle Yellow River Valley",
  "Niger Inland Delta",        "Orkhon Valley",            
  "Paris Basin",               "Sogdiana",                 
  "Susiana",                   "Upper Egypt"
)

expect_equal(as.character(data$NGA), NGA_column)

#Test signifiance of doctrinal ritual preceding moralizing gods
print(t.test(data$PreMGRitual))

dir_init("./output")

file.copy("./temp/PrePostComparison.csv", "./output")

print("comparison table PrePostComparison.csv created")

