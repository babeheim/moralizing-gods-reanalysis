
rm(list = ls())
source("../project_support.r")

dir_init("./temp")


# Create a pre/post comparison table

polities <- read.csv('./input/polities.csv', header=TRUE)

#New scripts for automated analysis of rates of change in social complexity pre/post moralising gods/doctrinal mode/writing

dat <- read.table("./input/PC1_traj_merged.csv", sep=",", header=TRUE)

all(dim(dat) == c(864, 49))

dat$NGA<-as.character(dat$NGA)
# convert NGA to character from factor
NGAs <- levels(polities$NGA)
# extract NGAs
NGAs <- NGAs[NGAs != "Crete"]    #### Remove new NGAs
NGAs <- NGAs[NGAs != "Galilee"]

# Calculate the rates of appearance of analysis variable
# for main analysis this is Moralising Gods
#Overall rate (beginning to end of polity)
out <- matrix(NA, nrow=0, ncol=4)
for(i in 1:length(NGAs)){
  # split data by NGA
  dt <- dat[dat$NGA == NGAs[i],]
  # extract earliest time point in NGA
  Earliest<-subset(dt,Time==min(Time))  
  # extract latest time point in NGA
  Latest<-subset(dt,Time==max(Time))
  # extract all data with present analysis variable
  # for the main analysis this is Moralising Gods, but for confirmatory/ additional analysis
  # Replace "MoralisingGods" with "DoctrinalMode" or "Writing" to do these analyses
  MG<-subset(dt,MoralisingGods=="1")
  # extract the first time point the analysis variable appears at
  MGAppear<-subset(MG, Time==min(Time))
  # calculate rates
  # NGA
  # PreRate = Difference in mean between analysis variable appearance and first time point of NGA divided by the difference in centuries
  # PostRate = Difference in mean between analysis variable appearance and last time point of NGA divided by the difference in centuries
  # MGUncertainty = Time range (years) of polity with the first appearance of analysis variable
  rates<-cbind(MGAppear$NGA,(MGAppear$Mean-Earliest$Mean)/(MGAppear$Time-Earliest$Time),((Latest$Mean-MGAppear$Mean)/(Latest$Time-MGAppear$Time)),(MGAppear$End-MGAppear$Start))
  out <- rbind(out,rates)
}
colnames(out)<-c("NGA","PreRate","PostRate","MGUncertainty")
#mean(out$MGUncertainty) #Use this while replacing "MoralisingGods" as above to get uncertainty values for time-series


#Exporting/importing to force it to read as numeric (there is probably a more elegant way to do this)
write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE) 
out<-read.table("./temp/FullRates.csv", sep=",", header=TRUE)
# calculate the difference between the Pre- and Post- Rates
out$Difference<-out[,2]-out[,3]
write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE)

expect_equal(dim(out), c(20, 5))

# Full time windows (up to 10,000 years before and after moralizing gods)
# extract NGAs
NGAs <- levels(out$NGA)
# create empty data frame
out <- matrix(NA, nrow=0, ncol=5)
for(i in 1:length(NGAs)){
  # split data by NGA
  dt <- dat[dat$NGA == NGAs[i],]
  # extract all data with present analysis variable
  # for the main analysis this is Moralising Gods
  # Replace "MoralisingGods" with "DoctrinalMode" or "Writing" to do these analyses
  MG<-subset(dt,MoralisingGods=="1")
  #library(dplyr)
  #MG<-as.data.frame(MG %>% group_by(PolID) %>% sample_n(size = 1)) #randomly samples so there is only one century per polity
  # extract the first time point the analysis variable appears at
  MGAppear<-subset(MG, Time==min(Time))
  for(j in 1: 100){
    # extract the first time point the analysis variable appears at
    Earliest<-subset(dt, Time==MGAppear$Time-j*100) 
    # extract the last time point the analysis variable appears at
    Latest<-subset(dt, Time==MGAppear$Time+j*100)
    # calculate rates
    # NGA
    # PreRate = if there is no appearance of analysis variable in earliest time == NA, else rate =  Difference in mean between analysis variable appearance and first time point of NGA divided by the difference in centuries
    # PostRate =  if there is no appearance of analysis variable in latest time == NA, else rate = Difference in mean between analysis variable appearance and last time point of NGA divided by the difference in centuries
    # MGUncertainty = Time range (years) of polity with the first appearance of analysis variable
    # TimeWindow = Century
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
# calculate the difference between the Pre- and Post- Rates
out$Difference<-out[,3]-out[,2]

# calculate the difference between the previous time window and the current time window
# this marks the time window 100 with -9900
for(i in 2:length(out[,5])){
  out[i,7]<-out[i,5]-out[i-1,5]
}
#getting rid of bug when the final row repeats in each NGA
out <-subset(out, out[,7]!=0) 

write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE)

all(dim(out) == c(1999, 7))

# Change this to modify time-window restriction from 700 years pre/post moralizing gods (<750) or # out to use full time-window
out <-subset(out, out[,5]<2050)
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


# Figure 2b Nature
#histogram of differences
png("./temp/rates_of_change.png", res = 300, units = "in", height = 5, width = 8)
hist(1000*out[,6],col="gray",breaks="FD",xlim=c(-15,5),ylim=c(0,80),xlab="Rates of change in social complexity before vs. after moralizing gods (SC/ky)", ylab="Frequency")
abline(v = 0,col="black")

dev.off()

#significance test
print(t.test(out[,3], out[,2],paired=TRUE))

# Calculate the Mean Pre:post ratio
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
  # split rates data by NGA
  dt <- out[out$NGA == NGAs[i],]
  # split imputed data by NGA
  ot <- dat[dat$NGA == NGAs[i],]
  # extract all data with present analysis variable
  MG<-subset(ot,MoralisingGods=="1") #Replace "MoralisingGods" with "DoctrinalMode" or "Writing" to do these analyses
  #library(dplyr)
  #MG<-as.data.frame(MG %>% group_by(PolID) %>% sample_n(size = 1)) #randomly samples so there is only one century per polity
  MGAppear<-subset(MG, Time==min(Time))
  # extract DoctrinalMode == present
  DM<-subset(ot,DoctrinalMode=="1") #Replace "MoralisingGods" with "DoctrinalMode" or "Writing" to do these analyses
  # extract the first time point the DoctrinalMode appears at
  DMAppear<-subset(DM, Time==min(Time))
  # extract writing == present
  WR<-subset(ot,Writing=="1") #Replace "MoralisingGods" with "DoctrinalMode" or "Writing" to do these analyses
  # extract the first time point the WRAppear appears at
  WRAppear<-subset(WR, Time==min(Time))
  # Extract rates, confidence intervals, p, df, n, t, etc. 
  # NGA
  # PostRate = mean(PostRate)
  # PreRate = mean(PreRate)
  # PostConfInt = 1.96*std.error(PostRate)
  # PreConfInt = 1.96*std.error(PreRate)
  # t-test between pre- and post-rates
  # p = p-value
  # df
  # n = number of rows/variables
  # t = t-value
  # PreMGRitual = Difference in time between MG and DoctrinalMode appearing
  # PreMGWriting = Difference in time between writing and MG appearing
  # RangeMGAppear = Range of MG appearing
  # RangeDMAppear = Range of DoctrinalMode appearing
  # RangeWRAppear = Range of writing appearing
  # MGAppear = century MG appear
  my.values<-c(
    NGAs[i],
    mean(dt[,3]),
    mean(dt[,2]),
    1.96 * std.error(dt[,3]),
    1.96 * std.error(dt[,2]),
    t.test(dt[,3],dt[,2])$p.value,
    t.test(dt[,3],dt[,2])$parameter,
    length(dt[,2]),
    t.test(dt[,3],dt[,2])$statistic,
    DMAppear$Time - MGAppear$Time,
    WRAppear$Time - MGAppear$Time,
    MGAppear$End - MGAppear$Start,
    DMAppear$End - DMAppear$Start,
    WRAppear$End - WRAppear$Start,
    MGAppear$Time
  )
  data <- rbind(data,my.values)
}
colnames(data)<-c("NGA","PostRate","PreRate","PostConfInt","PreConfInt","p","df","n","t","PreMGRitual","PreMGWriting","RangeMGAppear","RangeDMAppear","RangeWRAppear","MGAppear")
write.csv(data, file="./temp/PrePostComparison.csv",  row.names=FALSE)
data<-read.table("./temp/PrePostComparison.csv", sep=",", header=TRUE)
data<-as.data.frame(data)
# calculate the difference between Pre- and Post-rates
data$Difference<-data[,2]-data[,3]
# multiply PostRate, PreRate, PostConfInt, PreConfInt, Difference by 1000
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

# Our data show that doctrinal rituals standardized by routinization (that is, those performed weekly or daily) or institutionalized policing (religions with multiple hierarchical levels) significantly predate moralizing gods, by an average of 1,100 years (t = 2.8, d.f. = 11, P = 0.018; Fig. 2a).
#Test signifiance of doctrinal ritual preceding moralizing gods
print(t.test(data$PreMGRitual))

dir_init("./output")

file.copy("./temp/PrePostComparison.csv", "./output")

print("comparison table PrePostComparison.csv created")

