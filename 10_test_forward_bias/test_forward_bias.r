

# Note that comments in the text starting with #??OUR_COMMENT are our new comments.
# We left all the all comments in the text as well. 

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
### 1. Original Whitehouse et al. code ####

#### 1.1 Load data ####

# This section loads data on polities, their social complexity (SC),
# and the presence moralizing gods (MG).

{
  rm(list = ls())
  
  source("../project_support.r")
  
  dir_init("./temp")
  polities <- read.csv("./input/polities.csv", header = TRUE, stringsAsFactors = TRUE)
  
  #New scripts for automated analysis of rates of change in social complexity pre/post
  # moralising gods/doctrinal mode/writing
  
  #dat <- read.table("PC1_traj_merged.csv", sep=",", header=TRUE, stringsAsFactors = TRUE) #?? everything is here
  dat <- read.csv("./input/PC1_traj_merged.csv", header = TRUE, stringsAsFactors = TRUE)
  
  dat$NGA<-as.character(dat$NGA)
  NGAs <- levels(polities$NGA)
  NGAs <- NGAs[NGAs != "Crete"]    #### Remove new NGAs
  NGAs <- NGAs[NGAs != "Galilee"]
  
 
}

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#### 1.2 Pre-Post MG ####

#Overall rate (beginning to end of polity)
out <- matrix(NA, nrow=0, ncol=4)
for(i in 1:length(NGAs)){
  dt <- dat[dat$NGA == NGAs[i],]
  Earliest<-subset(dt,Time==min(Time))	
  Latest<-subset(dt,Time==max(Time))	
  MG<-subset(dt,MoralisingGods=="1") #Replace "MoralisingGods" with "DoctrinalMode" or
  #"Writing" to do these analyses
  MGAppear<-subset(MG, Time==min(Time))
  rates<-cbind(MGAppear$NGA,(MGAppear$Mean-Earliest$Mean)/(MGAppear$Time-Earliest$Time),
               ((Latest$Mean-MGAppear$Mean)/(Latest$Time-MGAppear$Time)),
               (MGAppear$End-MGAppear$Start))
  out <- rbind(out,rates)
}
colnames(out)<-c("NGA","PreRate","PostRate","MGUncertainty")
#mean(out$MGUncertainty) #Use this while replacing "MoralisingGods" as above to get 
#uncertainty values for time-series

write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE) #Exporting/importing to force it to read
#as numeric (there is probably a more elegant way to do this)

out<-read.table("./temp/FullRates.csv", sep=",", header=TRUE, stringsAsFactors = TRUE)

out$Difference<-out[,2]-out[,3]
write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#### 1.3 Pre-Post MG by century ####

#Full time windows (up to 10,000 years before and after moralizing gods)
NGAs <- levels(out$NGA)
out <- matrix(NA, nrow=0, ncol=5)
for(i in 1:length(NGAs)){
  dt <- dat[dat$NGA == NGAs[i],]
  MG<-subset(dt,MoralisingGods=="1") #Replace "MoralisingGods" with "DoctrinalMode" 
  #or "Writing" to do these analyses
  #library(dplyr)
  #MG<-as.data.frame(MG %>% group_by(PolID) %>% sample_n(size = 1)) 
  #randomly samples so there is only one century per polity
  MGAppear<-subset(MG, Time==min(Time))
  for(j in 1: 100){
    Earliest<-subset(dt, Time==MGAppear$Time-j*100) ##?OUR_COMMENT:: Potentially confusing label,
    # it's not earliest, this is just the number of centuries before MG
    Latest<-subset(dt, Time==MGAppear$Time+j*100) ##?OUR_COMMENT:: Same as above
    rates<-cbind(MGAppear$NGA,ifelse(class(Earliest$Time)=="NULL","NA",
                                     (MGAppear$Mean-Earliest$Mean)/(MGAppear$Time-Earliest$Time)),
                 ifelse(class(Latest$Time)=="NULL","NA",
                        ((Latest$Mean-MGAppear$Mean)/(Latest$Time-MGAppear$Time))),
                 (MGAppear$End-MGAppear$Start),j*100)
    out <- rbind(out,rates)
  }
  out <- rbind(out,rates)
}
colnames(out)<-c("NGA","PreRate","PostRate","MGUncertainty","TimeWindow")

#mean(out$MGUncertainty) #Use this while replacing "MoralisingGods" as above to get uncertainty
# values for time-series

write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE) #Exporting/importing to force it to
# read as numeric (there is probably a more elegant way to do this)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#### 1.4 Testing Pre-Post MG ####

#### 1.4.1 2000 years ####
out<-read.table("./temp/FullRates.csv", sep=",", header=TRUE, stringsAsFactors = TRUE)

out$Difference<-out[,3]-out[,2]

##?OUR_COMMENT:: Let's save this variable, we will need it later for plotting the histogram
# differences
d0 <- out$Difference*1000 

##?OUR_COMMENT:: It is strange that this loop starts with "2:length" rather than "1:length".
  # It practically deletes the first row of data from the data set and there is no obvious
  # rationale for this step. Probaly, this is to avoid "NA" results for i == 1, which would give
  # out[1,5]-out[0,5] and produce an error. Note that this is random and does not depend on
  # a specfic NGA. Likewise, it does not delete the first row for each NGA, only for the one that
  # happens to be first in the data set.
  # We will keep it here to reproduce the original analysis but fix
  # it in our subsequent analyses.

for(i in 2:length(out[,5])){
  out[i,7]<-out[i,5]-out[i-1,5]
}
out <-subset(out, out[,7]!=0) #getting rid of bug when the final row repeats in each NGA

write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE)

out <-subset(out, out[,5]<2050) #Change this to modify time-window restriction from 700 years
# pre/post moralizing gods (<750) or # out to use full time-window

write.csv(out, file="./temp/EqualRates.csv",  row.names=FALSE)

##?OUR_COMMENT:: Here is another problem. Whitehouse et al. claim that all 12 key NGAs have data
  # for the +/- 700 years analysis. Yet, looking into EqualRates.csv revealed that Deccan has
  # 'NA' for the Post-MG rate of SC change at TimeWindow = 700. Also, Paris Basin has Post-MG
  # rate of SC change 'NA' at TimeWindow = 500. Inspection of the main source file 'PC1_traj_merged.csv'
  # explained the problem - Deccan and Paris Basin have missing rows. For instance, there is no
  # year 400 in the Time column for the Deccan NGA. This looks like an error because if there would
  # be missing data for this century, the row should be there containing all "NA" rather than missing.
  # Again, for the sake of reproducing the original analysis, we will not correct the mistake
  # here but do it in our own analyses.

#bar chart paired
my.values<-mean(1000*out[,6],na.rm=TRUE)
err1<-1.96*std.error(1000*out[,6],na.rm=TRUE)
x <- barplot(my.values,names.arg="After - before moralizing gods",ylim = c(-1.3,1),
ylab="Rate of increase in social complexity (SC/ky)")
arrows(x,my.values-err1 ,x,my.values+err1, code=3, angle=90, length=.1)
abline(h=0)


##?OUR_COMMENT:: Here are the main results, that is +/- 2000 years Pre-/Post-MG.
# Note the DoF - 199, coming from only 12 NGAs with multiple timepoints (the time points are
# dependent).
print("Original t-test for MG at t = 0 with +/-2000 time-span")
print(t.test(out[,3], out[,2],paired=TRUE))

print("SC rate of change Pre/Post-MG")
print(mean(out[,2], na.rm = T)/mean(out[,3],na.rm = T))
print("########################################################################")
##?OUR_COMMENT:: How many rows of Pre-/Post-MG data?
NROW(out[(!is.na(out[,2]+out[,3])),])

##?OUR_COMMENT:: How many Pre-/Post-MG data per NGA?
out.na <- out[!is.na(out$Difference),]
out.na %>% group_by(NGA) %>% summarise(length(Difference))


##?OUR_COMMENT:: Was the rate of missigness equal across NGAs?
  # First select NGAs used in the t-test
NGAs <- c("Deccan", "Kachi Plain", "Kansai", "Konya Plain", "Latium", 
          "Middle Yellow River Valley", "Niger Inland Delta", "Orkhon Valley",
          "Paris Basin", "Sogdiana", "Susiana", "Upper Egypt")

miss <- matrix(NA,12,1)
for(i in 1:length(NGAs)){
  ot <- out[out$NGA == NGAs[i],]
  miss[i] <- NROW(ot[(is.na(ot[,2]+ot[,3])),])
  }
cbind(NGAs,miss) ##?OUR_COMMENT:: missing centuries per NGA


##?OUR_COMMENT:: Plot a histogram of the rate of SC change.

d <- d0
d <- as.data.frame(d)

ggplot() + 
  geom_histogram(data = d, aes(x = d), fill = "black", colour = "black", binwidth = 0.1,
                 alpha = 0.5) + 
  xlim(c(-2,2)) +
  ylim(c(0,100)) +
  geom_vline(xintercept = 0, color = "black", size=1, alpha = 1) +   
  labs(x="Rate of SC change", y="Frequency") +
  theme_bw() + 
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = rel(1.5)),        
    axis.line = element_line(colour = "black"),
    legend.position = c(0.95,0.95),
    legend.justification = c("right", "top"),
    legend.key.size = unit(0.8, "cm"),
    legend.text = element_text(size = rel(1)),
    axis.title = element_text(size = rel(1.1)),
    axis.text.y= element_text(size = rel(1.1)),
    axis.text.x= element_text(size = rel(1.1)),
    plot.margin=unit(c(1,1,1,1),"cm")) 


##?OUR_COMMENT:: If needed, save the plot
ggsave("./temp/Histogram_MG0.png",width = 4,height = 2.7, dpi = 300)

#______________________________________________________________________________________________
#### 1.4.2 700 years ####
##?OUR_COMMENT:: The same analysis for +/- 700 years

out<-read.table("./temp/FullRates.csv", sep=",", header=TRUE, stringsAsFactors = TRUE)

out$Difference<-out[,3]-out[,2]

for(i in 2:length(out[,5])){
  out[i,7]<-out[i,5]-out[i-1,5]
}
out <-subset(out, out[,7]!=0) #getting rid of bug when the final row repeats in each NGA

write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE)

out <-subset(out, out[,5]<750) #Change this to modify time-window restriction from 700 years
# pre/post moralizing gods (<750) or # out to use full
# time-window

write.csv(out, file="./temp/EqualRates.csv",  row.names=FALSE)

#bar chart paired
#my.values<-mean(1000*out[,6],na.rm=TRUE)
#err1<-1.96*std.error(1000*out[,6],na.rm=TRUE)
#x <- barplot(my.values,names.arg="After - before moralizing gods",ylim = c(-1.3,1),
#ylab="Rate of increase in social complexity (SC/ky)")
#arrows(x,my.values-err1 ,x,my.values+err1, code=3, angle=90, length=.1)
#abline(h=0)
print("Original t-test for MG at t = 0 with +/-700 time-span")
print(t.test(out[,3], out[,2],paired=TRUE))

print("SC rate of change Pre/Post-MG")
print(mean(out[,2], na.rm = T)/mean(out[,3],na.rm = T))
print("########################################################################")
##?OUR_COMMENT:: How many Pre-/Post-MG data per NGA? Note the problem with Deccan and Paris Basin
out.na <- out[!is.na(out$Difference),]
out.na %>% group_by(NGA) %>% summarise(length(Difference))


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
### 1.5. Dating uncertainty ####

##?OUR_COMMENT:: Whitehouse et al. provided robustness check against dating uncertainty, randomly sampling time
  # within the polity in which MG first appeared. How many polities could potentially shift MGs back in time?

matched <- matrix(NA, nrow=12, ncol=1)
unmatched <- matrix(NA, nrow=12, ncol=1)
for(i in 1:length(NGAs)){
  dt <- dat[dat$NGA == NGAs[i],]
  firstMG.row <- match(1, dt$MoralisingGods)  ##?OUR_COMMENT:: find first MG appearance
  pol1 <- dt$PolID[firstMG.row]  ##?OUR_COMMENT:: find polity with first MG
  if(dt$PolID[firstMG.row-1]==pol1){ ##?OUR_COMMENT:: MG appeared during polity
    matched[i] <- as.character(pol1)}
  else if(dt$PolID[firstMG.row-1]!=pol1){ ##?OUR_COMMENT:: MG appeared at beginning of polity
    unmatched[i] <- as.character(pol1)}
}

##?OUR_COMMENT:: For 11 NGAs, MGs appeared during the polity's first century
cbind(matched,unmatched)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## 2. Excluding NGAs with conquest ####

##?OUR_COMMENT:: From now on, we present our own re-analysis of the original script.
# We will keep up with the t-test, but also fix the errors mentioned in Section 1.


##?OUR_COMMENT:: In this section, we will check the robustness of Whitehouse et al. results
# by excluding NGAs that acquired MGs trhough conquest by larger empires (Deccan, Kachi Plain,
# Sogdiana) 


##?OUR_COMMENT:: First, correct for missing centuries at Deccan and Paris Basin

dat.cor <- dat
NGAs <- c("Deccan", "Kachi Plain", "Kansai", "Konya Plain", "Latium",
          "Middle Yellow River Valley", "Niger Inland Delta", "Orkhon Valley",
          "Paris Basin", "Sogdiana", "Susiana", "Upper Egypt")
dat.cor <- dat.cor[dat.cor$NGA %in% NGAs,]

dat.cor$Time[dat.cor$NGA == "Deccan" & (dat.cor$Time > 300)] <- 
  dat.cor$Time[dat.cor$NGA == "Deccan" & (dat.cor$Time > 300)] - 100
dat.cor$Time[dat.cor$NGA == "Deccan" & (dat.cor$Time > 900)] <- 
  dat.cor$Time[dat.cor$NGA == "Deccan" & (dat.cor$Time > 900)] - 100

dat.cor$Time[dat.cor$NGA == "Paris Basin" & (dat.cor$Time > 300)] <- 
  dat.cor$Time[dat.cor$NGA == "Paris Basin" & (dat.cor$Time > 300)] - 100

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
### 2.1 +-700 years time span ####


NGAs <- c("Kansai", "Konya Plain", "Latium", "Paris Basin", 
                 "Middle Yellow River Valley","Niger Inland Delta", "Orkhon Valley",
                 "Susiana", "Upper Egypt")

out <- matrix(NA, nrow=0, ncol=5)
for(i in 1:length(NGAs)){
  dt <- dat.cor[dat.cor$NGA == NGAs[i],]
  MG<-subset(dt,MoralisingGods=="1") #Replace "MoralisingGods" with "DoctrinalMode" 
  #or "Writing" to do these analyses
  #library(dplyr)
  #MG<-as.data.frame(MG %>% group_by(PolID) %>% sample_n(size = 1)) 
  #randomly samples so there is only one century per polity
  MGAppear<-subset(MG, Time==min(Time))
  for(j in 1: 100){
    Earliest<-subset(dt, Time==MGAppear$Time-j*100) ##?OUR_COMMENT:: Potentially confusing label,
    # it's not earliest, this is just the number of centuries before MG
    Latest<-subset(dt, Time==MGAppear$Time+j*100) ##?OUR_COMMENT:: Same as above
    rates<-cbind(as.character(MGAppear$NGA),
                 ifelse(class(Earliest$Time)=="NULL","NA",
                        (MGAppear$Mean-Earliest$Mean)/(MGAppear$Time-Earliest$Time)),
                 ifelse(class(Latest$Time)=="NULL","NA",
                        ((Latest$Mean-MGAppear$Mean)/(Latest$Time-MGAppear$Time))),
                 (MGAppear$End-MGAppear$Start),j*100)
    out <- rbind(out,rates)
  }
  out <- rbind(out,rates)
}
colnames(out)<-c("NGA","PreRate","PostRate","MGUncertainty","TimeWindow")

#mean(out$MGUncertainty) #Use this while replacing "MoralisingGods" as above to get uncertainty
# values for time-series

write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE) #Exporting/importing to force it to
# read as numeric (there is probably a more elegant way to do this)

out<-read.table("./temp/FullRates.csv", sep=",", header=TRUE, stringsAsFactors = TRUE)

out$Difference<-out[,3]-out[,2]

for(i in 1:length(out[,5])){
  ##?OUR_COMMENT:: First, fix the deleted first row discussed on lines 163:170
  if(i == 1){
    out[i,7]<- -9900
  }
  else if(i >1){
    out[i,7]<-out[i,5]-out[i-1,5]}}
out <-subset(out, out[,7]!=0) #getting rid of bug when the final row repeats in each NGA

write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE)

out <-subset(out, out[,5]<750) #Change this to modify time-window restriction from 700 years
# pre/post moralizing gods (<750) or # out to use full time-window

write.csv(out, file="./temp/EqualRates_uncq.csv",  row.names=FALSE)

print("t-test for unconquered NGAs with +/-700 time-span")
print(t.test(out[,3], out[,2],paired=TRUE))

print("SC rate of change Pre/Post-MG")
print(mean(out[,2], na.rm = T)/mean(out[,3],na.rm = T))
print("########################################################################")
##?OUR_COMMENT:: How many rows of Pre-/Post-MG data?
NROW(out[(!is.na(out[,2]+out[,3])),])

##?OUR_COMMENT:: How many Pre-/Post-MG data per NGA?
out.na <- out[!is.na(out$Difference),]
out.na %>% group_by(NGA) %>% summarise(length(Difference))


### 2.2 +-2000 years time span ####

NGAs <- c("Kansai", "Konya Plain", "Latium", "Paris Basin", 
          "Middle Yellow River Valley","Niger Inland Delta", "Orkhon Valley",
          "Susiana", "Upper Egypt")

out <- matrix(NA, nrow=0, ncol=5)
for(i in 1:length(NGAs)){
  dt <- dat.cor[dat.cor$NGA == NGAs[i],]
  MG<-subset(dt,MoralisingGods=="1") #Replace "MoralisingGods" with "DoctrinalMode" 
  #or "Writing" to do these analyses
  #library(dplyr)
  #MG<-as.data.frame(MG %>% group_by(PolID) %>% sample_n(size = 1)) 
  #randomly samples so there is only one century per polity
  MGAppear<-subset(MG, Time==min(Time))
  for(j in 1: 100){
    Earliest<-subset(dt, Time==MGAppear$Time-j*100) ##?OUR_COMMENT:: Potentially confusing label,
    # it's not earliest, this is just the number of centuries before MG
    Latest<-subset(dt, Time==MGAppear$Time+j*100) ##?OUR_COMMENT:: Same as above
    rates<-cbind(as.character(MGAppear$NGA),
                 ifelse(class(Earliest$Time)=="NULL","NA",
                        (MGAppear$Mean-Earliest$Mean)/(MGAppear$Time-Earliest$Time)),
                 ifelse(class(Latest$Time)=="NULL","NA",
                        ((Latest$Mean-MGAppear$Mean)/(Latest$Time-MGAppear$Time))),
                 (MGAppear$End-MGAppear$Start),j*100)
    out <- rbind(out,rates)
  }
  out <- rbind(out,rates)
}
colnames(out)<-c("NGA","PreRate","PostRate","MGUncertainty","TimeWindow")

#mean(out$MGUncertainty) #Use this while replacing "MoralisingGods" as above to get uncertainty
# values for time-series

write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE) #Exporting/importing to force it to
# read as numeric (there is probably a more elegant way to do this)

out<-read.table("./temp/FullRates.csv", sep=",", header=TRUE, stringsAsFactors = TRUE)

out$Difference<-out[,3]-out[,2]

for(i in 1:length(out[,5])){
  ##?OUR_COMMENT:: First, fix the deleted first row discussed on lines 163:170
  if(i == 1){
    out[i,7]<- -9900
  }
  else if(i >1){
    out[i,7]<-out[i,5]-out[i-1,5]}}
out <-subset(out, out[,7]!=0) #getting rid of bug when the final row repeats in each NGA

write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE)

out <-subset(out, out[,5]<2050) #Change this to modify time-window restriction from 700 years
# pre/post moralizing gods (<750) or # out to use full time-window

write.csv(out, file="./temp/EqualRates_uncq.csv",  row.names=FALSE)

print("t-test for unconquered NGAs with +/-2000 time-span")
print(t.test(out[,3], out[,2],paired=TRUE))

print("SC rate of change Pre/Post-MG")
print(mean(out[,2], na.rm = T)/mean(out[,3],na.rm = T))
print("########################################################################")
##?OUR_COMMENT:: How many rows of Pre-/Post-MG data?
NROW(out[(!is.na(out[,2]+out[,3])),])

##?OUR_COMMENT:: How many Pre-/Post-MG data per NGA?
out.na <- out[!is.na(out$Difference),]
out.na %>% group_by(NGA) %>% summarise(length(Difference))


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## 3. Excluding NGAs with mission and conquest ####

##?OUR_COMMENT:: First, correct for missing centuries at Deccan and Paris Basin

dat.cor <- dat
NGAs <- c("Deccan", "Kachi Plain", "Kansai", "Konya Plain", "Latium",
          "Middle Yellow River Valley", "Niger Inland Delta", "Orkhon Valley",
          "Paris Basin", "Sogdiana", "Susiana", "Upper Egypt")
dat.cor <- dat.cor[dat.cor$NGA %in% NGAs,]

dat.cor$Time[dat.cor$NGA == "Deccan" & (dat.cor$Time > 300)] <- 
  dat.cor$Time[dat.cor$NGA == "Deccan" & (dat.cor$Time > 300)] - 100
dat.cor$Time[dat.cor$NGA == "Deccan" & (dat.cor$Time > 900)] <- 
  dat.cor$Time[dat.cor$NGA == "Deccan" & (dat.cor$Time > 900)] - 100

dat.cor$Time[dat.cor$NGA == "Paris Basin" & (dat.cor$Time > 300)] <- 
  dat.cor$Time[dat.cor$NGA == "Paris Basin" & (dat.cor$Time > 300)] - 100


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
### 3.1 +- 2000 years ####

NGAs <- c("Kansai", "Konya Plain", "Latium", "Paris Basin", 
                 "Middle Yellow River Valley",
                 "Susiana", "Upper Egypt")

out <- matrix(NA, nrow=0, ncol=5)
for(i in 1:length(NGAs)){
  dt <- dat.cor[dat.cor$NGA == NGAs[i],]
  MG<-subset(dt,MoralisingGods=="1") #Replace "MoralisingGods" with "DoctrinalMode" 
  #or "Writing" to do these analyses
  #library(dplyr)
  #MG<-as.data.frame(MG %>% group_by(PolID) %>% sample_n(size = 1)) 
  #randomly samples so there is only one century per polity
  MGAppear<-subset(MG, Time==min(Time))
  for(j in 1: 100){
    Earliest<-subset(dt, Time==MGAppear$Time-j*100) ##?OUR_COMMENT:: Potentially confusing label,
    # it's not earliest, this is just the number of centuries before MG
    Latest<-subset(dt, Time==MGAppear$Time+j*100) ##?OUR_COMMENT:: Same as above
    rates<-cbind(as.character(MGAppear$NGA),
                 ifelse(class(Earliest$Time)=="NULL","NA",
                        (MGAppear$Mean-Earliest$Mean)/(MGAppear$Time-Earliest$Time)),
                 ifelse(class(Latest$Time)=="NULL","NA",
                        ((Latest$Mean-MGAppear$Mean)/(Latest$Time-MGAppear$Time))),
                 (MGAppear$End-MGAppear$Start),j*100)
    out <- rbind(out,rates)
  }
  out <- rbind(out,rates)
}
colnames(out)<-c("NGA","PreRate","PostRate","MGUncertainty","TimeWindow")

#mean(out$MGUncertainty) #Use this while replacing "MoralisingGods" as above to get uncertainty
# values for time-series

write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE) #Exporting/importing to force it to
# read as numeric (there is probably a more elegant way to do this)

out<-read.table("./temp/FullRates.csv", sep=",", header=TRUE, stringsAsFactors = TRUE)

out$Difference<-out[,3]-out[,2]

for(i in 1:length(out[,5])){
  ##?OUR_COMMENT:: First, fix the deleted first row discussed on lines 163:170
  if(i == 1){
    out[i,7]<- -9900
  }
  else if(i >1){
    out[i,7]<-out[i,5]-out[i-1,5]}}
out <-subset(out, out[,7]!=0) #getting rid of bug when the final row repeats in each NGA

write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE)

out <-subset(out, out[,5]<750) #Change this to modify time-window restriction from 700 years
# pre/post moralizing gods (<750) or # out to use full time-window

write.csv(out, file="./temp/EqualRates_uncq.csv",  row.names=FALSE)

print("t-test for 6 NGAs where MG emerge at t = 0 for time span +-700 years")
print(t.test(out[,3], out[,2],paired=TRUE))

print("SC rate of change Pre/Post-MG")
print(mean(out[,2], na.rm = T)/mean(out[,3],na.rm = T))
print("########################################################################")
##?OUR_COMMENT:: How many rows of Pre-/Post-MG data?
NROW(out[(!is.na(out[,2]+out[,3])),])

##?OUR_COMMENT:: How many Pre-/Post-MG data per NGA?
out.na <- out[!is.na(out$Difference),]
out.na %>% group_by(NGA) %>% summarise(length(Difference))


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
### 3.2 +- 700 years ####

NGAs <- c("Kansai", "Konya Plain", "Latium", "Paris Basin", 
          "Middle Yellow River Valley",
          "Susiana", "Upper Egypt")

out <- matrix(NA, nrow=0, ncol=5)
for(i in 1:length(NGAs)){
  dt <- dat.cor[dat.cor$NGA == NGAs[i],]
  MG<-subset(dt,MoralisingGods=="1") #Replace "MoralisingGods" with "DoctrinalMode" 
  #or "Writing" to do these analyses
  #library(dplyr)
  #MG<-as.data.frame(MG %>% group_by(PolID) %>% sample_n(size = 1)) 
  #randomly samples so there is only one century per polity
  MGAppear<-subset(MG, Time==min(Time))
  for(j in 1: 100){
    Earliest<-subset(dt, Time==MGAppear$Time-j*100) ##?OUR_COMMENT:: Potentially confusing label,
    # it's not earliest, this is just the number of centuries before MG
    Latest<-subset(dt, Time==MGAppear$Time+j*100) ##?OUR_COMMENT:: Same as above
    rates<-cbind(as.character(MGAppear$NGA),
                 ifelse(class(Earliest$Time)=="NULL","NA",
                        (MGAppear$Mean-Earliest$Mean)/(MGAppear$Time-Earliest$Time)),
                 ifelse(class(Latest$Time)=="NULL","NA",
                        ((Latest$Mean-MGAppear$Mean)/(Latest$Time-MGAppear$Time))),
                 (MGAppear$End-MGAppear$Start),j*100)
    out <- rbind(out,rates)
  }
  out <- rbind(out,rates)
}
colnames(out)<-c("NGA","PreRate","PostRate","MGUncertainty","TimeWindow")

#mean(out$MGUncertainty) #Use this while replacing "MoralisingGods" as above to get uncertainty
# values for time-series

write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE) #Exporting/importing to force it to
# read as numeric (there is probably a more elegant way to do this)

out<-read.table("./temp/FullRates.csv", sep=",", header=TRUE, stringsAsFactors = TRUE)

out$Difference<-out[,3]-out[,2]

for(i in 1:length(out[,5])){
  ##?OUR_COMMENT:: First, fix the deleted first row discussed on lines 163:170
  if(i == 1){
    out[i,7]<- -9900
  }
  else if(i >1){
    out[i,7]<-out[i,5]-out[i-1,5]}}
out <-subset(out, out[,7]!=0) #getting rid of bug when the final row repeats in each NGA

write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE)

out <-subset(out, out[,5]<2050) #Change this to modify time-window restriction from 700 years
# pre/post moralizing gods (<750) or # out to use full time-window

write.csv(out, file="./temp/EqualRates_uncq.csv",  row.names=FALSE)

print("t-test for 6 NGAs where MG emerge at t = 0 for time span +-2000 years")
print(t.test(out[,3], out[,2],paired=TRUE))

print("SC rate of change Pre/Post-MG")
print(mean(out[,2], na.rm = T)/mean(out[,3],na.rm = T))
print("########################################################################")
##?OUR_COMMENT:: How many rows of Pre-/Post-MG data?
NROW(out[(!is.na(out[,2]+out[,3])),])

##?OUR_COMMENT:: How many Pre-/Post-MG data per NGA?
out.na <- out[!is.na(out$Difference),]
out.na %>% group_by(NGA) %>% summarise(length(Difference))




#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## 4. Forward-bias correction for +- 2000 years ####

##?OUR_COMMENT:: In this section, we will try to correct for forward bias,
# i.e., the fact that if we take the first known occurrence of MGs, it is very unlikely
# that it is actually the oldest occurrence. To do so, we need to exclude NGAs 
# that acquired MGs trhough conquest by larger empires (Deccan, Kachi Plain,
# Sogdiana) or through mission (Kansai and Niger Inland Delta).

##?OUR_COMMENT:: First, correct for missing centuries at Deccan and Paris Basin

dat.cor <- dat
NGAs <- c("Deccan", "Kachi Plain", "Kansai", "Konya Plain", "Latium",
          "Middle Yellow River Valley", "Niger Inland Delta", "Orkhon Valley",
          "Paris Basin", "Sogdiana", "Susiana", "Upper Egypt")
dat.cor <- dat.cor[dat.cor$NGA %in% NGAs,]

dat.cor$Time[dat.cor$NGA == "Deccan" & (dat.cor$Time > 300)] <- 
  dat.cor$Time[dat.cor$NGA == "Deccan" & (dat.cor$Time > 300)] - 100
dat.cor$Time[dat.cor$NGA == "Deccan" & (dat.cor$Time > 900)] <- 
  dat.cor$Time[dat.cor$NGA == "Deccan" & (dat.cor$Time > 900)] - 100

dat.cor$Time[dat.cor$NGA == "Paris Basin" & (dat.cor$Time > 300)] <- 
  dat.cor$Time[dat.cor$NGA == "Paris Basin" & (dat.cor$Time > 300)] - 100

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
### 4.1 Full sample - MG t-100 ####

NGAs <- c("Deccan", "Kachi Plain", "Kansai", "Konya Plain", "Latium",
          "Middle Yellow River Valley", "Niger Inland Delta", "Orkhon Valley",
          "Paris Basin", "Sogdiana", "Susiana", "Upper Egypt")

FB <- 100 ##?OUR_COMMENT:: this is Forward Bias in years - we will subtract this number from
#MG appearence times. # Can be changed to e.g. 300

##?OUR_COMMENT:: Adapting the original script...
out <- matrix(NA, nrow=0, ncol=5)
for(i in 1:length(NGAs)){
  dt <- dat.cor[dat.cor$NGA == NGAs[i],]
  MG<-subset(dt,MoralisingGods=="1")
  ##?OUR_COMMENT:: This is our addition - we have to create a new subset MGAppear.FB,
  # which will contain Mean SC for the time of MG's first appearance minus forward bias(FB)
  MGAppear<-subset(MG, Time==min(Time))
  t <- min(MG$Time)-FB
  MGAppear.FB<-subset(dt, Time==t)
  ##?OUR_COMMENT:: Now calculate the rate of SC change for the time-shifted data
  for(j in 1: 100){
    Earliest<-subset(dt, Time==MGAppear.FB$Time-j*100) 
    Latest<-subset(dt, Time==MGAppear.FB$Time+j*100)
    rates<-cbind(as.character(MGAppear$NGA), ##?OUR_COMMENT:: had to add 'as.character()' here
                 ifelse(class(Earliest$Time)=="NULL","NA",
                        (MGAppear.FB$Mean-Earliest$Mean)/(MGAppear.FB$Time-Earliest$Time)),
                 ifelse(class(Latest$Time)=="NULL","NA",
                        ((Latest$Mean-MGAppear.FB$Mean)/(Latest$Time-MGAppear.FB$Time))),
                 (MGAppear.FB$End-MGAppear.FB$Start),
                 j*100)
    out <- rbind(out,rates)
  }
  out <- rbind(out,rates) ##?OUR_COMMENT:: This line is responsible for the "bug" described
  # a few lines below. It could be removed, but let's stick with
  # the original script.
}
colnames(out)<-c("NGA","PreRate","PostRate","MGUncertainty","TimeWindow")

write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE) #Exporting/importing to force it to read
#as numeric (there is probably a more
#elegant way to do this)

out<-read.table("./temp/FullRates.csv", sep=",", header=TRUE, stringsAsFactors = TRUE)

out$Difference<-out[,3]-out[,2]

for(i in 1:length(out[,5])){
  ##?OUR_COMMENT:: First, fix the deleted first row discussed on lines 163:170
  if(i == 1){
    out[i,7]<- -9900
  }
  else if(i >1){
    out[i,7]<-out[i,5]-out[i-1,5]}}
out <-subset(out, out[,7]!=0) #getting rid of bug when the final row repeats in each NGA

write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE)

out<-read.table("./temp/FullRates.csv", sep=",", header=TRUE, stringsAsFactors = TRUE)

out <-subset(out, out[,5]<2050) #Change this to modify time-window restriction from 700 years
#pre/post moralizing gods (<750) or 2050 out to use full time-window
write.csv(out, file="./temp/EqualRates_FB_100.csv",  row.names=FALSE)


##?OUR_COMMENT:: Result for Pre/Post-MG with MG appearance shifted
#  100 years back in time.
print("t-test for 12 NGAs shifting MG 't - 100' for time span +-2000 years")
print(t.test(out[,3], out[,2],paired=TRUE))

print("SC rate of change Pre/Post-MG")
print(mean(out[,2], na.rm = T)/mean(out[,3],na.rm = T))
print("########################################################################")
##?OUR_COMMENT:: How many rows of Pre-/Post-MG data?
NROW(out[(!is.na(out[,2]+out[,3])),])

##?OUR_COMMENT:: How many Pre-/Post-MG data per NGA?
out.na <- out[!is.na(out$Difference),]
out.na %>% group_by(NGA) %>% summarise(length(Difference))

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
### 4.2 Full sample - MG t-300 ####

NGAs <- c("Deccan", "Kachi Plain", "Kansai", "Konya Plain", "Latium",
          "Middle Yellow River Valley", "Niger Inland Delta", "Orkhon Valley",
          "Paris Basin", "Sogdiana", "Susiana", "Upper Egypt")

FB <- 300 ##?OUR_COMMENT:: this is Forward Bias in years - we will subtract this number from
#MG appearence times.

##?OUR_COMMENT:: Adapting the original script...
out <- matrix(NA, nrow=0, ncol=5)
for(i in 1:length(NGAs)){
  dt <- dat.cor[dat.cor$NGA == NGAs[i],]
  MG<-subset(dt,MoralisingGods=="1")
  ##?OUR_COMMENT:: This is our addition - we have to create a new subset MGAppear.FB,
  # which will contain Mean SC for the time of MG's first appearance minus forward bias(FB)
  MGAppear<-subset(MG, Time==min(Time))
  t <- min(MG$Time)-FB
  MGAppear.FB<-subset(dt, Time==t)
  ##?OUR_COMMENT:: Now calculate the rate of SC change for the time-shifted data
  for(j in 1: 100){
    Earliest<-subset(dt, Time==MGAppear.FB$Time-j*100) 
    Latest<-subset(dt, Time==MGAppear.FB$Time+j*100)
    rates<-cbind(as.character(MGAppear$NGA), ##?OUR_COMMENT:: had to add 'as.character()' here
                 ifelse(class(Earliest$Time)=="NULL","NA",
                        (MGAppear.FB$Mean-Earliest$Mean)/(MGAppear.FB$Time-Earliest$Time)),
                 ifelse(class(Latest$Time)=="NULL","NA",
                        ((Latest$Mean-MGAppear.FB$Mean)/(Latest$Time-MGAppear.FB$Time))),
                 (MGAppear.FB$End-MGAppear.FB$Start),
                 j*100)
    out <- rbind(out,rates)
  }
  out <- rbind(out,rates) ##?OUR_COMMENT:: This line is responsible for the "bug" described
  # a few lines below. It could be removed, but let's stick with
  # the original script.
}
colnames(out)<-c("NGA","PreRate","PostRate","MGUncertainty","TimeWindow")

write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE) #Exporting/importing to force it to read
#as numeric (there is probably a more
#elegant way to do this)

out<-read.table("./temp/FullRates.csv", sep=",", header=TRUE, stringsAsFactors = TRUE)

out$Difference<-out[,3]-out[,2]


for(i in 1:length(out[,5])){
  ##?OUR_COMMENT:: First, fix the deleted first row discussed on lines 163:170
  if(i == 1){
    out[i,7]<- -9900
  }
  else if(i >1){
    out[i,7]<-out[i,5]-out[i-1,5]}}
out <-subset(out, out[,7]!=0) #getting rid of bug when the final row repeats in each NGA

write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE)

out<-read.table("./temp/FullRates.csv", sep=",", header=TRUE, stringsAsFactors = TRUE)

out <-subset(out, out[,5]<2050) #Change this to modify time-window restriction from 700 years
#pre/post moralizing gods (<750) or 2050 out to use full time-window
write.csv(out, file="./temp/EqualRates_FB_300.csv",  row.names=FALSE)


##?OUR_COMMENT:: Result for Pre/Post-MG with MG appearance shifted
#  300 years back in time.
print("t-test for 12 NGAs shifting MG 't - 300' for time span +-2000 years")
print(t.test(out[,3], out[,2],paired=TRUE))

print("SC rate of change Pre/Post-MG")
print(mean(out[,2], na.rm = T)/mean(out[,3],na.rm = T))
print("########################################################################")
##?OUR_COMMENT:: How many Pre-/Post-MG data per NGA?
out.na <- out[!is.na(out$Difference),]
out.na %>% group_by(NGA) %>% summarise(length(Difference))



#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
### 4.3 Half sample - MG t-100 ####

NGAs <- c("Konya Plain", "Latium", "Paris Basin", 
          "Middle Yellow River Valley", 
          "Susiana", "Upper Egypt")

FB <- 100 ##?OUR_COMMENT:: this is Forward Bias in years - we will subtract this number from
#MG appearence times. # Can be changed to e.g. 300

##?OUR_COMMENT:: Adapting the original script...
out <- matrix(NA, nrow=0, ncol=5)
for(i in 1:length(NGAs)){
  dt <- dat.cor[dat.cor$NGA == NGAs[i],]
  MG<-subset(dt,MoralisingGods=="1")
  ##?OUR_COMMENT:: This is our addition - we have to create a new subset MGAppear.FB,
  # which will contain Mean SC for the time of MG's first appearance minus forward bias(FB)
  MGAppear<-subset(MG, Time==min(Time))
  t <- min(MG$Time)-FB
  MGAppear.FB<-subset(dt, Time==t)
  ##?OUR_COMMENT:: Now calculate the rate of SC change for the time-shifted data
  for(j in 1: 100){
    Earliest<-subset(dt, Time==MGAppear.FB$Time-j*100) 
    Latest<-subset(dt, Time==MGAppear.FB$Time+j*100)
    rates<-cbind(as.character(MGAppear$NGA), ##?OUR_COMMENT:: had to add 'as.character()' here
                 ifelse(class(Earliest$Time)=="NULL","NA",
                        (MGAppear.FB$Mean-Earliest$Mean)/(MGAppear.FB$Time-Earliest$Time)),
                 ifelse(class(Latest$Time)=="NULL","NA",
                        ((Latest$Mean-MGAppear.FB$Mean)/(Latest$Time-MGAppear.FB$Time))),
                 (MGAppear.FB$End-MGAppear.FB$Start),
                 j*100)
    out <- rbind(out,rates)
  }
  out <- rbind(out,rates) ##?OUR_COMMENT:: This line is responsible for the "bug" described
  # a few lines below. It could be removed, but let's stick with
  # the original script.
}
colnames(out)<-c("NGA","PreRate","PostRate","MGUncertainty","TimeWindow")

write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE) #Exporting/importing to force it to read
#as numeric (there is probably a more
#elegant way to do this)

out<-read.table("./temp/FullRates.csv", sep=",", header=TRUE, stringsAsFactors = TRUE)

out$Difference<-out[,3]-out[,2]

for(i in 1:length(out[,5])){
  ##?OUR_COMMENT:: First, fix the deleted first row discussed on lines 163:170
  if(i == 1){
    out[i,7]<- -9900
  }
  else if(i >1){
    out[i,7]<-out[i,5]-out[i-1,5]}}
out <-subset(out, out[,7]!=0) #getting rid of bug when the final row repeats in each NGA

write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE)

out<-read.table("./temp/FullRates.csv", sep=",", header=TRUE, stringsAsFactors = TRUE)

out <-subset(out, out[,5]<2050) #Change this to modify time-window restriction from 700 years
#pre/post moralizing gods (<750) or 2050 out to use full time-window
write.csv(out, file="./temp/EqualRates_uncq_FB_100.csv",  row.names=FALSE)


##?OUR_COMMENT:: Result for Pre/Post-MG with MG appearance shifted
#  100 years back in time.
print("t-test for 6 NGAs where MG emerge at t - 100 for time span +-2000 years")
print(t.test(out[,3], out[,2],paired=TRUE))


print("SC rate of change Pre/Post-MG")
print(mean(out[,2], na.rm = T)/mean(out[,3],na.rm = T))
print("########################################################################")
##?OUR_COMMENT:: How many rows of Pre-/Post-MG data?
NROW(out[(!is.na(out[,2]+out[,3])),])

##?OUR_COMMENT:: How many Pre-/Post-MG data per NGA?
out.na <- out[!is.na(out$Difference),]
out.na %>% group_by(NGA) %>% summarise(length(Difference))

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
### 4.4 Half sample - MG t-300 ####
NGAs <- c("Konya Plain", "Latium", "Paris Basin", 
          "Middle Yellow River Valley",  "Niger Inland Delta",
          "Susiana", "Upper Egypt")

FB <- 300 ##?OUR_COMMENT:: this is Forward Bias in years - we will subtract this number from
#MG appearence times.

##?OUR_COMMENT:: Adapting the original script...
out <- matrix(NA, nrow=0, ncol=5)
for(i in 1:length(NGAs)){
  dt <- dat.cor[dat.cor$NGA == NGAs[i],]
  MG<-subset(dt,MoralisingGods=="1")
  ##?OUR_COMMENT:: This is our addition - we have to create a new subset MGAppear.FB,
  # which will contain Mean SC for the time of MG's first appearance minus forward bias(FB)
  MGAppear<-subset(MG, Time==min(Time))
  t <- min(MG$Time)-FB
  MGAppear.FB<-subset(dt, Time==t)
  ##?OUR_COMMENT:: Now calculate the rate of SC change for the time-shifted data
  for(j in 1: 100){
    Earliest<-subset(dt, Time==MGAppear.FB$Time-j*100) 
    Latest<-subset(dt, Time==MGAppear.FB$Time+j*100)
    rates<-cbind(as.character(MGAppear$NGA), ##?OUR_COMMENT:: had to add 'as.character()' here
                 ifelse(class(Earliest$Time)=="NULL","NA",
                        (MGAppear.FB$Mean-Earliest$Mean)/(MGAppear.FB$Time-Earliest$Time)),
                 ifelse(class(Latest$Time)=="NULL","NA",
                        ((Latest$Mean-MGAppear.FB$Mean)/(Latest$Time-MGAppear.FB$Time))),
                 (MGAppear.FB$End-MGAppear.FB$Start),
                 j*100)
    out <- rbind(out,rates)
  }
  out <- rbind(out,rates) ##?OUR_COMMENT:: This line is responsible for the "bug" described
  # a few lines below. It could be removed, but let's stick with
  # the original script.
}
colnames(out)<-c("NGA","PreRate","PostRate","MGUncertainty","TimeWindow")

write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE) #Exporting/importing to force it to read
#as numeric (there is probably a more
#elegant way to do this)

out<-read.table("./temp/FullRates.csv", sep=",", header=TRUE, stringsAsFactors = TRUE)

out$Difference<-out[,3]-out[,2]


for(i in 1:length(out[,5])){
  ##?OUR_COMMENT:: First, fix the deleted first row discussed on lines 163:170
  if(i == 1){
    out[i,7]<- -9900
  }
  else if(i >1){
    out[i,7]<-out[i,5]-out[i-1,5]}}
out <-subset(out, out[,7]!=0) #getting rid of bug when the final row repeats in each NGA

write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE)

out<-read.table("./temp/FullRates.csv", sep=",", header=TRUE, stringsAsFactors = TRUE)

out <-subset(out, out[,5]<2050) #Change this to modify time-window restriction from 700 years
#pre/post moralizing gods (<750) or 2050 out to use full time-window
write.csv(out, file="./temp/EqualRates_uncq_FB_300.csv",  row.names=FALSE)


##?OUR_COMMENT:: Result for Pre/Post-MG with MG appearance shifted
#  300 years back in time.
print("t-test for 6 NGAs where MG emerge at t - 300  for time span +-2000 years")
print(t.test(out[,3], out[,2],paired=TRUE))

print("SC rate of change Pre/Post-MG")
print(mean(out[,2], na.rm = T)/mean(out[,3],na.rm = T))
print("########################################################################")
##?OUR_COMMENT:: How many Pre-/Post-MG data per NGA?
out.na <- out[!is.na(out$Difference),]
out.na %>% group_by(NGA) %>% summarise(length(Difference))




#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## 5. Forward-bias for +- 700 years ####

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
### 5.1 Full sample - MG t-100 ####

NGAs <- c("Deccan", "Kachi Plain", "Kansai", "Konya Plain", "Latium",
          "Middle Yellow River Valley", "Niger Inland Delta", "Orkhon Valley",
          "Paris Basin", "Sogdiana", "Susiana", "Upper Egypt")

FB <- 100 ##?OUR_COMMENT:: this is Forward Bias in years - we will subtract this number from
#MG appearence times. # Can be changed to e.g. 300

##?OUR_COMMENT:: Adapting the original script...
out <- matrix(NA, nrow=0, ncol=5)
for(i in 1:length(NGAs)){
  dt <- dat.cor[dat.cor$NGA == NGAs[i],]
  MG<-subset(dt,MoralisingGods=="1")
  ##?OUR_COMMENT:: This is our addition - we have to create a new subset MGAppear.FB,
  # which will contain Mean SC for the time of MG's first appearance minus forward bias(FB)
  MGAppear<-subset(MG, Time==min(Time))
  t <- min(MG$Time)-FB
  MGAppear.FB<-subset(dt, Time==t)
  ##?OUR_COMMENT:: Now calculate the rate of SC change for the time-shifted data
  for(j in 1: 100){
    Earliest<-subset(dt, Time==MGAppear.FB$Time-j*100) 
    Latest<-subset(dt, Time==MGAppear.FB$Time+j*100)
    rates<-cbind(as.character(MGAppear$NGA), ##?OUR_COMMENT:: had to add 'as.character()' here
                 ifelse(class(Earliest$Time)=="NULL","NA",
                        (MGAppear.FB$Mean-Earliest$Mean)/(MGAppear.FB$Time-Earliest$Time)),
                 ifelse(class(Latest$Time)=="NULL","NA",
                        ((Latest$Mean-MGAppear.FB$Mean)/(Latest$Time-MGAppear.FB$Time))),
                 (MGAppear.FB$End-MGAppear.FB$Start),
                 j*100)
    out <- rbind(out,rates)
  }
  out <- rbind(out,rates) ##?OUR_COMMENT:: This line is responsible for the "bug" described
  # a few lines below. It could be removed, but let's stick with
  # the original script.
}
colnames(out)<-c("NGA","PreRate","PostRate","MGUncertainty","TimeWindow")

write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE) #Exporting/importing to force it to read
#as numeric (there is probably a more
#elegant way to do this)

out<-read.table("./temp/FullRates.csv", sep=",", header=TRUE, stringsAsFactors = TRUE)

out$Difference<-out[,3]-out[,2]

for(i in 1:length(out[,5])){
  ##?OUR_COMMENT:: First, fix the deleted first row discussed on lines 163:170
  if(i == 1){
    out[i,7]<- -9900
  }
  else if(i >1){
    out[i,7]<-out[i,5]-out[i-1,5]}}
out <-subset(out, out[,7]!=0) #getting rid of bug when the final row repeats in each NGA

write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE)

out<-read.table("./temp/FullRates.csv", sep=",", header=TRUE, stringsAsFactors = TRUE)

out <-subset(out, out[,5]<750) #Change this to modify time-window restriction from 700 years
#pre/post moralizing gods (<750) or 2050 out to use full time-window
write.csv(out, file="./temp/EqualRates_uncq_FB_100.csv",  row.names=FALSE)


##?OUR_COMMENT:: Result for Pre/Post-MG with MG appearance shifted
#  100 years back in time.
print("t-test for 12 NGAs shifting MG 't - 100' for time span +-700 years")
print(t.test(out[,3], out[,2],paired=TRUE))


print("SC rate of change Pre/Post-MG")
print(mean(out[,2], na.rm = T)/mean(out[,3],na.rm = T))
print("########################################################################")
##?OUR_COMMENT:: How many rows of Pre-/Post-MG data?
NROW(out[(!is.na(out[,2]+out[,3])),])

##?OUR_COMMENT:: How many Pre-/Post-MG data per NGA?
out.na <- out[!is.na(out$Difference),]
out.na %>% group_by(NGA) %>% summarise(length(Difference))

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
### 5.2 Full sample - MG t-300 ####
NGAs <- c("Deccan", "Kachi Plain", "Kansai", "Konya Plain", "Latium",
          "Middle Yellow River Valley", "Niger Inland Delta", "Orkhon Valley",
          "Paris Basin", "Sogdiana", "Susiana", "Upper Egypt")

FB <- 300 ##?OUR_COMMENT:: this is Forward Bias in years - we will subtract this number from
#MG appearence times. # Can be changed to e.g. 300

##?OUR_COMMENT:: Adapting the original script...
out <- matrix(NA, nrow=0, ncol=5)
for(i in 1:length(NGAs)){
  dt <- dat.cor[dat.cor$NGA == NGAs[i],]
  MG<-subset(dt,MoralisingGods=="1")
  ##?OUR_COMMENT:: This is our addition - we have to create a new subset MGAppear.FB,
  # which will contain Mean SC for the time of MG's first appearance minus forward bias(FB)
  MGAppear<-subset(MG, Time==min(Time))
  t <- min(MG$Time)-FB
  MGAppear.FB<-subset(dt, Time==t)
  ##?OUR_COMMENT:: Now calculate the rate of SC change for the time-shifted data
  for(j in 1: 100){
    Earliest<-subset(dt, Time==MGAppear.FB$Time-j*100) 
    Latest<-subset(dt, Time==MGAppear.FB$Time+j*100)
    rates<-cbind(as.character(MGAppear$NGA), ##?OUR_COMMENT:: had to add 'as.character()' here
                 ifelse(class(Earliest$Time)=="NULL","NA",
                        (MGAppear.FB$Mean-Earliest$Mean)/(MGAppear.FB$Time-Earliest$Time)),
                 ifelse(class(Latest$Time)=="NULL","NA",
                        ((Latest$Mean-MGAppear.FB$Mean)/(Latest$Time-MGAppear.FB$Time))),
                 (MGAppear.FB$End-MGAppear.FB$Start),
                 j*100)
    out <- rbind(out,rates)
  }
  out <- rbind(out,rates) ##?OUR_COMMENT:: This line is responsible for the "bug" described
  # a few lines below. It could be removed, but let's stick with
  # the original script.
}
colnames(out)<-c("NGA","PreRate","PostRate","MGUncertainty","TimeWindow")

write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE) #Exporting/importing to force it to read
#as numeric (there is probably a more
#elegant way to do this)

out<-read.table("./temp/FullRates.csv", sep=",", header=TRUE, stringsAsFactors = TRUE)

out$Difference<-out[,3]-out[,2]


for(i in 1:length(out[,5])){
  ##?OUR_COMMENT:: First, fix the deleted first row discussed on lines 163:170
  if(i == 1){
    out[i,7]<- -9900
  }
  else if(i >1){
    out[i,7]<-out[i,5]-out[i-1,5]}}
out <-subset(out, out[,7]!=0) #getting rid of bug when the final row repeats in each NGA

write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE)

out<-read.table("./temp/FullRates.csv", sep=",", header=TRUE, stringsAsFactors = TRUE)

out <-subset(out, out[,5]<750) #Change this to modify time-window restriction from 700 years
#pre/post moralizing gods (<750) or 2050 out to use full time-window
write.csv(out, file="./temp/EqualRates_uncq_FB_300.csv",  row.names=FALSE)


##?OUR_COMMENT:: Result for Pre/Post-MG with MG appearance shifted
#  300 years back in time.
print("t-test for 12 NGAs shifting MG 't - 300' for time span +-700 years")
print(t.test(out[,3], out[,2],paired=TRUE))

print("SC rate of change Pre/Post-MG")
print(mean(out[,2], na.rm = T)/mean(out[,3],na.rm = T))
print("########################################################################")
##?OUR_COMMENT:: How many Pre-/Post-MG data per NGA?
out.na <- out[!is.na(out$Difference),]
out.na %>% group_by(NGA) %>% summarise(length(Difference))


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
### 5.3 Half sample - MG t-100 ####

NGAs <- c("Konya Plain", "Latium", "Paris Basin", 
          "Middle Yellow River Valley",
          "Susiana", "Upper Egypt")

FB <- 100 ##?OUR_COMMENT:: this is Forward Bias in years - we will subtract this number from
#MG appearence times. # Can be changed to e.g. 300

##?OUR_COMMENT:: Adapting the original script...
out <- matrix(NA, nrow=0, ncol=5)
for(i in 1:length(NGAs)){
  dt <- dat.cor[dat.cor$NGA == NGAs[i],]
  MG<-subset(dt,MoralisingGods=="1")
  ##?OUR_COMMENT:: This is our addition - we have to create a new subset MGAppear.FB,
  # which will contain Mean SC for the time of MG's first appearance minus forward bias(FB)
  MGAppear<-subset(MG, Time==min(Time))
  t <- min(MG$Time)-FB
  MGAppear.FB<-subset(dt, Time==t)
  ##?OUR_COMMENT:: Now calculate the rate of SC change for the time-shifted data
  for(j in 1: 100){
    Earliest<-subset(dt, Time==MGAppear.FB$Time-j*100) 
    Latest<-subset(dt, Time==MGAppear.FB$Time+j*100)
    rates<-cbind(as.character(MGAppear$NGA), ##?OUR_COMMENT:: had to add 'as.character()' here
                 ifelse(class(Earliest$Time)=="NULL","NA",
                        (MGAppear.FB$Mean-Earliest$Mean)/(MGAppear.FB$Time-Earliest$Time)),
                 ifelse(class(Latest$Time)=="NULL","NA",
                        ((Latest$Mean-MGAppear.FB$Mean)/(Latest$Time-MGAppear.FB$Time))),
                 (MGAppear.FB$End-MGAppear.FB$Start),
                 j*100)
    out <- rbind(out,rates)
  }
  out <- rbind(out,rates) ##?OUR_COMMENT:: This line is responsible for the "bug" described
  # a few lines below. It could be removed, but let's stick with
  # the original script.
}
colnames(out)<-c("NGA","PreRate","PostRate","MGUncertainty","TimeWindow")

write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE) #Exporting/importing to force it to read
#as numeric (there is probably a more
#elegant way to do this)

out<-read.table("./temp/FullRates.csv", sep=",", header=TRUE, stringsAsFactors = TRUE)

out$Difference<-out[,3]-out[,2]

for(i in 1:length(out[,5])){
  ##?OUR_COMMENT:: First, fix the deleted first row discussed on lines 163:170
  if(i == 1){
    out[i,7]<- -9900
  }
  else if(i >1){
    out[i,7]<-out[i,5]-out[i-1,5]}}
out <-subset(out, out[,7]!=0) #getting rid of bug when the final row repeats in each NGA

write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE)

out<-read.table("./temp/FullRates.csv", sep=",", header=TRUE, stringsAsFactors = TRUE)

out <-subset(out, out[,5]<750) #Change this to modify time-window restriction from 700 years
#pre/post moralizing gods (<750) or 2050 out to use full time-window
write.csv(out, file="./temp/EqualRates_uncq_FB_100.csv",  row.names=FALSE)


##?OUR_COMMENT:: Result for Pre/Post-MG with MG appearance shifted
#  100 years back in time.
print("t-test for 6 NGAs where MG emerge at t - 100 for time span +-700 years")
print(t.test(out[,3], out[,2],paired=TRUE))


print("SC rate of change Pre/Post-MG")
print(mean(out[,2], na.rm = T)/mean(out[,3],na.rm = T))
print("########################################################################")
##?OUR_COMMENT:: How many rows of Pre-/Post-MG data?
NROW(out[(!is.na(out[,2]+out[,3])),])

##?OUR_COMMENT:: How many Pre-/Post-MG data per NGA?
out.na <- out[!is.na(out$Difference),]
out.na %>% group_by(NGA) %>% summarise(length(Difference))

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
### 5.4 Half sample - MG t-300 ####
NGAs <- c("Konya Plain", "Latium", "Paris Basin", 
          "Middle Yellow River Valley",
          "Susiana", "Upper Egypt")

FB <- 300 ##?OUR_COMMENT:: this is Forward Bias in years - we will subtract this number from
#MG appearence times. # Can be changed to e.g. 300

##?OUR_COMMENT:: Adapting the original script...
out <- matrix(NA, nrow=0, ncol=5)
for(i in 1:length(NGAs)){
  dt <- dat.cor[dat.cor$NGA == NGAs[i],]
  MG<-subset(dt,MoralisingGods=="1")
  ##?OUR_COMMENT:: This is our addition - we have to create a new subset MGAppear.FB,
  # which will contain Mean SC for the time of MG's first appearance minus forward bias(FB)
  MGAppear<-subset(MG, Time==min(Time))
  t <- min(MG$Time)-FB
  MGAppear.FB<-subset(dt, Time==t)
  ##?OUR_COMMENT:: Now calculate the rate of SC change for the time-shifted data
  for(j in 1: 100){
    Earliest<-subset(dt, Time==MGAppear.FB$Time-j*100) 
    Latest<-subset(dt, Time==MGAppear.FB$Time+j*100)
    rates<-cbind(as.character(MGAppear$NGA), ##?OUR_COMMENT:: had to add 'as.character()' here
                 ifelse(class(Earliest$Time)=="NULL","NA",
                        (MGAppear.FB$Mean-Earliest$Mean)/(MGAppear.FB$Time-Earliest$Time)),
                 ifelse(class(Latest$Time)=="NULL","NA",
                        ((Latest$Mean-MGAppear.FB$Mean)/(Latest$Time-MGAppear.FB$Time))),
                 (MGAppear.FB$End-MGAppear.FB$Start),
                 j*100)
    out <- rbind(out,rates)
  }
  out <- rbind(out,rates) ##?OUR_COMMENT:: This line is responsible for the "bug" described
  # a few lines below. It could be removed, but let's stick with
  # the original script.
}
colnames(out)<-c("NGA","PreRate","PostRate","MGUncertainty","TimeWindow")

write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE) #Exporting/importing to force it to read
#as numeric (there is probably a more
#elegant way to do this)

out<-read.table("./temp/FullRates.csv", sep=",", header=TRUE, stringsAsFactors = TRUE)

out$Difference<-out[,3]-out[,2]


for(i in 1:length(out[,5])){
  ##?OUR_COMMENT:: First, fix the deleted first row discussed on lines 163:170
  if(i == 1){
    out[i,7]<- -9900
  }
  else if(i >1){
    out[i,7]<-out[i,5]-out[i-1,5]}}
out <-subset(out, out[,7]!=0) #getting rid of bug when the final row repeats in each NGA

write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE)

out<-read.table("./temp/FullRates.csv", sep=",", header=TRUE, stringsAsFactors = TRUE)

out <-subset(out, out[,5]<750) #Change this to modify time-window restriction from 700 years
#pre/post moralizing gods (<750) or 2050 out to use full time-window
write.csv(out, file="./temp/EqualRates_uncq_FB_300.csv",  row.names=FALSE)


##?OUR_COMMENT:: Result for Pre/Post-MG with MG appearance shifted
#  300 years back in time.
print("t-test for 6 NGAs where MG emerge at t - 300 for time span +-700 years")
print(t.test(out[,3], out[,2],paired=TRUE))

print("SC rate of change Pre/Post-MG")
print(mean(out[,2], na.rm = T)/mean(out[,3],na.rm = T))
print("########################################################################")
##?OUR_COMMENT:: How many Pre-/Post-MG data per NGA?
out.na <- out[!is.na(out$Difference),]
out.na %>% group_by(NGA) %>% summarise(length(Difference))

dir_init("./output")

files <- list.files("./temp", full.names = TRUE)
file.copy(files, "./output")

