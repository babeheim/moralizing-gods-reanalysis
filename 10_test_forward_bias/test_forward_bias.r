
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
  polities <- read.csv("./input/polities.csv", header = TRUE)
  
  #New scripts for automated analysis of rates of change in social complexity pre/post
  # moralising gods/doctrinal mode/writing
  
  #dat <- read.table("PC1_traj_merged.csv", sep=",", header=TRUE) #?? everything is here
  dat <- read.csv("./input/PC1_traj_merged.csv", header = TRUE)
  
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

out<-read.table("./temp/FullRates.csv", sep=",", header=TRUE)

out$Difference<-out[,2]-out[,3]
write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#### 1.3 Pre-Post MG by century ####

##?OUR_COMMENT:: We have to add World.Region and language Family to this data set, 
# so we can use it later in the models.

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
out<-read.table("./temp/FullRates.csv", sep=",", header=TRUE)

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
  # be missing data for this century, the row shoul be there containing all "NA" rather than missing.
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
print("Original t-test for MG = 0 with +/-2000 time-span")
print(t.test(out[,3], out[,2],paired=TRUE))

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

out<-read.table("./temp/FullRates.csv", sep=",", header=TRUE)

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
print("Original t-test for MG = 0 with +/-700 time-span")
print(t.test(out[,3], out[,2],paired=TRUE))

mean(out[,2], na.rm = T)/mean(out[,3],na.rm = T)

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
## 2. Forward-bias correction t-100 ####

##?OUR_COMMENT:: From now on, we present our own re-analysis of the original script.

##?OUR_COMMENT:: In this section, we will keep up with the t-test,
# but try to correct for forward bias, i.e., the fact that if we take the first known
# occurrence of MGs, it is very unlikely that it is actually the oldest occurrence. We will also
# fix the errors mentioned in Section 1.


##?OUR_COMMENT:: Correct for missing centuries at Deccan and Paris Basin

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
### 2.1 Shift MG for t-100 ####

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
  
  out<-read.table("./temp/FullRates.csv", sep=",", header=TRUE)
  
  out$Difference<-out[,3]-out[,2]
  
  ##?OUR_COMMENT:: save the difference into separate variable needed for histogram comparisons;
  # note this variable automatically codes FB size
  assign(paste0("d",FB),out$Difference*1000 )
  
  
  for(i in 1:length(out[,5])){
    ##?OUR_COMMENT:: First, fix the deleted first row discussed on lines 163:170
    if(i == 1){
      out[i,7]<- -9900
    }
    else if(i >1){
    out[i,7]<-out[i,5]-out[i-1,5]}}
  out <-subset(out, out[,7]!=0) #getting rid of bug when the final row repeats in each NGA
  
  write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE)
  
  out<-read.table("./temp/FullRates.csv", sep=",", header=TRUE)
  
  out <-subset(out, out[,5]<2050) #Change this to modify time-window restriction from 700 years
  #pre/post moralizing gods (<750) or 2050 out to use full time-window
  write.csv(out, file="./temp/EqualRates_FB_100.csv",  row.names=FALSE)
  

##?OUR_COMMENT:: Result for +/- 2000 years Pre/Post-MG with MG appearance shifted
#  100 years back in time.
print("Corrected t-test for MG-100 with +/-2000 time-span")
print(t.test(out[,3], out[,2],paired=TRUE))

##?OUR_COMMENT:: This result shows that if we decrease each MG time by 100 years, MG positively
# predict SC.

mean(out[,3], na.rm = T)/mean(out[,2],na.rm = T)

##?OUR_COMMENT:: How many Pre-/Post-MG data per NGA?
out.na <- out[!is.na(out$Difference),]
out.na %>% group_by(NGA) %>% summarise(length(Difference))


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## 2.2 Main plot t-100 ####

##?OUR_COMMENT:: Here, we report the same plot as Whitehouse et al.'s Fig. 2a,
# only moving the assumed appearance of MGs 100 years back.

NGAs <- c("Deccan", "Kachi Plain", "Kansai", "Konya Plain", "Latium",
          "Middle Yellow River Valley", "Niger Inland Delta", "Orkhon Valley",
          "Paris Basin", "Sogdiana", "Susiana", "Upper Egypt")
dat.cora <- matrix(NA, nrow=0, ncol=5)

for(i in 1:length(NGAs)){
  dt <- out[out$NGA == NGAs[i],]
  ot <- dat.cor[dat.cor$NGA == NGAs[i],]
  
  ##?OUR_COMMENT:: correct for FB
  firstMG.row <- match(1, ot$MoralisingGods)  ##?OUR_COMMENT:: find first MG appearance
  ot$MoralisingGods[firstMG.row-(FB/100)] <- 1  ##?OUR_COMMENT:: shift first MG X lines back
  
  MG<-subset(ot,MoralisingGods=="1") #Replace "MoralisingGods" with "DoctrinalMode" or 
  #"Writing" to do these analyses
  #library(dplyr)
  #MG<-as.data.frame(MG %>% group_by(PolID) %>% sample_n(size = 1))
  #randomly samples so there is only one century per polity
  MGAppear<-subset(MG, Time==min(Time))
  assign(paste0("t",i),MGAppear$Time)
  DM<-subset(ot,DoctrinalMode=="1") #Replace "MoralisingGods" with "DoctrinalMode" or 
  #"Writing" to do these analyses
  DMAppear<-subset(DM, Time==min(Time))
  WR<-subset(ot,Writing=="1") #Replace "MoralisingGods" with "DoctrinalMode" or "Writing" 
  #to do these analyses
  WRAppear<-subset(WR, Time==min(Time))
  my.values<-c(NGAs[i],ifelse(class(mean(dt[,3],na.rm=TRUE))=="NULL","NA",
                              mean(dt[,3],na.rm=TRUE)),
               ifelse(class(mean(dt[,2],na.rm=TRUE))=="NULL","NA",
                              mean(dt[,2],na.rm=TRUE)),
               ifelse(class(std.error(dt[,3],na.rm=TRUE))=="NULL","NA",
                      1.96*std.error(dt[,3],na.rm=TRUE)),
               ifelse(class(std.error(dt[,2],na.rm=TRUE))=="NULL","NA",
                      1.96*std.error(dt[,2],na.rm=TRUE)))
  dat.cora <- rbind(dat.cora,my.values)
}

colnames(dat.cora)<-c("NGA","PostRate","PreRate","PostConfInt","PreConfInt")
write.csv(dat.cora, file="./temp/PrePostComparisonFull.csv",  row.names=FALSE)
dat.cora<-read.table("./temp/PrePostComparisonFull.csv", sep=",", header=TRUE)
dat.cora<-as.data.frame(dat.cora)
dat.cora$Difference<-dat.cora[,2]-dat.cora[,3]
dat.cora[,2:6]<-dat.cora[,2:6]*1000
write.csv(dat.cora, file="./temp/PrePostComparisonFull.csv",  row.names=FALSE)

#Full values for matching pre-/post-NGAs
out <-subset(out, out[,6]<1000) #Removing windows without matching pre-/post-MG rates
write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE)
out <- read.table("./temp/FullRates.csv", sep=",", header=TRUE)

dat.cora <- matrix(NA, nrow=0, ncol=15)

for(i in 1:length(NGAs)){
  dt <- out[out$NGA == NGAs[i],]
  ot <- dat.cor[dat.cor$NGA == NGAs[i],]
  
  ##?OUR_COMMENT:: correct for FB
  firstMG.row <- match(1, ot$MoralisingGods)
  ot$MoralisingGods[firstMG.row-(FB/100)] <- 1
  
  MG<-subset(ot,MoralisingGods=="1") #Replace "MoralisingGods" with "DoctrinalMode" or
  #"Writing" to do these analyses
  #library(dplyr)
  #MG<-as.data.frame(MG %>% group_by(PolID) %>% sample_n(size = 1))
  #randomly samples so there is only one century per polity
  MGAppear<-subset(MG, Time==min(Time))
  DM<-subset(ot,DoctrinalMode=="1") #Replace "MoralisingGods" with "DoctrinalMode" or
  #"Writing" to do these analyses
  DMAppear<-subset(DM, Time==min(Time))
  WR<-subset(ot,Writing=="1") #Replace "MoralisingGods" with "DoctrinalMode" or "Writing"
  #to do these analyses
  WRAppear<-subset(WR, Time==min(Time))
  my.values <- c(NGAs[i],
                 mean(dt[,3]),
                 mean(dt[,2]),
                 1.96*std.error(dt[,3]),
                 1.96*std.error(dt[,2]),
                 t.test(dt[,3],dt[,2])$p.value,
                 t.test(dt[,3],dt[,2])$parameter,
                 length(dt[,2]),
                 t.test(dt[,3],
                        dt[,2])$statistic,
                 DMAppear$Time-MGAppear$Time,
                 WRAppear$Time-MGAppear$Time,
                 MGAppear$End-MGAppear$Start,
                 DMAppear$End-DMAppear$Start,
                 WRAppear$End-WRAppear$Start,
                 MGAppear$Time)
  dat.cora <- rbind(dat.cora,my.values)
}

colnames(dat.cora)<-c("NGA","PostRate","PreRate","PostConfInt","PreConfInt",
                  "p","df","n","t","PreMGRitual","PreMGWriting","RangeMGAppear",
                  "RangeDMAppear","RangeWRAppear","MGAppear")
write.csv(dat.cora, file="./temp/PrePostComparison.csv",  row.names=FALSE)
dat.cora<-read.table("./temp/PrePostComparison.csv", sep=",", header=TRUE)
dat.cora<-as.data.frame(dat.cora)
dat.cora$Difference<-dat.cora[,2]-dat.cora[,3]
dat.cora[,c(2:5,16)]<-dat.cora[,c(2:5,16)]*1000
write.csv(dat.cora, file="./temp/PrePostComparison.csv",  row.names=FALSE)

######Normalize time-series centered around moralising god appearance

out <- matrix(NA, nrow=0, ncol=0)
for(i in 1:length(NGAs)){
  dt <- dat.cor[dat.cor$NGA == NGAs[i],]
  Earliest<-subset(dt,Time==min(Time))	
  Latest<-subset(dt,Time==max(Time))	
  
  ##?OUR_COMMENT:: correct for FB
  firstMG.row <- match(1, dt$MoralisingGods)
  dt$MoralisingGods[firstMG.row-(FB/100)] <- 1
  
  
  MG<-subset(dt,MoralisingGods=="1") #Replace "MoralisingGods" with "DoctrinalMode" or
  #"Writing" to do these analyses
  MGAppear<-subset(MG, Time==min(Time))
  dt$Time.norm<-dt$Time-(MGAppear$Time)
  out <- rbind(out,dt)
}

out$MGUncertainty<-out$End-out$Start
write.csv(out, file="./temp/TimeNorm.csv",  row.names=FALSE) 

#Merge Normalized times
dat.cor <- read.table("./temp/TimeNorm.csv", sep=",", header=TRUE)
out<-unique(dat.cor$Time.norm)

#library(dplyr) ##Bugs when I load dplyr and plyr (even when only loading dplyr after plyr)!

for(i in 1:length(NGAs)){
  dt <- dat.cor[dat.cor$NGA == NGAs[i],]
  out<-merge(out,dt[,c("Time.norm","Mean")],by.x="x",by.y="Time.norm",all=TRUE)
  out<-rename(out, c("Mean" = NGAs[i]))
}

MeanPCs<-out[,2:(1+length(NGAs))]
out$Mean <- apply(MeanPCs,1,mean,na.rm=TRUE)
out$Lower <- out$Mean - 1.96*apply(MeanPCs,1, std.error,na.rm=TRUE)
out$Upper <- out$Mean + 1.96*apply(MeanPCs,1, std.error,na.rm=TRUE)
write.csv(out, file="./temp/SCNorm.csv",  row.names=FALSE) 

{
  #plot normalized times
  FullImpdat.cor<-read.csv('./temp/SCNorm.csv', header=TRUE)
  
  
  #data <- FullImpdat.cor
  y <- out$Mean
  x <- out$x
  l <- out$Lower
  u <- out$Upper
  #earliest<-read.csv('EarliestAppearance.csv', header=TRUE)
  #a<-earliest[,5]
  #b<-earliest[,6]
  #c<-earliest[,2]
  #d<-earliest[,4]
  
  #ylim <- c(min(na.rm=TRUE,y), max(na.rm=TRUE,y)) #This sets y-axis from the
  # minimum to maximum values throughout the whole global dataset
  ylim <- c(0,1) #This sets x-axis to min (0) to max (1) social complexity
  #xlim <- c(min(na.rm=TRUE,x), max(na.rm=TRUE,x)) #Use this instead to set x-axis from earliest
  #to latest time in dataset
  xlim <- c(-2000,2000) #This sets x-axis to 1,000 years before and after the appearance of
  #moralising gods
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
  col2<-rgb(207,58,58,max=255,alpha=200)
  col3<-rgb(24,176,232,max=255,alpha=200)
}

png("./temp/Fig2A_MG100.png",width = 6,height = 5.5,units = 'in', res = 300)

plot(x, y, ylim=ylim, xlim=xlim, pch=pch, cex=cex, xaxt=xaxt,
     ann=ann, yaxt=yaxt,type=type,xaxs=xaxs,yaxs=yaxs)
panel.first = rect(
  c(mean(dat.cora$PreMGRitual)-(1.96*std.error(dat.cora$PreMGRitual)),
    0-mean(dat.cora$RangeMGAppear)/2), -1e6,
  c(mean(dat.cora$PreMGRitual)+ (1.96*std.error(dat.cora$PreMGRitual)),
    0+mean(dat.cora$RangeMGAppear)/2), 1e6,
  col=col2, border=NA)

#panel.first = rect(c(mean(dat.cora$PreMGRitual)-(1.96*std.error(dat.cora$PreMGRitual)),
#mean(dat.cora$PreMGWriting)-(1.96*std.error(dat.cora$PreMGWriting)), 0), -1e6, c(mean(dat.cora$PreMGRitual)+
#(1.96*std.error(dat.cora$PreMGRitual)), mean(dat.cora$PreMGWriting)+(1.96*std.error(dat.cora$PreMGWriting)),
#0+mean(dat.cora$RangeMGAppear)), 1e6, col=c(col1,col2,col3), border=NA)
#This was earlier version that included writing
#polygon(c(x, rev(x)) , c(u, rev(l)) , col = 'grey' , border = NA) #out for now because of bugs
lines(x, y,type="l") 
lines(x, u,type="l",lty="dotted") 
lines(x, l,type="l",lty="dotted")
#text(100, 0.2, "Moralizing gods",
#     cex = 1, srt=90)
#text(-1200, 0.1, "Doctrinal rituals",
#     cex = 1)

#abline(h=0.6,lty="dashed")
dev.off()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## 2.3 Histogram for t-100 ####

##?OUR_COMMENT:: Let's plot the histogram with the rate of SC change:

##?OUR_COMMENT:: select the size of FB
ifelse(FB == 300, d <- d300,d <- d100)
d <- as.data.frame(d)

ifelse(FB == 300, c <- "#18b0e8",c <- "#cf3a3a")


ggplot() + 
  geom_histogram(data = d, aes(x = d), fill = c, colour = c, binwidth = 0.1,
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
ggsave("./temp/Histogram_MG100.png",width = 4,height = 2.7, dpi = 300)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## 3. Forward-bias correction t-300 ####

##?OUR_COMMENT:: In this section, we will keep up with the t-test,
# but try to correct for forward bias, i.e., the fact that if we take the first known
# occurrence of MG, it is very unlikely that it is actually the oldest occurrence.

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
### 3.1 Shift MG for t-300 ####

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
      rates<-cbind(as.character(MGAppear$NGA),
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
  
  out<-read.table("./temp/FullRates.csv", sep=",", header=TRUE)
  
  out$Difference<-out[,3]-out[,2]
  
  ##?OUR_COMMENT:: save the difference into separate variable needed for histogram comparisons;
  # note this variable automatically codes FB size
  assign(paste0("d",FB),out$Difference*1000 )
  
  
  for(i in 1:length(out[,5])){
    ##?OUR_COMMENT:: First, fix the deleted first row discussed on lines 163:170
    if(i == 1){
      out[i,7]<- -9900
    }
    else if(i >1){
      out[i,7]<-out[i,5]-out[i-1,5]}}
  out <-subset(out, out[,7]!=0) #getting rid of bug when the final row repeats in each NGA
  
  write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE)
  
  out<-read.table("./temp/FullRates.csv", sep=",", header=TRUE)
  
  out <-subset(out, out[,5]<2050) #Change this to modify time-window restriction from 700 years
  #pre/post moralizing gods (<750) or 2050 out to use full time-window
  write.csv(out, file="./temp/EqualRates_FB_300.csv",  row.names=FALSE)
  
  
##?OUR_COMMENT:: Result for +/- 2000 years Pre/Post-MG with MG appearance shifted
#  100 years back in time.
print("Corrected t-test for MG-300 with +/-2000 time-span")
print(t.test(out[,3], out[,2],paired=TRUE))

mean(out[,3], na.rm = T)/mean(out[,2],na.rm = T)

##?OUR_COMMENT:: How many Pre-/Post-MG data per NGA?
out.na <- out[!is.na(out$Difference),]
out.na %>% group_by(NGA) %>% summarise(length(Difference))


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## 3.2 Main plot t-300 ####

##?OUR_COMMENT:: Here, we report the same plot as Whitehouse et al.'s Fig. 2a,
# only moving the assumed appearance of MGs 100 years back.

NGAs <- c("Deccan", "Kachi Plain", "Kansai", "Konya Plain", "Latium",
          "Middle Yellow River Valley", "Niger Inland Delta", "Orkhon Valley",
          "Paris Basin", "Sogdiana", "Susiana", "Upper Egypt")
dat.cora <- matrix(NA, nrow=0, ncol=5)

for(i in 1:length(NGAs)){
  dt <- out[out$NGA == NGAs[i],]
  ot <- dat.cor[dat.cor$NGA == NGAs[i],]
  
  ##?OUR_COMMENT:: correct for FB
  firstMG.row <- match(1, ot$MoralisingGods)  ##?OUR_COMMENT:: find first MG appearance
  ot$MoralisingGods[firstMG.row-(FB/100)] <- 1  ##?OUR_COMMENT:: shift first MG X lines back
  
  MG<-subset(ot,MoralisingGods=="1") #Replace "MoralisingGods" with "DoctrinalMode" or 
  #"Writing" to do these analyses
  #library(dplyr)
  #MG<-as.data.frame(MG %>% group_by(PolID) %>% sample_n(size = 1))
  #randomly samples so there is only one century per polity
  MGAppear<-subset(MG, Time==min(Time))
  assign(paste0("t",i),MGAppear$Time)
  DM<-subset(ot,DoctrinalMode=="1") #Replace "MoralisingGods" with "DoctrinalMode" or 
  #"Writing" to do these analyses
  DMAppear<-subset(DM, Time==min(Time))
  WR<-subset(ot,Writing=="1") #Replace "MoralisingGods" with "DoctrinalMode" or "Writing" 
  #to do these analyses
  WRAppear<-subset(WR, Time==min(Time))
  my.values<-c(NGAs[i],ifelse(class(mean(dt[,3],na.rm=TRUE))=="NULL","NA",
                              mean(dt[,3],na.rm=TRUE)),
               ifelse(class(mean(dt[,2],na.rm=TRUE))=="NULL","NA",
                      mean(dt[,2],na.rm=TRUE)),
               ifelse(class(std.error(dt[,3],na.rm=TRUE))=="NULL","NA",
                      1.96*std.error(dt[,3],na.rm=TRUE)),
               ifelse(class(std.error(dt[,2],na.rm=TRUE))=="NULL","NA",
                      1.96*std.error(dt[,2],na.rm=TRUE)))
  dat.cora <- rbind(dat.cora,my.values)
}

colnames(dat.cora)<-c("NGA","PostRate","PreRate","PostConfInt","PreConfInt")
write.csv(dat.cora, file="./temp/PrePostComparisonFull.csv",  row.names=FALSE)
dat.cora<-read.table("./temp/PrePostComparisonFull.csv", sep=",", header=TRUE)
dat.cora<-as.data.frame(dat.cora)
dat.cora$Difference<-dat.cora[,2]-dat.cora[,3]
dat.cora[,2:6]<-dat.cora[,2:6]*1000
write.csv(dat.cora, file="./temp/PrePostComparisonFull.csv",  row.names=FALSE)

#Full values for matching pre-/post-NGAs
out <-subset(out, out[,6]<1000) #Removing windows without matching pre-/post-MG rates
write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE)
out <- read.table("./temp/FullRates.csv", sep=",", header=TRUE)

dat.cora <- matrix(NA, nrow=0, ncol=15)

for(i in 1:length(NGAs)){
  dt <- out[out$NGA == NGAs[i],]
  ot <- dat.cor[dat.cor$NGA == NGAs[i],]
  
  ##?OUR_COMMENT:: correct for FB
  firstMG.row <- match(1, ot$MoralisingGods)
  ot$MoralisingGods[firstMG.row-(FB/100)] <- 1
  
  MG<-subset(ot,MoralisingGods=="1") #Replace "MoralisingGods" with "DoctrinalMode" or
  #"Writing" to do these analyses
  #library(dplyr)
  #MG<-as.data.frame(MG %>% group_by(PolID) %>% sample_n(size = 1))
  #randomly samples so there is only one century per polity
  MGAppear<-subset(MG, Time==min(Time))
  DM<-subset(ot,DoctrinalMode=="1") #Replace "MoralisingGods" with "DoctrinalMode" or
  #"Writing" to do these analyses
  DMAppear<-subset(DM, Time==min(Time))
  WR<-subset(ot,Writing=="1") #Replace "MoralisingGods" with "DoctrinalMode" or "Writing"
  #to do these analyses
  WRAppear<-subset(WR, Time==min(Time))
  my.values <- c(NGAs[i],
                 mean(dt[,3]),
                 mean(dt[,2]),
                 1.96*std.error(dt[,3]),
                 1.96*std.error(dt[,2]),
                 t.test(dt[,3],dt[,2])$p.value,
                 t.test(dt[,3],dt[,2])$parameter,
                 length(dt[,2]),
                 t.test(dt[,3],
                        dt[,2])$statistic,
                 DMAppear$Time-MGAppear$Time,
                 WRAppear$Time-MGAppear$Time,
                 MGAppear$End-MGAppear$Start,
                 DMAppear$End-DMAppear$Start,
                 WRAppear$End-WRAppear$Start,
                 MGAppear$Time)
  dat.cora <- rbind(dat.cora,my.values)
}

colnames(dat.cora)<-c("NGA","PostRate","PreRate","PostConfInt","PreConfInt",
                  "p","df","n","t","PreMGRitual","PreMGWriting","RangeMGAppear",
                  "RangeDMAppear","RangeWRAppear","MGAppear")
write.csv(dat.cora, file="./temp/PrePostComparison.csv",  row.names=FALSE)
dat.cora<-read.table("./temp/PrePostComparison.csv", sep=",", header=TRUE)
dat.cora<-as.data.frame(dat.cora)
dat.cora$Difference<-dat.cora[,2]-dat.cora[,3]
dat.cora[,c(2:5,16)]<-dat.cora[,c(2:5,16)]*1000
write.csv(dat.cora, file="./temp/PrePostComparison.csv",  row.names=FALSE)

######Normalize time-series centered around moralising god appearance

out <- matrix(NA, nrow=0, ncol=0)
for(i in 1:length(NGAs)){
  dt <- dat.cor[dat.cor$NGA == NGAs[i],]
  Earliest<-subset(dt,Time==min(Time))	
  Latest<-subset(dt,Time==max(Time))	
  
  ##?OUR_COMMENT:: correct for FB
  firstMG.row <- match(1, dt$MoralisingGods)
  dt$MoralisingGods[firstMG.row-(FB/100)] <- 1
  
  
  MG<-subset(dt,MoralisingGods=="1") #Replace "MoralisingGods" with "DoctrinalMode" or
  #"Writing" to do these analyses
  MGAppear<-subset(MG, Time==min(Time))
  dt$Time.norm<-dt$Time-(MGAppear$Time)
  out <- rbind(out,dt)
}

out$MGUncertainty<-out$End-out$Start
write.csv(out, file="./temp/TimeNorm.csv",  row.names=FALSE) 

#Merge Normalized times
dat.cor <- read.table("./temp/TimeNorm.csv", sep=",", header=TRUE)
out<-unique(dat.cor$Time.norm)

#library(dplyr) ##Bugs when I load dplyr and plyr (even when only loading dplyr after plyr)!

for(i in 1:length(NGAs)){
  dt <- dat.cor[dat.cor$NGA == NGAs[i],]
  out<-merge(out,dt[,c("Time.norm","Mean")],by.x="x",by.y="Time.norm",all=TRUE)
  out<-rename(out, c("Mean" = NGAs[i]))
}

MeanPCs<-out[,2:(1+length(NGAs))]
out$Mean <- apply(MeanPCs,1,mean,na.rm=TRUE)
out$Lower <- out$Mean - 1.96*apply(MeanPCs,1, std.error,na.rm=TRUE)
out$Upper <- out$Mean + 1.96*apply(MeanPCs,1, std.error,na.rm=TRUE)
write.csv(out, file="./temp/SCNorm.csv",  row.names=FALSE) 

{
  #plot normalized times
  FullImpdat.cor<-read.csv('./temp/SCNorm.csv', header=TRUE)
  
  
  #data <- FullImpdat.cor
  y <- out$Mean
  x <- out$x
  l <- out$Lower
  u <- out$Upper
  #earliest<-read.csv('EarliestAppearance.csv', header=TRUE)
  #a<-earliest[,5]
  #b<-earliest[,6]
  #c<-earliest[,2]
  #d<-earliest[,4]
  
  #ylim <- c(min(na.rm=TRUE,y), max(na.rm=TRUE,y)) #This sets y-axis from the
  # minimum to maximum values throughout the whole global datset
  ylim <- c(0,1) #This sets x-axis to min (0) to max (1) social complexity
  #xlim <- c(min(na.rm=TRUE,x), max(na.rm=TRUE,x)) #Use this instead to set x-axis from earliest
  #to latest time in dataset
  xlim <- c(-2000,2000) #This sets x-axis to 1,000 years before and after the appearance of
  #moralising gods
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
  col2<-rgb(207,58,58,max=255,alpha=200)
  col3<-rgb(24,176,232,max=255,alpha=200)
}
png("./temp/Fig2A_MG300.png",width = 6,height = 5.5,units = 'in', res = 300)

##?OUR_COMMENT:: Had to take the doctrinal rituals shading out, it coverd the shifted MGs

plot(x, y, ylim=ylim, xlim=xlim, pch=pch, cex=cex, xaxt=xaxt,
     ann=ann, yaxt=yaxt,type=type,xaxs=xaxs,yaxs=yaxs)
panel.first = rect(
  c(0-mean(dat.cora$RangeMGAppear)/2), -1e6,
  c(0+mean(dat.cora$RangeMGAppear)/2), 1e6,
  col=col3, border=NA)

#panel.first = rect(c(mean(dat.cora$PreMGRitual)-(1.96*std.error(dat.cora$PreMGRitual)),
#mean(dat.cora$PreMGWriting)-(1.96*std.error(dat.cora$PreMGWriting)), 0), -1e6, c(mean(dat.cora$PreMGRitual)+
#(1.96*std.error(dat.cora$PreMGRitual)), mean(dat.cora$PreMGWriting)+(1.96*std.error(dat.cora$PreMGWriting)),
#0+mean(dat.cora$RangeMGAppear)), 1e6, col=c(col1,col2,col3), border=NA)
#This was earlier version that included writing
#polygon(c(x, rev(x)) , c(u, rev(l)) , col = 'grey' , border = NA) #out for now because of bugs
lines(x, y,type="l") 
lines(x, u,type="l",lty="dotted") 
lines(x, l,type="l",lty="dotted")
#text(100, 0.2, "Moralizing gods",
#     cex = 1, srt=90)
#text(-1200, 0.1, "Doctrinal rituals",
#     cex = 1)
#abline(h=0.6,lty="dashed")
dev.off()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
### 3.3 Histogram for t-300####


##?OUR_COMMENT:: Let's plot the histogram with the rate of SC change:

##?OUR_COMMENT:: select the size of FB
ifelse(FB == 300, d <- d300,d <- d100)
d <- as.data.frame(d)

ifelse(FB == 300, c <- "#18b0e8",c <- "#cf3a3a")


ggplot() + 
  geom_histogram(data = d, aes(x = d), fill = c, colour = c, binwidth = 0.1,
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
ggsave("./temp/Histogram_MG300.png",width = 4,height = 2.7, dpi = 300)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## 4. Excluding conquered NGAs ####

##?OUR_COMMENT:: In this section, we will check the robustness of Whitehouse et al. results
# by excluding NGAs that acquited MGs trhough conquest by larger empires (Deccan, Kachi Plain,
# Sogdiana).

##?OUR_COMMENT:: Correct for missing centuries at Deccan and Paris Basin

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
### 4.1 Uncoquered MG t = 0 ####

NGAs.unconq <- c("Kansai", "Konya Plain", "Latium", "Paris Basin", 
                 "Middle Yellow River Valley", "Niger Inland Delta", "Orkhon Valley",
                 "Susiana", "Upper Egypt")

out <- matrix(NA, nrow=0, ncol=5)
for(i in 1:length(NGAs.unconq)){
  dt <- dat.cor[dat.cor$NGA == NGAs.unconq[i],]
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

out<-read.table("./temp/FullRates.csv", sep=",", header=TRUE)

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

print("Unconquered NGAs t-test for MG = 0 with +/-2000 time-span")
print(t.test(out[,3], out[,2],paired=TRUE))

mean(out[,2], na.rm = T)/mean(out[,3],na.rm = T)

##?OUR_COMMENT:: How many rows of Pre-/Post-MG data?
NROW(out[(!is.na(out[,2]+out[,3])),])

##?OUR_COMMENT:: How many Pre-/Post-MG data per NGA?
out.na <- out[!is.na(out$Difference),]
out.na %>% group_by(NGA) %>% summarise(length(Difference))

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
### 4.2 Uncoquered MG t-100 ####
NGAs <- c("Kansai", "Konya Plain", "Latium", "Paris Basin", 
          "Middle Yellow River Valley", "Niger Inland Delta", "Orkhon Valley",
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

out<-read.table("./temp/FullRates.csv", sep=",", header=TRUE)

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

out<-read.table("./temp/FullRates.csv", sep=",", header=TRUE)

out <-subset(out, out[,5]<2050) #Change this to modify time-window restriction from 700 years
#pre/post moralizing gods (<750) or 2050 out to use full time-window
write.csv(out, file="./temp/EqualRates_uncq_FB_100.csv",  row.names=FALSE)


##?OUR_COMMENT:: Result for +/- 2000 years Pre/Post-MG with MG appearance shifted
#  100 years back in time.
print("Unconquered NGAs t-test for MG - 100 with +/-2000 time-span")
print(t.test(out[,3], out[,2],paired=TRUE))


mean(out[,2], na.rm = T)/mean(out[,3],na.rm = T)

##?OUR_COMMENT:: How many rows of Pre-/Post-MG data?
NROW(out[(!is.na(out[,2]+out[,3])),])

##?OUR_COMMENT:: How many Pre-/Post-MG data per NGA?
out.na <- out[!is.na(out$Difference),]
out.na %>% group_by(NGA) %>% summarise(length(Difference))

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
### 4.3 Uncoquered MG t-300 ####
NGAs <- c("Kansai", "Konya Plain", "Latium", "Paris Basin", 
          "Middle Yellow River Valley", "Niger Inland Delta", "Orkhon Valley",
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

out<-read.table("./temp/FullRates.csv", sep=",", header=TRUE)

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

out<-read.table("./temp/FullRates.csv", sep=",", header=TRUE)

out <-subset(out, out[,5]<2050) #Change this to modify time-window restriction from 700 years
#pre/post moralizing gods (<750) or 2050 out to use full time-window
write.csv(out, file="./temp/EqualRates_uncq_FB_300.csv",  row.names=FALSE)


##?OUR_COMMENT:: Result for +/- 2000 years Pre/Post-MG with MG appearance shifted
#  100 years back in time.
print("Unconquered NGAs t-test for MG - 300 with +/-2000 time-span")
print(t.test(out[,3], out[,2],paired=TRUE))

mean(out[,3], na.rm = T)/mean(out[,2],na.rm = T)

##?OUR_COMMENT:: How many Pre-/Post-MG data per NGA?
out.na <- out[!is.na(out$Difference),]
out.na %>% group_by(NGA) %>% summarise(length(Difference))

dir_init("./output")

files <- list.files("./temp", full.names = TRUE)
file.copy(files, "./output")

