
# Note that comments in the text starting with #??OUR_COMMENT are our new comments.
# We left all the all comments in the text as well. 

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
### 1. Load libraries and data ####


# This section loads data on polities, their social complexity (SC),
# and the presence moralizing gods (MG).

{
  rm(list = ls())
  
  source("../project_support.r")
  
  dir_init("./hierarchical_models_output")
polities <- read.csv("./input/polities.csv", header = TRUE)

#New scripts for automated analysis of rates of change in social complexity pre/post
# moralising gods/doctrinal mode/writing

#dat <- read.table("PC1_traj_merged.csv", sep=",", header=TRUE) #?? everything is here
dat <- read.csv("./input/PC1_traj_merged.csv", stringsAsFactors = FALSE)

dat$NGA<-as.character(dat$NGA)
NGAs <- levels(polities$NGA)
NGAs <- NGAs[NGAs != "Crete"]    #### Remove new NGAs
NGAs <- NGAs[NGAs != "Galilee"]

##?OUR_COMMENT:: We have to include these variables for later use in the multi-level
# models (they will serve as nesting factors) 
dat$World.Region <- as.factor(dat$World.Region)
dat$Family <- as.factor(dat$Family)
}

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
### 1.1. MG data for +/- 2000 SC analysis ####

{
  dat.MG <- dat
  NGAs <- c("Deccan", "Kachi Plain", "Kansai", "Konya Plain", "Latium",
            "Middle Yellow River Valley", "Niger Inland Delta", "Orkhon Valley",
            "Paris Basin", "Sogdiana", "Susiana", "Upper Egypt")
  dat.MG <- dat.MG[dat.MG$NGA %in% NGAs,]
  
  ##?OUR_COMMENT:: turn the absence of evidence into the evidence of absence
  dat.MG$MoralisingGods[is.na(dat.MG$MoralisingGods)] <- 0
  dat.MG$MoralisingGods <- as.factor(dat.MG$MoralisingGods)
  
  
  ##?OUR_COMMENT:: Standardize Time for each NGA
  
  for(i in 1:length(NGAs)){
    dt <- dat.MG[dat.MG$NGA == NGAs[i],]
    mt <- min(dt$Time)
    if(mt >= 0){
      dat.MG$Time[dat.MG$NGA == NGAs[i]] <- dat.MG$Time[dat.MG$NGA == NGAs[i]]-mt}
    else if(mt < 0){ 
      dat.MG$Time[dat.MG$NGA == NGAs[i]] <- dat.MG$Time[dat.MG$NGA == NGAs[i]]+abs(mt)
    }
  }
  
  
  ##?OUR_COMMENT::
  # We have to censor the data to the +/- 2000 years as in the original paper.
  
  # Before doing that, however, we have to fix what appears to be a mistake at the Deccan
  # and Paris Basin NGAs. Deccan misses two centuries of data (i.e., there are no rows for these
  # centuries rather than having a row with NAs) and Paris Basin misses one century.
  # We have to recode the later centuries to get rid of the missing centuries,
  # which would be problematic in the later analysis.
  
  dat.MG$Time[dat.MG$NGA == "Deccan" & (dat.MG$Time > 3000)] <- 
    dat.MG$Time[dat.MG$NGA == "Deccan" & (dat.MG$Time > 3000)] - 100
  dat.MG$Time[dat.MG$NGA == "Deccan" & (dat.MG$Time > 3600)] <- 
    dat.MG$Time[dat.MG$NGA == "Deccan" & (dat.MG$Time > 3600)] - 100
  
  dat.MG$Time[dat.MG$NGA == "Paris Basin" & (dat.MG$Time > 3500)] <- 
    dat.MG$Time[dat.MG$NGA == "Paris Basin" & (dat.MG$Time > 3500)] - 100
  
  
  ##?OUR_COMMENT:: Ok, we can continue
  censor = 2050
  
  for(i in 1:length(NGAs)){
    dt <- dat.MG[dat.MG$NGA == NGAs[i],]
    MG<-subset(dt,MoralisingGods=="1")
    MGAppear<-subset(MG, Time==min(Time))
    t <- MGAppear$Time
    t1 <- t-censor
    t2 <- t+censor-100  ##?OUR_COMMENT:: Note that in order to have 20 centuries Pre- and Post-MG,
    # we have to extract 0-19 centuries for Pre-MG and 20-39 centuries
    # for Post-MG (rather than 20-40, which would be 21 centuries)
    if (min(dat.MG$Time[dat.MG$NGA == NGAs[i]], na.rm = T) < t1){
      dat.MG$Time[dat.MG$NGA == NGAs[i] & (dat.MG$Time < t1)] <- NA}
    if (max(dat.MG$Time[dat.MG$NGA == NGAs[i]], na.rm = T) > t2){ 
      dat.MG$Time[dat.MG$NGA == NGAs[i] & (dat.MG$Time > t2)] <- NA}
  }
  dat.MG <- dat.MG[!is.na(dat.MG$Time),]
  
  ##?OUR_COMMENT:: Make Time series at each NGA to start with 0
  for(i in 1:length(NGAs)){
    dt <- dat.MG[dat.MG$NGA == NGAs[i],]
    m <- max(dt$Time, na.rm=T)
    dat.MG$Time[dat.MG$NGA == NGAs[i]] <- dat.MG$Time[dat.MG$NGA == NGAs[i]] - (m+100)
    dt <- dat.MG[dat.MG$NGA == NGAs[i],]
    m <- min(dt$Time, na.rm=T)
    dat.MG$Time[dat.MG$NGA == NGAs[i]] <- dat.MG$Time[dat.MG$NGA == NGAs[i]] + abs(m)}
  
  
  ##?OUR_COMMENT:: We have to impute new rows now, so each NGA would have MG == 1 at Time == 2000.
  # We will not make this step in the +/- 700 years analysis because all NGAs have data
  # for that time-range.
  
  # Deccan
  New <- dat.MG[1,]
  New[,] <- NA 
  dat.1 <- dat.MG[1:39,]
  dat.2 <- dat.MG[40:nrow(dat.MG),]
  dat.MG <- rbind(dat.1, New, dat.2)
  dat.MG[40,1] <- dat.MG[39,1]
  dat.MG[40,4] <- dat.MG[39,4] + 100
  
  # Kansai
  New <- dat.MG[1:11,]
  New[1:11,] <- NA 
  dat.1 <- dat.MG[1:80,]
  dat.2 <- dat.MG[81:nrow(dat.MG),]
  dat.MG <- rbind(dat.1, New, dat.2)
  dat.MG[81:91,1] <- dat.2[1,1]
  
  New <- dat.MG[1:7,]
  New[1:7,] <- NA 
  dat.1 <- dat.MG[1:113,]
  dat.2 <- dat.MG[114:nrow(dat.MG),]
  dat.MG <- rbind(dat.1, New, dat.2)
  dat.MG[114:120,1] <- dat.1[113,1]
  
  dat.MG$Time[dat.MG$NGA=="Kansai"] <- seq(0,3900,100)
  
  
  # Niger Inland Delta
  New <- dat.MG[1:7,]
  New[1:7,] <- NA 
  dat.1 <- dat.MG[1:240,]
  dat.2 <- dat.MG[241:nrow(dat.MG),]
  dat.MG <- rbind(dat.1, New, dat.2)
  dat.MG[241:248,1] <- dat.2[1,1]
  
  
  New <- dat.MG[1:12,]
  New[1:12,] <- NA 
  dat.1 <- dat.MG[1:268,]
  dat.2 <- dat.MG[269:nrow(dat.MG),]
  dat.MG <- rbind(dat.1, New, dat.2)
  dat.MG[269:280,1] <- dat.1[268,1]
  
  dat.MG$Time[dat.MG$NGA=="Niger Inland Delta"] <- seq(0,3900,100)
  
  #Orkhon Valley
  
  New <- dat.MG[1,]
  New[,] <- NA 
  dat.1 <- dat.MG[1:280,]
  dat.2 <- dat.MG[281:nrow(dat.MG),]
  dat.MG <- rbind(dat.1, New, dat.2)
  dat.MG[281,1] <- dat.MG[282,1]
  dat.MG[281,4] <- 0
  
  
  New <- dat.MG[1:5,]
  New[1:5,] <- NA 
  dat.1 <- dat.MG[1:315,]
  dat.2 <- dat.MG[316:nrow(dat.MG),]
  dat.MG <- rbind(dat.1, New, dat.2)
  dat.MG[316:320,1] <- dat.1[315,1]
  
  dat.MG$Time[dat.MG$NGA=="Orkhon Valley"] <- seq(0,3900,100)
  
  # Paris Basin
  New <- dat.MG[1:2,]
  New[1:2,] <- NA 
  dat.1 <- dat.MG[1:358,]
  dat.2 <- dat.MG[359:nrow(dat.MG),]
  dat.MG <- rbind(dat.1, New, dat.2)
  dat.MG[359:360,1] <- dat.1[358,1]
  dat.MG[359,4] <- 3800
  dat.MG[360,4] <- 3900
  
  # Upper Egypt
  New <- dat.MG[1:5,]
  New[1:5,] <- NA 
  dat.1 <- dat.MG[1:440,]
  dat.2 <- dat.MG[441:nrow(dat.MG),]
  dat.MG <- rbind(dat.1, New, dat.2)
  dat.MG[441:445,1] <- dat.2[1,1]
  
  dat.MG$Time[dat.MG$NGA=="Upper Egypt"] <- seq(0,3900,100)
  


##?OUR_COMMENT::OK, now all sites acquire MG at year 2000, so our data have 20 centuries before
# and 20 centuries after MG, although we miss data from some sites where there are not so
# long records. This should be taken care of in our multi-level models.
# As a consequence, we can use more data than the t-test model (which had to discard data not
# having both Pre- and Post-MG SC change for specific time-range).

##?OUR_COMMENT:: Now, let's assume that MGs were present at all sites after their first detection
#(assumption of the t-test analysis presented in 'forward_bias_correction.R')

dat.MG$MoralisingGods[dat.MG$Time>1900] <- 1

dat.MG$NGA <- as.factor(dat.MG$NGA)

##?OUR_COMMENT:: Save this data set for later use
dat.MG.2000 <- dat.MG
}

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
### 1.2. MG data for +/- 700 SC analysis ####

{
dat.MG <- dat
NGAs <- c("Deccan", "Kachi Plain", "Kansai", "Konya Plain", "Latium",
          "Middle Yellow River Valley", "Niger Inland Delta", "Orkhon Valley",
          "Paris Basin", "Sogdiana", "Susiana", "Upper Egypt")
dat.MG <- dat.MG[dat.MG$NGA %in% NGAs,]

##?OUR_COMMENT:: turn the absence of evidence into the evidence of absence
dat.MG$MoralisingGods[is.na(dat.MG$MoralisingGods)] <- 0
dat.MG$MoralisingGods <- as.factor(dat.MG$MoralisingGods)


##?OUR_COMMENT:: Standardize Time for each NGA to start with 0
# (no matter whether the dates are BCE or CE)

for(i in 1:length(NGAs)){
  dt <- dat.MG[dat.MG$NGA == NGAs[i],]
  mt <- min(dt$Time)
  if(mt >= 0){
    dat.MG$Time[dat.MG$NGA == NGAs[i]] <- dat.MG$Time[dat.MG$NGA == NGAs[i]]-mt}
  else if(mt < 0){ 
    dat.MG$Time[dat.MG$NGA == NGAs[i]] <- dat.MG$Time[dat.MG$NGA == NGAs[i]]+abs(mt)
  }
}

dat.MG$NGA <- as.factor(dat.MG$NGA)

##?OUR_COMMENT::
# We have to censor the data to the +/- 700 years interval Pre- and Post-MG
# as in the original t-test analysis.

# Before doing that, however, we have to fix what appears to be a mistake at the Deccan
# and Paris Basin NGAs. Deccan misses two centuries of data (i.e., there are no rows for
# these centuries rather than having a row with NAs) and Paris Basin misses one century.
# We have to recode the later centuries to get rid of the missing centuries, which would
# be problematic in the later analysis.

dat.MG$Time[dat.MG$NGA == "Deccan" & (dat.MG$Time > 3000)] <- 
  dat.MG$Time[dat.MG$NGA == "Deccan" & (dat.MG$Time > 3000)] - 100
dat.MG$Time[dat.MG$NGA == "Deccan" & (dat.MG$Time > 3600)] <- 
  dat.MG$Time[dat.MG$NGA == "Deccan" & (dat.MG$Time > 3600)] - 100

dat.MG$Time[dat.MG$NGA == "Paris Basin" & (dat.MG$Time > 3500)] <- 
  dat.MG$Time[dat.MG$NGA == "Paris Basin" & (dat.MG$Time > 3500)] - 100


##?OUR_COMMENT:: Ok, we can continue
censor = 750 ##?OUR_COMMENT:: censor the data to +/- 7 centuries

for(i in 1:length(NGAs)){
  dt <- dat.MG[dat.MG$NGA == NGAs[i],]
  MG<-subset(dt,MoralisingGods=="1")
  MGAppear<-subset(MG, Time==min(Time))
  t <- MGAppear$Time
  t1 <- t-censor
  t2 <- t+censor-100 ##?OUR_COMMENT:: Note that to have 7 centuries Pre- and Post-MG,
  # we have to extract 0-6 centuries for Pre-MG and 7-13 centuries
  # for Post-MG (rather than 7-14, which would be 8 centuries)
  if (min(dat.MG$Time[dat.MG$NGA == NGAs[i]], na.rm = T) < t1){
    dat.MG$Time[dat.MG$NGA == NGAs[i] & (dat.MG$Time < t1)] <- NA}
  if (max(dat.MG$Time[dat.MG$NGA == NGAs[i]], na.rm = T) > t2){ 
    dat.MG$Time[dat.MG$NGA == NGAs[i] & (dat.MG$Time > t2)] <- NA}
}
dat.MG <- dat.MG[!is.na(dat.MG$Time),]

##?OUR_COMMENT:: Make Time series at each NGA to start with 0 and end with 1300
# (14 centuries of data = +/- 700 years)

for(i in 1:length(NGAs)){
  dt <- dat.MG[dat.MG$NGA == NGAs[i],]
  m <- max(dt$Time, na.rm=T)
  dat.MG$Time[dat.MG$NGA == NGAs[i]] <- dat.MG$Time[dat.MG$NGA == NGAs[i]] - (m+100)
  dat.MG$Time[dat.MG$NGA == NGAs[i]] <- dat.MG$Time[dat.MG$NGA == NGAs[i]] + 1400}

##?OUR_COMMENT:: Now all sites acquire MG at year 700, so our data have 7 centuries
#before and 7 centuries after MG

##?OUR_COMMENT:: Now, let's assume that MGs were present at all sites after their first
# detection (assumption of the t-test analysis presented in 'forward_bias_correction.R')

dat.MG$MoralisingGods[dat.MG$Time>600] <- 1

dat.MG.700 <- dat.MG
}


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## 2. SC growth plot ####

### 2.1. SC by MGs +/- 700 ####

##?OUR_COMMENT:: This plot is based on Whitehouse et al. raw data. 
# We used discrete data (data are sampled by centuries) and plotted mean data with
# standard errors for each century.
# However, please keep in mind that these are just raw data, not taking into account the
# various nesting effects. We will explore those below.

dat.p <- dat.MG.700
  
ggplot(dat.p, aes(Time, Mean, color=MoralisingGods)) + 
  stat_summary(fun.data=mean_se, geom="pointrange") + 
  scale_color_manual(values = alpha(c("coral1","aquamarine3"), .8), labels = c("NA", "Present"),
                     name = "Moralizing Gods") + 
  scale_x_continuous(limits = c(0,1400),
                     breaks = c(0,200,400,600,700,800,1000,1200,1400),
                     labels = c("-700","-500","-300","-100","0","100","300","500","700")) +
  geom_vline(xintercept = 650,
             color = "grey", size=9, alpha = 0.3) + 
  annotate("text", label = "", x = 650, y = 0.35, size = 4, colour = "black", angle = 90) +
  labs(y="Social Complexity", x="Time (years before/after MG)") + 
  ylim(c(0.1,1)) +
  ggtitle("") +
  theme_bw() + 
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = rel(1.5)),        
    axis.line = element_line(colour = "black"),
    legend.position = c(0.4,0.95),
    legend.justification = c("right", "top"),
    legend.key.size = unit(0.8, "cm"),
    legend.title = element_text(size = rel(1.5)),
    legend.text = element_text(size = rel(1.2)),
    axis.title = element_text(size = rel(1.5)),
    axis.text.y= element_text(size = rel(1.5)),
    axis.text.x= element_text(size = rel(1.5)),
    plot.margin=unit(c(1,1,1,1),"cm")) 

##?OUR_COMMENT:: Save the plot if needed
ggsave("./hierarchical_models_output/PrePost_MG.png",width = 6,height = 6, dpi = 300)

##?OUR_COMMENT:: Compute the between-century rate of SC change before MGs
MG.change <- ((mean(dat.MG.700$Mean[dat.MG.700$Time == (700)], na.rm = T)/
                 mean(dat.MG.700$Mean[dat.MG.700$Time == (600)],na.rm = T))-1)*100

change <- matrix(NA,6,1)
for(i in 1:6){
  change[i] <- ((mean(dat.MG.700$Mean[dat.MG.700$Time == (700-i*100)], na.rm = T)/
                   mean(dat.MG.700$Mean[dat.MG.700$Time == (600-i*100)],na.rm = T))-1)*100
}

cbind(MG.change,mean(change, na.rm = T))


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
### 2.2. Unconcquered SC by MGs +/- 700 ####

NGAs.unconq <- c("Kansai", "Konya Plain", "Latium", "Paris Basin", 
                 "Middle Yellow River Valley", "Niger Inland Delta", "Orkhon Valley",
                 "Susiana", "Upper Egypt")
dat.p <- dat.MG.700[dat.MG.700$NGA %in% NGAs.unconq,]

ggplot(dat.p, aes(Time, Mean, color=MoralisingGods)) + 
  stat_summary(fun.data=mean_se, geom="pointrange", shape=18) + 
  scale_color_manual(values = alpha(c("coral1","aquamarine3"), .8), labels = c("NA", "Present"),
                     name = "Moralizing Gods") + 
  scale_x_continuous(limits = c(0,1400),
                     breaks = c(0,200,400,600,700,800,1000,1200,1400),
                     labels = c("-700","-500","-300","-100","0","100","300","500","700")) +
  geom_vline(xintercept = 650,
             color = "grey", size=6, alpha = 0.3) + 
  annotate("text", label = "MGs emerge", x = 650, y = 0.3, size = 4, colour = "black",
           angle = 90, alpha = 0.7) +
  labs(y="Social Complexity", x="Time (years before/after MG)") + 
  ylim(c(0.1,1)) +
  ggtitle("") +
  theme_bw() + 
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = rel(1.5)),        
    axis.line = element_line(colour = "black"),
    legend.position = "",
    legend.justification = c("right", "top"),
    legend.key.size = unit(0.8, "cm"),
    legend.title = element_text(size = rel(1.5)),
    legend.text = element_text(size = rel(1.2)),
    axis.title = element_text(size = rel(1.5)),
    axis.text.y= element_text(size = rel(1.5)),
    axis.text.x= element_text(size = rel(1.5)),
    plot.margin=unit(c(1,1,1,1),"cm")) 

##?OUR_COMMENT:: Save the plot if needed
ggsave("./hierarchical_models_output/PrePost_MG_unconq.png",width = 5.5,height = 3.7, dpi = 300)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
### 2.3. Concquered SC by MGs +/- 700 ####

NGAs.conq <- c("Deccan","Kachi Plain", 
              "Sogdiana")
dat.p <- dat.MG.700[dat.MG.700$NGA %in% NGAs.conq,]

ggplot(dat.p, aes(Time, Mean, color=MoralisingGods)) + 
  stat_summary(fun.data=mean_se, geom="pointrange", shape=17) + 
  scale_color_manual(values = alpha(c("coral1","aquamarine3"), .8), labels = c("NA", "Present"),
                     name = "Moralizing Gods") + 
  scale_x_continuous(limits = c(0,1400),
                     breaks = c(0,200,400,600,700,800,1000,1200,1400),
                     labels = c("-700","-500","-300","-100","0","100","300","500","700")) +
  geom_vline(xintercept = 650,
             color = "grey", size=6, alpha = 0.3) + 
  annotate("text", label = "MGs via conquest", x = 650, y = 0.38, size = 4, colour = "black",
           angle = 90, alpha = 0.3) +
  labs(y="Social Complexity", x="Time (years before/after MG)") + 
  ylim(c(0.1,1)) +
  ggtitle("") +
  theme_bw() + 
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = rel(1.5)),        
    axis.line = element_line(colour = "black"),
    legend.position = "",
    legend.justification = c("right", "top"),
    legend.key.size = unit(0.8, "cm"),
    legend.title = element_text(size = rel(1.5)),
    legend.text = element_text(size = rel(1.2)),
    axis.title = element_text(size = rel(1.5)),
    axis.text.y= element_text(size = rel(1.5)),
    axis.text.x= element_text(size = rel(1.5)),
    plot.margin=unit(c(1,1,1,1),"cm")) 

##?OUR_COMMENT:: Save the plot if needed
ggsave("./hierarchical_models_output/PrePost_MG_conq.png",width = 5.5,height = 3.7, dpi = 300)


##?OUR_COMMENT:: Calculate SC increases for conquered NGAs
dat.MG.700$Mean[dat.MG.700$NGA=="Deccan"]

mean(dat.MG.700$Mean[dat.MG.700$Time==700],na.rm = T)/
  mean(dat.MG.700$Mean[dat.MG.700$Time==600],na.rm = T)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## 3. LMMs of rate of SC change ####

### 3.1. Calculate SC rate ####

#Full time windows (up to 10,000 years before and after moralizing gods)
NGAs <- c("Big Island Hawaii", "Cambodian Basin", "Central Java", "Chuuk Islands",
          "Deccan",    "Ghanaian Coast", "Iceland", "Kachi Plain", "Kansai", "Kapuasi Basin",
          "Konya Plain", "Latium", "Middle Yellow River Valley", "Niger Inland Delta",
          "Orkhon Valley",  "Paris Basin", "Sogdiana", "Susiana", "Upper Egypt",
          "Yemeni Coastal Plain")

out <- matrix(NA, nrow=0, ncol=7)
for(i in 1:length(NGAs)){
  dt <- dat[dat$NGA == NGAs[i],]
  MG<-subset(dt,MoralisingGods=="1") #Replace "MoralisingGods" with "DoctrinalMode" or
                                      #"Writing" to do these analyses
  #library(dplyr)
  #MG<-as.data.frame(MG %>% group_by(PolID) %>% sample_n(size = 1))
  #randomly samples so there is only one century per polity
  MGAppear<-subset(MG, Time==min(Time))
  for(j in 1: 100){
    Earliest<-subset(dt, Time==MGAppear$Time-j*100) 
    Latest<-subset(dt, Time==MGAppear$Time+j*100)
    rates <- cbind(MGAppear$NGA,
                   ifelse(class(Earliest$Time)=="NULL","NA",
                  (MGAppear$Mean-Earliest$Mean)/(MGAppear$Time-Earliest$Time)),
                  ifelse(class(Latest$Time)=="NULL","NA",
                  ((Latest$Mean-MGAppear$Mean)/(Latest$Time-MGAppear$Time))),
                  (MGAppear$End-MGAppear$Start),j*100,MGAppear$World.Region,MGAppear$Family)
    out <- rbind(out,rates)
  }
  out <- rbind(out,rates)
}
colnames(out)<-c("NGA","PreRate","PostRate","MGUncertainty","TimeWindow","Region", "Lang")

#mean(out$MGUncertainty) #Use this while replacing "MoralisingGods" as above to get
#uncertainty values for time-series

write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE) #Exporting/importing to force it to
          # read as numeric (there is probably a more elegant way to do this)

out<-read.table("./temp/FullRates.csv", sep=",", header=TRUE)

out$Difference<-out[,3]-out[,2]

for(i in 2:length(out[,5])){
  out[i,8]<-out[i,5]-out[i-1,5]
}
out <-subset(out, out[,8]!=0) #getting rid of bug when the final row repeats in each NGA

write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE)

out <-subset(out, out[,5]<2050) #Change this to modify time-window restriction from 700 years
                                #pre/post moralizing gods (<750) or # out to use full time-window


##?OUR_COMMENT:: Create stacked data set for figures
out.s <- rbind(out,out)
out.s$Rate <- out.s$PostRate
out.s$prepost <- 1

out.s$Rate[1:length(out$PreRate)] <- out$PreRate
out.s$prepost[1:length(out$PreRate)] <- 0
out.s$prepost <- factor(out.s$prepost)

out.s$NGA <- as.factor(out.s$NGA)
out.s$Region <- as.factor(out.s$Region)
out.s$Lang <- as.factor(out.s$Lang)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
### 3.2. Plot nesting within NGAs ####
##?OUR_COMMENT:: Now explore data nesting within NGAs, which was suspected during
  # the t-test analysis.

NGAs <- c("Deccan", "Kachi Plain", "Kansai", "Konya Plain", "Latium",
          "Middle Yellow River Valley", "Niger Inland Delta", "Orkhon Valley",
          "Paris Basin", "Sogdiana", "Susiana", "Upper Egypt")

out.s <- out.s[out.s$NGA %in% NGAs,]

png("./hierarchical_models_output/NGA_nesting.png",width = 8,height = 5.5,units = 'in', res = 300)
pirateplot(formula = Rate ~ NGA,
           data = out.s,
           theme = 0,
           main = "",
           pal = "info2", # southpark color palette
           bean.b.o = 1, # Bean fill
           bean.f.o = .1, # Bean fill
           point.o = .9, # Points
           inf.f.o = .0, # Inference fill
           inf.b.o = 0, # Inference border
           avg.line.o = 0.7, # AverMMAT line
           bar.f.o = .0, # Bar
           inf.f.col = "white", # Inf fill col
           inf.b.col = "black", # Inf border col
           avg.line.col = "black", # avg line col
           point.pch = 21,
           point.bg = "white",
           point.cex = 1,
           point.lwd = 1,
           bean.lwd = 2,
           ylim = c(-0.002,0.003),
           ylab = "Rate of SC change"
           )
abline(0.0002511869,0.00, )
dev.off()


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
### 3.3. Linear models of the rate of change ####
##?OUR_COMMENT::
  # Starting with a straightforward linear model, we should first account for the nesting
  # of time-points within NGAs, then the nesting of time-points within world regions, and finally
  # let the pre-/post-MG effect vary by each NGA.

  # Note that we also wanted to explore whether nesting within language families will have
  # an effect, but the model gives a singular fit probably because region == language family
  # for the 12 NGAs we have in the analysis. We could have also tried the nesting of NGAs within
  # world.regions, but the estimated interaction on varying intercepts would be probably very
  # low and possibly cause convergence issues.


summary(lm1 <- lmer(Rate ~ prepost + (1|NGA),data = out.s))

summary(lm2 <- lmer(Rate ~ prepost + (1|NGA) + (1|Region),data = out.s))

summary(lm3 <- lmer(Rate ~ prepost + (prepost|NGA) + (1|Region),data = out.s))


##?OUR_COMMENT:: The Pre-/Post effects reported in Whitehouse et al. held up after accounting for
  # the interdependencies among time-points. 
  # The nesting within NGAs and world regions did not explain much of the variability,
  # but in combination with the varying effects of pre-post showed a minor influence.

##?OUR_COMMENT:: Let's check the model fit.
print("For explanation of this warning, see the code")
simulationOutput = simulateResiduals(lm3)
#png("./hierarchical_models_output/qq1.png",width = 8,height = 5, units = 'in', res = 300)
plot(simulationOutput)
#dev.off()
testUniformity(simulationOutput = simulationOutput)
testDispersion(simulationOutput = simulationOutput)


##?OUR_COMMENT:: Normal distribution doesn't look like an appropriate assumption here.
  # This is probably caused by the fact that we are looking already at the rate of SC change
  # rather than modeling SC at each time point with a proper hierarchical model.
  # Let's try to model just SC, not the computed rate of change with the factorial MG variable.
  # However, to do that, we have to treat MG == NA as MG == 0,
  # otherwise there is nothing to model.


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## 4. GLMMs of SC ####

##?OUR_COMMENT:: First, let's look at the SC distribution.
dat.p <- na.omit(subset(dat, select = c(MoralisingGods, Mean)))
dat.p$MoralisingGods <- as.factor(dat.p$MoralisingGods)

ggplot() + 
  geom_density(data = dat.p, aes(x = Mean, fill = MoralisingGods, colour = MoralisingGods,
                                 linetype=MoralisingGods),size = 0.5) +
  scale_color_manual(values = alpha(c("lightskyblue4","lightskyblue4"), 1), guide = FALSE) + 
  scale_fill_manual(breaks = c(0,1),values = alpha(c("coral2","turquoise3"),0.6),
                     labels = c("Pre-MG", "Post-MG"), name = "") +
  scale_linetype_manual(values = c("solid","dashed"), guide = FALSE) +
  xlim(c(0,1)) +
  labs(x="Social Complexity", y="Density") +
  theme_bw() + 
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = rel(1.5)),        
    axis.line = element_line(colour = "black"),
    legend.position = c(0.25,0.95),
    legend.justification = c("right", "top"),
    legend.key.size = unit(0.8, "cm"),
    legend.text = element_text(size = rel(1.5)),
    axis.title = element_text(size = rel(1.5)),
    axis.text.y= element_text(size = rel(1.5)),
    axis.text.x= element_text(size = rel(1.5)),
    plot.margin=unit(c(1,1,1,1),"cm")) 


##?OUR_COMMENT:: Even after accounting for the Pre-/Post-MG difference,
  # the residuals will most likely not be normally distributed.
  # Given how the data were produced, the most suited distributional assumption is the beta model.
  # (data scaled on [0,1] interval)


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
### 4.1 +/- 2000 years analysis ####

#### 4.1.1. Compare beta with linear ####

##?OUR_COMMENT:: First, let's compare the gaussian and beta models
dat.MGn <- na.omit(subset(dat.MG.2000, select = c(Mean, Time, NGA)))

summary(lm1 <- glmmadmb(Mean ~ Time + (1|NGA), data = dat.MGn, family = 'gaussian'))

summary(glm.b1 <- glmmadmb(Mean ~ Time + (1|NGA), data = dat.MGn, family = 'beta'))

a <- AICtab(lm1, glm.b1, sort = F)
t <- matrix(c(a$dAIC[1],a$dAIC[2],a$df[1],a$df[2]),ncol=2,byrow=F)
colnames(t) <- c("dAIC","df")
rownames(t) <- c("lm1","glm.b1")
t

#______________________________________________________________________________________________
#### 4.1.2. Beta models ####

##?OUR_COMMENT:: glmmadmb doesn't handle NAs well, let's select the variables we will
  # work with and delete NAs manually
dat.MG1 <- na.omit(subset(dat.MG.2000, select = c(Mean, Time, MoralisingGods, NGA, World.Region)))

##?OUR_COMMENT:: The coeficients for Time will be very small, let's transform them to look
  # at average change per millenium rather than century
dat.MG1$Time <- dat.MG1$Time/1000

summary(glm.b1 <- glmmadmb(Mean ~ Time*MoralisingGods + (1|NGA), data = dat.MG1, family = 'beta'))

summary(glm.b2 <- glmmadmb(Mean ~ Time*MoralisingGods + (1|NGA) + (1|World.Region),
                           data = dat.MG1, family = 'beta'))

summary(glm.b3 <- glmmadmb(Mean ~ Time*MoralisingGods + (Time|NGA) + (1|World.Region),
                           data = dat.MG1, family = 'beta'))

##?OUR_COMMENT:: Compare models' AIC
a <- AICtab(glm.b3,glm.b2,glm.b1, sort = F)
t <- matrix(c(a$dAIC[1],a$dAIC[2],a$dAIC[3],a$df[1],a$df[2],a$df[3]),ncol=2,byrow=F)
colnames(t) <- c("dAIC","df")
rownames(t) <- c("glm.b3","glm.b2","glm.b1")
t

##?OUR_COMMENT:: Transform coefficients and compute 95%CI
cf <- summary(glm.b3)$coefficients[,1]
se <- sqrt(diag(vcov(glm.b3)))
ci <- (cbind(est = (plogis(cf)-0.5), LL = (plogis(cf - 1.96 * se)-0.5),
             UL = (plogis(cf + 1.96 * se)-0.5)))
ci[1,] <- (ci[1,]+0.5)

print("Beta model for MG = 0 with +/-2000 time-span")
print(ci)
ci[2,1]+ci[4,1] ##?OUR_COMMENT:: Post-MG slope

##?OUR_COMMENT:: Assess goodness of fit. 
  # We have to use the glmmTMB package because glmmADMB is not compatible with DHARMa.

summary(glm.b3 <- glmmTMB(Mean ~ Time*MoralisingGods + (Time|NGA) + (1|World.Region),
                          data = dat.MG.2000, family = 'beta'))

print("For explanation of this warning, see the code")
simulationOutput = simulateResiduals(glm.b3)
#png("./hierarchical_models_output/qq3.png",width = 8,height = 5, units = 'in', res = 300)
plot(simulationOutput)
#dev.off()
testUniformity(simulationOutput = simulationOutput)
testDispersion(simulationOutput = simulationOutput)


##?OUR_COMMENT:: The results of the multilevel growth curve model revealed that by using a more
  # appropriate statistical approach, the original t-test result presented by Whitehouse et al.
  # does not hold anymore. While the SC growth after MG appearance is more negative compared
  # to the SC growth before MG, this difference is negligible.

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
### 4.2 +/- 700 years analysis ####

#### 4.2.1. Compare beta with linear ####

##?OUR_COMMENT:: Here, we use AIC to compare gaussian and beta assumptions.

summary(lm1 <- glmmadmb(Mean ~ Time + (1|NGA), data = dat.MG.700, family = 'gaussian'))

summary(glm.b1 <- glmmadmb(Mean ~ Time + (1|NGA), data = dat.MG.700, family = 'beta'))

a <- AICtab(lm1, glm.b1, sort = F)
t <- matrix(c(a$dAIC[1],a$dAIC[2],a$df[1],a$df[2]),ncol=2,byrow=F)
colnames(t) <- c("dAIC","df")
rownames(t) <- c("lm1","glm.b1")
t

##?OUR_COMMENT::The beta model seems to fit the data better.

#______________________________________________________________________________________________
#### 4.2.2. Beta models ####

##?OUR_COMMENT:: This is the main multi-level growth curve analysis assuming a beta distribution
# of residuals. We fit three models, analogically to the hierarchical models presented in
# section 3 (see that section for the rationale of our modeling building procedure). 

# Note that it would be possible to also add a quadratic effect of time or to compare the
# shape of the growth curve (linear vs. quadratic) between conditions.
# However, this would venture too far beyond the original assumptions of
# Whitehouse et al.'s paper.


summary(glm.b1 <- glmmadmb(Mean ~ Time*MoralisingGods + (1|NGA), data = dat.MG.700, family = 'beta'))

summary(glm.b2 <- glmmadmb(Mean ~ Time*MoralisingGods + (1|NGA) + (1|World.Region),
                           data = dat.MG.700, family = 'beta'))

summary(glm.b3 <- glmmadmb(Mean ~ Time*MoralisingGods + (Time|NGA) + (1|World.Region),
                           data = dat.MG.700, family = 'beta'))

a <- AICtab(glm.b3,glm.b2,glm.b1, sort = F)
t <- matrix(c(a$dAIC[1],a$dAIC[2],a$dAIC[3],a$df[1],a$df[2],a$df[3]),ncol=2,byrow=F)
colnames(t) <- c("dAIC","df")
rownames(t) <- c("glm.b3","glm.b2","glm.b1")
t

##?OUR_COMMENT:: The third model fits the data best. 

##?OUR_COMMENT:: Transform coefficients and compute 95%CI
cf <- summary(glm.b3)$coefficients[,1]
se <- sqrt(diag(vcov(glm.b3)))
ci <- (cbind(est = (plogis(cf)-0.5), LL = (plogis(cf - 1.96 * se)-0.5),
             UL = (plogis(cf + 1.96 * se)-0.5)))
ci[1,] <- (ci[1,]+0.5)

print("Beta model for MG = 0 with +/-700 time-span")
print(ci)

##?OUR_COMMENT:: Let's explore goodness-of-fit measures.

##?OUR_COMMENT:: We have to use the glmmTMB package because glmmADMV is not compatible with
# DHARMa

summary(glm.b3 <- glmmTMB(Mean ~ Time*MoralisingGods + (Time|NGA) + (1|World.Region),
                          data = dat.MG.700, family = 'beta'))

print("For explanation of this warning, see the code")
simulationOutput = simulateResiduals(glm.b3)
#png("./temp/qq2.png",width = 7,height = 5, units = 'in', res = 300)
plot(simulationOutput)
#dev.off()
testUniformity(simulationOutput = simulationOutput)
testDispersion(simulationOutput = simulationOutput)


##?OUR_COMMENT:: The glm.b3 model appears to be a pretty good fit based on the DHARMa estimation,
# although we are still missing something and the data seem to be underdispersed. 
# Note that DHARMa can't yet handle random effects in GLMM in the residual vs. predicted plot).


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
### 5. Beta models with forward bias correction ####

#### 5.1. Shift MGs 100 years back ####

##?OUR_COMMENT:: In this final set of analyses, we combine our modeling approach with
  #the correction for forward bias.

##?OUR_COMMENT:: Here we shift MG 100 years back; i.e., on the time line from 0 to 3900 years
  # where MGs appear in the middle (2000), we now say appear in year 1900
dat.FB1 <- dat.MG.2000

dat.FB1$MoralisingGods[dat.FB1$Time>1800] <- 1 

dat.FB1$NGA <- as.factor(dat.FB1$NGA)

dat.FB1 <- na.omit(subset(dat.FB1, select = c(Mean, Time, MoralisingGods, NGA, World.Region)))

##?OUR_COMMENT:: The coeficients for Time will be very small,
  # let's transform them to look at average change per millenium rather than century
dat.FB1$Time <- dat.FB1$Time/1000

summary(glm.b1 <- glmmadmb(Mean ~ Time*MoralisingGods + (1|NGA), data = dat.FB1, family = 'beta'))

summary(glm.b2 <- glmmadmb(Mean ~ Time*MoralisingGods + (1|NGA) + (1|World.Region),
                           data = dat.FB1, family = 'beta'))

summary(glm.b3 <- glmmadmb(Mean ~ Time*MoralisingGods + (Time|NGA) + (1|World.Region),
                           data = dat.FB1, family = 'beta'))

a <- AICtab(glm.b3,glm.b2,glm.b1, sort = F)
t <- matrix(c(a$dAIC[1],a$dAIC[2],a$dAIC[3],a$df[1],a$df[2],a$df[3]),ncol=2,byrow=F)
colnames(t) <- c("dAIC","df")
rownames(t) <- c("glm.b3","glm.b2","glm.b1")
t

##?OUR_COMMENT:: Transform coefficients and compute 95%CI
cf <- summary(glm.b3)$coefficients[,1]
se <- sqrt(diag(vcov(glm.b3)))
ci <- (cbind(est = (plogis(cf)-0.5), LL = (plogis(cf - 1.96 * se)-0.5),
             UL = (plogis(cf + 1.96 * se)-0.5)))
ci[1,] <- (ci[1,]+0.5)

print("Beta model for MG - 100 with +/-2000 time-span")
print(ci)
##?OUR_COMMENT:: Again, we see an opposite trend than reported by Whitehouse et al.

#______________________________________________________________________________________________
#### 5.2. Shift MGs 300 years back ####

##?OUR_COMMENT:: Now shift MGs 300 years back.

##?OUR_COMMENT:: Here we shift MG 300 years back; i.e., on the time line from 0 to 4000 years
  # where MGs appear in the middle (2000), we now say appear in year 1700
dat.FB3 <- dat.FB1

dat.FB3$MoralisingGods[dat.FB3$Time>1600] <- 1 

dat.FB3$NGA <- as.factor(dat.FB3$NGA)

dat.FB3 <- na.omit(subset(dat.FB3, select = c(Mean, Time, MoralisingGods, NGA, World.Region)))

##?OUR_COMMENT:: The coeficients for Time will be very small,
  #let's transform them to look at average change per millenium rather than century
dat.FB3$Time <- dat.FB3$Time/1000

summary(glm.b1 <- glmmadmb(Mean ~ Time*MoralisingGods + (1|NGA), data = dat.FB3, family = 'beta'))

summary(glm.b2 <- glmmadmb(Mean ~ Time*MoralisingGods + (1|NGA) + (1|World.Region),
                           data = dat.FB3, family = 'beta'))

summary(glm.b3 <- glmmadmb(Mean ~ Time*MoralisingGods + (Time|NGA) + (1|World.Region),
                           data = dat.FB3, family = 'beta'))


a <- AICtab(glm.b3,glm.b2,glm.b1, sort = F)
t <- matrix(c(a$dAIC[1],a$dAIC[2],a$dAIC[3],a$df[1],a$df[2],a$df[3]),ncol=2,byrow=F)
colnames(t) <- c("dAIC","df")
rownames(t) <- c("glm.b3","glm.b2","glm.b1")
t

##?OUR_COMMENT:: Transform coefficients and compute 95%CI
cf <- summary(glm.b3)$coefficients[,1]
se <- sqrt(diag(vcov(glm.b3)))
ci <- (cbind(est = (plogis(cf)-0.5), LL = (plogis(cf - 1.96 * se)-0.5),
             UL = (plogis(cf + 1.96 * se)-0.5)))
ci[1,] <- (ci[1,]+0.5)
print("Beta model for MG - 300 with +/-2000 time-span")
print(ci)


dir_init("./output")

files <- list.files("./hierarchical_models_output", full.names = TRUE)
file.copy(files, "./output")