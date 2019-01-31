#The following code was used for analyses in Whitehouse, Francois, Savage, et al., "Complex societies and doctrinal rituals precede moralizing gods throughout world history", with an "exportdat.csv" file scraped from the Seshat database on 19 Jan 2018.
#This software was developed by Peter Turchin and Patrick Savage (Copyright 29 Jan 2018). For the general methodology underlying these analyses please refer to: 
#Turchin, P. et al. Quantitative historical analysis uncovers a single dimension of complexity that structures global variation in human social organization. Proc. Natl. Acad. Sci. U. S. A. 115, E144-E151 (2018).
#Turchin, P. Fitting dynamical regression models to Seshat data. Cliodynamics 9, (2018).  

######
#CC By-NC SA License

#Copyright (c) 2018 Peter Turchin and Patrick E. Savage

#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software under the conditions of Creative Commons Attribution Non-Commercial (CC By-NC SA) licensing (https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode), subject to the following conditions:

#Please include the following text in any publication using this software:

#This research employed data from the Seshat Databank (seshatdatabank.info) under Creative Commons Attribution Non-Commercial (CC By-NC SA) licensing (https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode).

#and cite:

#1) Whitehouse, Francois, Savage, et al., "Complex societies and doctrinal rituals precede moralizing gods throughout world history" (2019) Nature.

#2) Turchin, P. et al. Quantitative historical analysis uncovers a single dimension of complexity that structures global variation in human social organization. Proc. Natl. Acad. Sci. U. S. A. 115, E144-E151 (2018).

#3) Turchin, P. Fitting dynamical regression models to Seshat data. Cliodynamics 9, (2018).

#4) Turchin P. et al. 2015. Seshat: The Global History Databank. Cliodynamics 6(1):77–107. 

#The views and conclusions contained in this document are those of the authors and should not be interpreted as representing the official positions, either expressed or implied, of the Seshat Databank, its collaborative scholarly community, or the Evolution Institute.

#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
######

#The following packages must be installed/loaded
library(maps)
library(plotrix)
library(plyr) #(But maybe shouldn't load this when loading dplyr for confirmatory analyses) 

# Make sure to use a full polities.csv file (if you re-run using the output the multiple imputation loop will run into bugs when it comes across NGAs with only one polity)

### NB: Using shortened variables.csv file that only uses variables of interest to save computational time. If using full file, may need to renumber variables in Aggr.MI

# Don't rerun this section unless it's a new scrape

setwd("/Users/pesavage/Documents/Research/Oxford Seshat/Data/MoralizingGods")
Section1 <- "Social Complexity variables"
Section2 <- "Ritual variables"
Section3 <- "Religion and Normative Ideology"
source("precheck.R")

############################################ Start here if precheck.R is done

Vars <- as.matrix(read.csv('variables.csv', header=TRUE))
Vars <- Vars[Vars[,6]==Section1 | Vars[,6]==Section2 | Vars[,6]==Section3,] # Reduce the variables list to the Section set above

Vars[,1] <- paste(Vars[,2],Vars[,1]) #Creating unique variable/section combinations

polities <- read.csv('polities.csv', header=TRUE)
polities <- polities[polities$NGA != "Crete",] #### Remove new NGAs
polities <- polities[polities$NGA != "Galilee",]
polities <- polities[polities$NGA != "Middle Ganga",]
polities <- polities[polities$PolID != "IdCJBun",] #### Remove low-coverage polities causing bugs in C. Java)
polities <- polities[polities$PolID != "IdKalin",] 
polities <- polities[polities$PolID != "EsHabsb",] #### Remove post-colonial polities (also see below for removing post-colonial polities from NGAs with only 2 polities)
polities <- polities[polities$PolID != "USIllinL",]
polities <- polities[polities$PolID != "GbEmpir",]
polities <- polities[polities$PolID != "InEInCo",]
polities <- polities[polities$PolID != "InBritP",]
polities <- polities[polities$PolID != "RuYakuL",]
polities <- polities[polities$PolID != "UsIroqL",]

write.csv(polities, file="polities.csv",  row.names=FALSE)
polities <- read.csv('polities.csv', header=TRUE)
NGAs <- levels(polities$NGA)

nrep <- 20
ImpDatRepl <- matrix(NA, nrow=0, ncol=0) 
for(irep in 1:nrep){
  print(irep)
  source("ConstrMI.R")
  source("AggrMI.R")
  source("ImputeMI.R")
  ones <- matrix(data=1,nrow=length(AggrDat[,1]),ncol=1)
  colnames(ones) <- "irep"
  ImpDat <- cbind(AggrDat[,1:4],ImpDat,(ones*irep),AggrDat[,14:32])
  ImpDatRepl <- rbind(ImpDatRepl,ImpDat)
}

####### Remove polity-dates that didn't yield 20 repl #and post-colonial polities that couldn't be removed from multiple imputation due to bugs with only 1 polity/NGA

polities <- read.csv('polities.csv', header=TRUE)
polities <- polities[polities$PolID != "InGaroL",] #removing here because it caused bugs earlier
write.csv(polities, file="polities.csv",  row.names=FALSE) 
polities <- polities[polities$PolID != "CnHChin",] #removing here because it caused bugs earlier
write.csv(polities, file="polities.csv",  row.names=FALSE) 
polities <- polities[polities$PolID != "PgOrokL",] #removing here because it caused bugs earlier
write.csv(polities, file="polities.csv",  row.names=FALSE) 

ImpDatRepl <- ImpDatRepl[ImpDatRepl$PolID != "InGaroL",] #removing here because it seemed to create bugs when you have only 1 polity in an NGA, so couldn't remove earlier
ImpDatRepl <- ImpDatRepl[ImpDatRepl$PolID != "CnHChin",] #removing here because it seemed to create bugs when you have only 1 polity in an NGA, so couldn't remove earlier
ImpDatRepl <- ImpDatRepl[ImpDatRepl$PolID != "PgOrokL",] #removing here because it seemed to create bugs when you have only 1 polity in an NGA, so couldn't remove earlier


dat_temp <- ImpDatRepl
for(i in 1:nrow(polities)){
  dat <- ImpDatRepl[as.character(ImpDatRepl[,2])==as.character(polities[i,2]),]
  if(nrow(dat)!=0){
    Time <- unique(dat$Time)
    for(j in 1:length(Time)){
      dt <- dat[dat$Time==Time[j],]
      if(nrow(dt) != nrep){
        print(nrow(dt))
        print(dt[1,1:3])
        dat_temp[as.character(dat_temp$PolID)==as.character(dat$PolID[1]) & dat_temp$Time==Time[j],14] <- -99999
      }
    }
  }
}
ImpDatRepl <- dat_temp[dat_temp$irep!=-99999,]



write.csv(ImpDatRepl, file="ImpDatRepl.csv",  row.names=FALSE)

#######################################################################
######### end of the new scrape section


# Run PCA on ImpDat
ImpDatRepl <- read.csv('ImpDatRepl.csv', header=TRUE)
ImpDat <- ImpDatRepl[ImpDatRepl$irep==1,]
nrep <- ImpDatRepl$irep[nrow(ImpDatRepl)]
PC1 <- matrix(NA,(length(ImpDat[,1])),0)
Rotations <- matrix(NA,0,9)
PropVar <- matrix(NA,0,9)
for(irep in 1:nrep){
  print(irep)
  ImpDat <- ImpDatRepl[ImpDatRepl$irep==irep,5:13] #Replace "5:13" with "5:8" and "9:13" for confirmatory analyses (changing "9" in "Rotations <- matrix(NA,0,9)" and "PropVar <- matrix(NA,0,9)" to "4" or "5" as appropriate)
  res <- prcomp(ImpDat, scale=TRUE)
  print(summary(res))
  PCA<-predict(res)
  
  PC1 <- cbind(PC1,PCA[,1])
  Rotations <- rbind(Rotations,res$rotation)
  PropVar <- rbind(PropVar,summary(res)$importance[2,])
}

write.csv(PC1, file="PC1.csv",  row.names=FALSE)
write.csv(Rotations, file="Rotations.csv",  row.names=TRUE)
write.csv(PropVar, file="PropVar.csv",  row.names=FALSE)

#rm(AggrDat, ImpDat, PCA, res)

##### gPC1.R
## Calculate PC1 trajectories +/- 2SD conf int [lower, mean, upper]
#setwd("C:/Users/Peter Turchin/Google Drive/2.Seshat/1.R/PCA-MI")
PC1 <- read.table("PC1.csv", sep=",", header=TRUE)

#PC1[,5] <- -(PC1[,5]) #NB: It is possible that the signs may be arbitrarily flipped in different PCA runs, in which case used this type of code to correct arbitrarily flipped signs


PC1 <- (PC1-min(PC1))/max(PC1-min(PC1)) # Scale PC1 to between 0-1
dat <- ImpDatRepl[ImpDatRepl$irep==1,]
dat[,1] <- as.character(dat[,1])
dat[,2] <- as.character(dat[,2])
dat$Mean <- apply(PC1,1,mean)
dat$Lower <- dat$Mean - 2*apply(PC1,1,sd)
dat$Upper <- dat$Mean + 2*apply(PC1,1,sd)

output <- read.table("MIoutput.csv", sep=",", header=TRUE)
dt <- output[,1:length(dat)]
dt[,1] <- as.character(dt[,1])
dt[,2] <- as.character(dt[,2])
dt[,4:length(dat)] <- NA
#names(dt) <- c("NGA","PolID","Time","Mean","Lower","Upper")
names(dt) <- names(dat)

for(i in 1:length(dat[,1])){
  j <- 1:length(dt[,1])
  j <- j[dt[,1]==dat[i,1] & dt[,2]==dat[i,2] & dt[,3]==dat[i,3]]
  dt[j,] <- dat[i,]
}
for(i in 1:length(dt[,1])){   if(is.na(dt[i,4])){dt[i,4:length(dat)] <- dt[i-1,4:length(dat)]} }
dat <- dt

write.csv(dat, file="PC1_traj.csv",  row.names=FALSE)



################################################

### Graphs PCA results for Figure 3 of the PCA article. (see Figure3.R for more code) 
###    and for the Suppl Materials: all regions
#setwd("C:/Users/Peter Turchin/Google Drive/2.Seshat/1.R/PCA-MI")
dat <- read.table("PC1_traj.csv", sep=",", header=TRUE)
NGAs <- c("Upper Egypt","Middle Yellow River Valley")

colors <- c("red","blue","darkgreen","brown","tan","black","orange","cyan","darkgrey")
colors <- c(colors,colors,colors)
shade_colors <- c("pink","lightblue","lightgreen","tan")
x_ends <- c(-5000,2000)
y_ends <- c(0,1)
plot(x_ends, y_ends, "n", xlim = x_ends, ylim=y_ends, xlab="year", ylab="PC1", main = "")
legend("topleft", NGAs, lty=1, lwd=3, bty="n", col=colors)

for(j in 1:length(NGAs)){
  NGA <- NGAs[j]
  gdat <- dat[dat[,1]==NGA,]
  polygon(x=c(gdat[,3],rev(gdat[,3])),y=c(gdat$Upper,rev(gdat$Lower)), col=shade_colors[j], border=NA)
  lines(gdat[,3],gdat$Mean, lwd=3, col=colors[j])
  #   lines(gdat[,3],gdat[,5], lwd=1, lty = 2, col=colors[j])
  #   lines(gdat[,3],gdat[,6], lwd=1, lty = 2, col=colors[j])
  #   text(x=x_ends[1], y=(y_ends[2] - .5*j), NGA, col=colors[j], pos=4)
}

################################################
#Merge PC1_traj.csv with other key info
dat <- read.table("PC1_traj.csv", sep=",", header=TRUE)
dat2 <- read.table("polities.csv", sep=",", header=TRUE)
dat <- merge(dat,dat2,by="PolID")
dat2 <- read.table("NGAcoords.csv", sep=",", header=TRUE)
dat <- merge(dat,dat2,by.x="NGA.y",by.y="NGA")
#library(plyr)
dat<-rename(dat, c("NGA.y" = "NGA"))
dat<-rename(dat, c("NGA.x" = "OriginalNGA"))
dat<-dat[order(dat$NGA, dat$Time),]
dat<-subset(dat,dat$Time>=dat$Start & dat$Time<=dat$End)
write.csv(dat, file="PC1_traj_merged.csv",  row.names=FALSE)

################################################
###Run pre-/post-moralizing gods analyses:
source("BigGodAnalysesEditedV2.R")

###Run logistic regression
source("RegrDat.R")

#   FULL MODEL
RD1 <- RegrDat[is.na(RegrDat$Lag1) == FALSE,]
RD2 <- RD1[is.na(RD1$Lag2) == FALSE,]
LogistRegrDat = as.data.frame(cbind(RD2$MG,RD2[,c(36,51:54)]))
source("LogistRegr.R")

######
##Creat Fig. 1 map using manually created map.csv file
source("Seshat_BG_map.r")
