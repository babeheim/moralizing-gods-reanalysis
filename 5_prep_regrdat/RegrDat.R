#Constructing RegrDat file for regression analyses
#setwd("/Users/pesavage/Documents/Research/Oxford Seshat/Data/SCBigGodsOct2017")

rm(list = ls())
source("../project_support.R")

# load("./HS.Rdata")

DistMatrix <- read.csv("./input/DistMatrix.csv", stringsAsFactors = FALSE)
DistMatrix <- as.matrix(DistMatrix)
colnames(DistMatrix) <- gsub("\\.", " ", colnames(DistMatrix))
rownames(DistMatrix) <- colnames(DistMatrix)

# this Rdata file contains:
# "coords"     "DistMatrix" "NGAs"       "polities"   "RegrDat"

# but polities is overwritten before being called
# RegrDat is overwritten before being called
# the DistMatrix is used!
# NGAs and Coords are never used
# so basically we just need the DistMatrix, not HS.Rdata!!

dir_init("./temp")

# polities here taken from impute_data...not! where did they change upstream??

polities <- read.csv('./input/polities.csv', header=TRUE, stringsAsFactors = FALSE)
# There is an issue with compatibility between my polities.csv file and
# the one Peter was  using for HS. Resolving this as follows based on debugging
# the problem with non-matching PolID, identified by using the following code at
# the end of the original code rerun phylogeny section:
#setdiff(levels(polities$PolID),levels(PolID))
#setdiff(levels(PolID),levels(polities$PolID))
#[1] "GapDec"  "GapKP1"  "IrQajar" "TrOttm5" "YeOttoL" "YeTahir" "YeZiyad"
polities <- polities[polities$PolID != "GapDec",]
polities <- polities[polities$PolID != "GapKP1",]
polities <- polities[polities$PolID != "IrQajar",]
polities <- polities[polities$PolID != "TrOttm5",]
polities <- polities[polities$PolID != "YeOttoL",] 
write.csv(polities, file="./temp/polities.csv",  row.names=FALSE) 
polities <- read.csv('./temp/polities.csv', header=TRUE, stringsAsFactors = FALSE)
##
data <- read.table("./input/PC1_traj_merged.csv", sep=",", header=TRUE, stringsAsFactors = FALSE) # from 3_run_pca
data[is.na(data)] <- 0 #This treats NA values as 0. Should checek later to see how much this affects results
data$MG<-data$MoralisingGods

refcols <- c("NGA", "PolID","Time","MG")
data <- data[, c(refcols, setdiff(names(data), refcols))]
names(data)

### Add lags
RegrDat<-data #change label for compatibility with earlier code
dat <- cbind(RegrDat,matrix(NA,nrow(RegrDat),2))
nm <- colnames(dat)
nm[c((length(nm)-1),length(nm))] <- c("Lag1","Lag2")
colnames(dat) <- nm
for(i in 2:nrow(dat)){
  if(dat$Time[i] == dat$Time[i-1]+100){dat$Lag1[i] <- dat$MG[i-1]}
}
for(i in 3:nrow(dat)){
  if(dat$Time[i] == dat$Time[i-2]+200){dat$Lag2[i] <- dat$MG[i-2]}
}
RegrDat <- dat

#### Calculate Space using the estimated d = 1000 km (suggested by Peter Turchin as default because this is approximate average distance between NGAs)
dpar <- 1000
Space <- RegrDat[,1:4]
Space[,4] <- NA
colnames(Space) <- c("NGA","PolID","Time","Space")

colMat <- colnames(DistMatrix)
rowMat <- rownames(DistMatrix)
for(i in 1:nrow(RegrDat)){
  t1 <- RegrDat$Time[i] - 100
  dat <- RegrDat[RegrDat$Time==t1,1:4] 
  if(nrow(dat) > 1){
    delta <- vector(length=nrow(dat))
    for(j in 1:nrow(dat)){
      dt <- DistMatrix[colMat==dat$NGA[j],]
      delta[j] <- dt[rowMat==RegrDat$NGA[i]]
    }
    s <- exp(-delta/dpar)*dat$MG
    s <- s[delta != 0]  ### Exclude i=j
    Space$Space[i] <- mean(s)
  }
}
Space$Space[is.na(Space$Space)] <- mean(Space$Space, na.rm = TRUE) # Substitute NAs with the mean
Space$Space <- Space$Space/max(Space$Space)   ##### Scale to max = 1
RegrDat <- cbind(RegrDat,Space$Space)
nm <- colnames(RegrDat)
nm[length(nm)] <- "Space"
colnames(RegrDat) <- nm

#### Calculate Language = matrix of linguistic distances
Phylogeny <- RegrDat[,1:4]
Phylogeny[,4] <- NA
colnames(Phylogeny) <- c("NGA","PolID","Time","Phylogeny")

for(i in 1:nrow(RegrDat)){
  t1 <- RegrDat$Time[i] - 100
  dat <- RegrDat[RegrDat$Time==t1,1:4]
  dat <- dat[dat$NGA != RegrDat$NGA[i],]   ### Exclude i = j
  PolID <- RegrDat$PolID[i]
  PolLang <- polities[polities$PolID==PolID,9:11]
  if(nrow(dat) > 1){
    weight <- vector(length=nrow(dat)) * 0
    for(j in 1:nrow(dat)){
      dt <- dat[j,]
      PolLang2 <- polities[polities$PolID==dt$PolID,9:11]
      if(PolLang[1,3]==PolLang2[1,3]){weight[j] <- 0.25}
      if(PolLang[1,2]==PolLang2[1,2]){weight[j] <- 0.5}
      if(PolLang[1,1]==PolLang2[1,1]){weight[j] <- 1}
    }
    s <- weight*dat$MG
    Phylogeny$Phylogeny[i] <- mean(s)
  }
}
Phylogeny$Phylogeny[is.na(Phylogeny$Phylogeny)] <- 0 # Substitute NAs with 0
RegrDat <- cbind(RegrDat,Phylogeny$Phylogeny)
nm <- colnames(RegrDat)
nm[length(nm)] <- "Phylogeny"
colnames(RegrDat) <- nm

dir_init("./output")

write.csv(RegrDat, "./output/RegrDat.csv", row.names = FALSE)
