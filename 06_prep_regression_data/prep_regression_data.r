
rm(list = ls())
source("../project_support.r")

dir_init("./temp")

# load("./HS.Rdata") <- only DistMatrix was needed from this

DistMatrix <- read.csv("./input/DistMatrix.csv", stringsAsFactors = FALSE)
DistMatrix <- as.matrix(DistMatrix)
colnames(DistMatrix) <- gsub("\\.", " ", colnames(DistMatrix))
rownames(DistMatrix) <- colnames(DistMatrix)

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

####

data <- read.table("./input/PC1_traj_merged.csv",
  sep=",", header=TRUE, stringsAsFactors = FALSE)

expect_equal(dim(data), c(864, 49))
expect_equal(length(unique(data$NGA)), 30)
expect_equal(sum(is.na(data$MoralisingGods)), 528)
expect_equal(sum(is.na(data$MoralisingHighGods)), 513)
expect_equal(sum(is.na(data$GeneralMoralisticPunishment)), 545)

data$MG <- data$MoralisingGods

# flags added to keep track of 0's imputed from NA's
data$MG_missing <- as.numeric(is.na(data$MG))
data$MHG_missing <- as.numeric(is.na(data$MHG))
data$BSP_missing <- as.numeric(is.na(data$BSP))

data[is.na(data)] <- 0 
# This treats NA values as 0. Should checek later to see how much this affects results

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

expect_equal(dim(RegrDat), c(864, 52 + 3)) # i added MG_missing, MHG_missing, BSP_missing
expect_equal(sum(is.na(RegrDat$Lag1)), c(32))
expect_equal(sum(is.na(RegrDat$Lag2)), c(63))

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

expect_true(abs(mean(RegrDat$Space) - 0.2) < 0.01)

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

expect_true(abs(mean(RegrDat$Phylogeny) - 0.032) < 0.0003)

dir_init("./output")

write.csv(RegrDat, "./output/RegrDat.csv", row.names = FALSE)

print("regression table RegrDat.csv created")
