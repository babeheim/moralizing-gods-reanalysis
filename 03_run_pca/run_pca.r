
rm(list = ls())
source("../project_support.r")

dir_init("./temp")

# Nature
# In order to create an overall measure of social complexity,we took a previously published approach based on principal component analysis (PCA) and applied it to the latest available data from Seshat. This method aggregates the 51 social complexity variables (Extended Data Table 5) into nine complexity characteristics and then analyses them using PCA.
# PCA is a commonly used tool for dimension reduction—in this case we have nine different aggregated variables that we want to reduce to a single variable that best captures social complexity. However, we obtain the same conclusions even without using PCA regardless of which of the nine complexity characteristics we choose as a proxy for social complexity.

# these different complexity characteristics turn out to be highly correlated and all load heavily onto a single principal component that captures 76% of the variance in the individual complexity characteristic variables.
# PNAS
# All nine CCs were highly and significantly correlated with each other. Correlation coefficients varied between 0.49 (government and writing) and 0.88 (polity population and polity territory). Only a single principal component, PC1, has an eigenvalue greater than 1. It explains 77.20.4 percent of variance. The proportion of variance explained by other principal components drops rapidly towards zero (e.g. PC2 explains only 6.00.4 percent). Furthermore, when we examine the “loadings” of the nine variables on PC1 (correlations between raw variables and PCs), we observe that all variables contribute about equally to PC1.
# Loadings on PC2 seem to capture a slight residual but negative relationship between “social scale” variables (capital and polity population, hierarchical levels, and polity territory) and informational/economic complexity (writing, texts, and money). This could reflect cases where these features have diffused from large-scale societies where they were originally developed to smaller societies, or cases where large-societies have reduced in size but still retained these features. However, it should be emphasized that PC2 is not well-supported so we should be careful not to over-interpret this result. Overall, these results support the idea that different aspects of social organization have co-evolved in predictable ways, and that social complexity is a concept that can be well-represented by a measure such as PC1.

# Run PCA on ImpDat
ImpDatRepl <- read.csv('./input/ImpDatRepl.csv', header=TRUE) # from 2_impute_data
# extract only first iteration 
ImpDat <- ImpDatRepl[ImpDatRepl$irep==1,]
# extract the number of repetitions
# max(ImpDatRepl$irep) would be simplier
nrep <- ImpDatRepl$irep[nrow(ImpDatRepl)]
# create empty matricies for PC1, rotations and Proportions of variance
PC1 <- matrix(NA,(length(ImpDat[,1])),0)
Rotations <- matrix(NA,0,9)
PropVar <- matrix(NA,0,9)
for(irep in 1:nrep){
  print(irep)
  # extract variables: PolPop PolTerr CapPop levels government infrastr writing texts money
  ImpDat <- ImpDatRepl[ImpDatRepl$irep==irep,5:13] #Replace "5:13" with "5:8" and "9:13" for confirmatory analyses (changing "9" in "Rotations <- matrix(NA,0,9)" and "PropVar <- matrix(NA,0,9)" to "4" or "5" as appropriate)
  # perform PCA
  res <- prcomp(ImpDat, scale=TRUE)
  # print summary of PCA output
  print(summary(res))
  # model predictions
  PCA<-predict(res)
  # extract PC1 (PCA$PC1)
  PC1 <- cbind(PC1,PCA[,1])
  # extract rotation
  Rotations <- rbind(Rotations,res$rotation)
  # extract Proportion of Variance explained (this is the second row)
  PropVar <- rbind(PropVar,summary(res)$importance[2,])
}

write.csv(PC1, file="./temp/PC1.csv",  row.names=FALSE)
write.csv(Rotations, file="./temp/Rotations.csv",  row.names=TRUE)
write.csv(PropVar, file="./temp/PropVar.csv",  row.names=FALSE)

expect_equal(dim(PC1), c(414, nrep))
expect_equal(dim(Rotations), c(9 * nrep, 9))
expect_equal(dim(PropVar), c(nrep, 9))

# all look right

#rm(AggrDat, ImpDat, PCA, res)

##### gPC1.R
## Calculate PC1 trajectories +/- 2SD conf int [lower, mean, upper]
#setwd("C:/Users/Peter Turchin/Google Drive/2.Seshat/1.R/PCA-MI")
PC1 <- read.table("./temp/PC1.csv", sep=",", header=TRUE)

#PC1[,5] <- -(PC1[,5]) #NB: It is possible that the signs may be arbitrarily flipped in different PCA runs, in which case used this type of code to correct arbitrarily flipped signs

# Scale PC1 to between 0-1
PC1 <- (PC1-min(PC1))/max(PC1-min(PC1)) 
dat <- ImpDatRepl[ImpDatRepl$irep==1,]
# convert NGA to character from factor
dat[,1] <- as.character(dat[,1])
# convert PolID to character from factor
dat[,2] <- as.character(dat[,2])
# calculate mean rowise (only from PC1)
dat$Mean <- apply(PC1,1,mean)
# calculate lower limit of standard deviation (only from PC1)
dat$Lower <- dat$Mean - 2*apply(PC1,1,sd)
# calculate upper limit of standard deviation (only from PC1)
dat$Upper <- dat$Mean + 2*apply(PC1,1,sd)

output <- read.table("./input/MIoutput.csv", sep=",", header=TRUE) # from 2_impute_data
# extract only the following columns
# "NGA" "PolID" "Date" "PropCoded" "PolPop" "PolTerr" "CapPop" "AdmLev" "MilLev" "ReligLev"    "SettlHier" " ProfOfficer" "ProfSoldier" "ProfPriest" "FullTBur" "ExamSyst" "MeritProm"   "GovtBldg" "Court" "LegCode" "Judge" "Lawyer" "Irrigation" "WaterSuppl" "Market" "FoodStor"    "Road" "Bridge" "Canal" "Port" "Mine" "Courier" "PostStation" "PostService" "Mnemonic"    "NonWRecord" 
dt <- output[,1:length(dat)]
# convert NGA to character from factor
dt[,1] <- as.character(dt[,1])
# convert PolID to character from factor
dt[,2] <- as.character(dt[,2])
# remove all values from columns other than NGA, PolID and Date
dt[,4:length(dat)] <- NA
#names(dt) <- c("NGA","PolID","Time","Mean","Lower","Upper")
names(dt) <- names(dat)

# Find years present in MIOutput which are missing in ImpDatRepl and fill with NA
for(i in 1:length(dat[,1])){
  # nrow(dt) - 414
  j <- 1:length(dt[,1])
  # match NGA, PolID and Date between ImpDatRepl.csv and MIoutput.csv
  j <- j[dt[,1]==dat[i,1] & dt[,2]==dat[i,2] & dt[,3]==dat[i,3]]
  dt[j,] <- dat[i,]
}
# for rows with a missing PropCoded in ImpDatRepl (only variables present are NGA, PolID and Date), replace missing values with those from closest previous date without missing values
for(i in 1:length(dt[,1])){   if(is.na(dt[i,4])){dt[i,4:length(dat)] <- dt[i-1,4:length(dat)]} }
dat <- dt

write.csv(dat, file="./temp/PC1_traj.csv",  row.names=FALSE)

expect_equal(dim(dat), c(817, 36))


################################################

### Graphs PCA results for Figure 3 of the PCA article. (see Figure3.R for more code) 
###    and for the Suppl Materials: all regions
#setwd("C:/Users/Peter Turchin/Google Drive/2.Seshat/1.R/PCA-MI")

png("./temp/figure3.png", units = "in", height = 5, width = 8, res = 300)
dat <- read.table("./temp/PC1_traj.csv", sep=",", header=TRUE)
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

dev.off()

################################################
#Merge PC1_traj.csv with other key info
dat <- read.table("./temp/PC1_traj.csv", sep=",", header=TRUE)
dat2 <- read.table("./input/polities.csv", sep=",", header=TRUE)
# merge PC1 trajectory with metadata from polities.csv
dat <- merge(dat,dat2,by="PolID")
# merge PC1 trajectory with metadata from NGA coordinates
dat2 <- read.table("./input/NGAcoords.csv", sep=",", header=TRUE)
dat <- merge(dat,dat2,by.x="NGA.y",by.y="NGA")
# rename NGA.y to NGA
dat<-plyr::rename(dat, c("NGA.y" = "NGA"))
# rename NGA.x to Originial NGA (the first NGA?)
dat<-plyr::rename(dat, c("NGA.x" = "OriginalNGA"))
# reorder by NGA and Time
dat<-dat[order(dat$NGA, dat$Time),]
# remove NGA, Polity, Time with centuries before the start date or after the end date (there are none)
dat<-subset(dat,dat$Time>=dat$Start & dat$Time<=dat$End)
write.csv(dat, file="./temp/PC1_traj_merged.csv",  row.names=FALSE)

expect_equal(dim(dat), c(864, 49))

dir_init("./output")

file.copy("./temp/PC1_traj_merged.csv", "./output")
