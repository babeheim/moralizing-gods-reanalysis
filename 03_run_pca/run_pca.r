
rm(list = ls())
source("../project_support.r")

dir_init("./temp")

# Run PCA on ImpDat
ImpDatRepl <- read.csv('./input/ImpDatRepl.csv', header=TRUE) # from 2_impute_data
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


PC1 <- (PC1-min(PC1))/max(PC1-min(PC1)) # Scale PC1 to between 0-1
dat <- ImpDatRepl[ImpDatRepl$irep==1,]
dat[,1] <- as.character(dat[,1])
dat[,2] <- as.character(dat[,2])
dat$Mean <- apply(PC1,1,mean)
dat$Lower <- dat$Mean - 2*apply(PC1,1,sd)
dat$Upper <- dat$Mean + 2*apply(PC1,1,sd)

output <- read.table("./input/MIoutput.csv", sep=",", header=TRUE) # from 2_impute_data
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
dat <- merge(dat,dat2,by="PolID")
dat2 <- read.table("./input/NGAcoords.csv", sep=",", header=TRUE)
dat <- merge(dat,dat2,by.x="NGA.y",by.y="NGA")
dat<-plyr::rename(dat, c("NGA.y" = "NGA"))
dat<-plyr::rename(dat, c("NGA.x" = "OriginalNGA"))
dat<-dat[order(dat$NGA, dat$Time),]
dat<-subset(dat,dat$Time>=dat$Start & dat$Time<=dat$End)
write.csv(dat, file="./temp/PC1_traj_merged.csv",  row.names=FALSE)

expect_equal(dim(dat), c(864, 49))

dir_init("./output")

file.copy("./temp/PC1_traj_merged.csv", "./output")
