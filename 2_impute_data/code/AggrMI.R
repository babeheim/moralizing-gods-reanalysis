# Aggregate.R for MI
# Aggregate data into few categories
# scale, levels, government = (professions, bureaucracy, law), infrastr, writing, texts, money
# setwd("C:/Users/Peter Turchin/Google Drive/2.Seshat/1.R/PCA-MI")

data <- read.table('MIoutput.csv', sep=",", header=TRUE)
data <- data[data$PropCoded > 30,] # Omit sparsely coded polities with PropCoded < 30%
row.names(data) <- NULL
dat <- data[,5:96]
for(i in 1:52){dat[,i] <- as.numeric(dat[,i])}

# Construct aggregated data
AggrDat <- matrix(NA, 0, 0)
for(i in 1:nrow(dat)){
   row <- dat[i,1:3]                                                # scale variables, not-logged
   dt <- dat[i,c(4:7,52)]
   for(j in length(dt)){ row <- cbind(row,mean(dt[is.na(dt)==FALSE])) } # levels + HeredStatus
   dt <- dat[i,8:18]
   for(j in length(dt)){ row <- cbind(row,mean(dt[is.na(dt)==FALSE])) } # government
   dt <- dat[i,19:30]
   for(j in length(dt)){ row <- cbind(row,mean(dt[is.na(dt)==FALSE])) } # infrastr
   dt <- dat[i,31:37]
   for(j in length(dt)){ row <- cbind(row,mean(dt[is.na(dt)==FALSE])) } # writing
   dt <- dat[i,38:45]
   for(j in length(dt)){ row <- cbind(row,mean(dt[is.na(dt)==FALSE])) } # texts
   money <- as.numeric(dat[i,46:51])*1:6                                # money
      money <- money[is.na(money)==FALSE]
      row <- cbind(row,NA)
      if(length(money)!=0){row[9] <- max(money) }
      dt <- dat[i,52:54]
   for(j in length(dt)){ row <- cbind(row,mean(dt[is.na(dt)==FALSE])) } # Largest: freq 
   dt <- dat[i,55:56]
   for(j in length(dt)){ row <- cbind(row,mean(dt[is.na(dt)==FALSE])) } # conformity
   dt <- dat[i,57:59]
   for(j in length(dt)){ row <- cbind(row,mean(dt[is.na(dt)==FALSE])) } # Widespread: freq
   dt <- dat[i,60:61]
   for(j in length(dt)){ row <- cbind(row,mean(dt[is.na(dt)==FALSE])) } # conformity
   dt <- dat[i,62:64]
   for(j in length(dt)){ row <- cbind(row,mean(dt[is.na(dt)==FALSE])) } # Frequent: freq
   dt <- dat[i,65:66]
   for(j in length(dt)){ row <- cbind(row,mean(dt[is.na(dt)==FALSE])) } # conformity
   dt <- dat[i,67:69]
   for(j in length(dt)){ row <- cbind(row,mean(dt[is.na(dt)==FALSE])) } # Euphoric: freq
   dt <- dat[i,70:71]
   for(j in length(dt)){ row <- cbind(row,mean(dt[is.na(dt)==FALSE])) } # conformity
   dt <- dat[i,72:74]
   for(j in length(dt)){ row <- cbind(row,mean(dt[is.na(dt)==FALSE])) } # Dysphoric: freq
   dt <- dat[i,75:76]
   for(j in length(dt)){ row <- cbind(row,mean(dt[is.na(dt)==FALSE])) } # conformity
   dt <- dat[i,81]
   for(j in length(dt)){ row <- cbind(row,mean(dt[is.na(dt)==FALSE])) } # Big Gods: MHG
   dt <- dat[i,c(83:84,87)]
   for(j in length(dt)){ row <- cbind(row,mean(dt[is.na(dt)==FALSE])) } # Big Gods: BSP
   dt <- dat[i,33]
   for(j in length(dt)){ row <- cbind(row,mean(dt[is.na(dt)==FALSE])) } # written records
      dt <- dat[i,6]
   for(j in length(dt)){ row <- cbind(row,mean(dt[is.na(dt)==FALSE])) } # religious hierarchy
   
   rowNA <- row
   rowNA[is.na(row)] <- -99999
   if(length(AggrDat)==0){AggrDat <- rbind(AggrDat,rowNA)}
   if(length(AggrDat)>0)
   {if(all(AggrDat[length(AggrDat[,1]),]==rowNA)==FALSE){AggrDat <- rbind(AggrDat,rowNA)} } 
} 
for(i in 1:length(AggrDat[,1])){
   for(j in 1:length(AggrDat[1,])){
      if(AggrDat[i,j]==-99999){AggrDat[i,j] <- NA}
   }}
AggrDat <- cbind( data[row.names(AggrDat),c(1,2,3,4)] ,AggrDat)
AggrDat[,1] <- as.character(AggrDat[,1])
AggrDat[,2] <- as.character(AggrDat[,2])
colnames(AggrDat) <- c("NGA", "PolID","Time", "PropCoded", "PolPop","PolTerr","CapPop","levels", "government", "infrastr", "writing", "texts", "money",   "FreqLR", "ConformityLR", "FreqWR",  "ConformityWR", "FreqFR", "ConformityFR", "FreqER", "ConformityER", "FreqDR", "ConformityDR", "MHG", "BSP", "WrittenRecords","ReligiousHier")
row.names(AggrDat) <- NULL

AggrDat$GeneralMoralisticPunishment<-ifelse(AggrDat$BSP>0.1,1,0)
AggrDat$MoralisingHighGods<-ifelse(AggrDat$MHG==3,1,0)
AggrDat$MoralisingGods<-ifelse(AggrDat$GeneralMoralisticPunishment>0.1 | AggrDat$MoralisingHighGods>0.1 ,1,0)
AggrDat$DoctrinalMode<-ifelse(AggrDat$FreqLR>4.5 | AggrDat$FreqWR>4.5 | AggrDat$FreqFR>4.5 | AggrDat$FreqER>4.5 | AggrDat$FreqDR>4.5 | AggrDat$ReligiousHier>=2,1,0)
AggrDat$Writing<-ifelse(AggrDat$WrittenRecords>0.1,1,0)

write.csv(AggrDat, file="MIAggrDat.csv",  row.names=FALSE)

rm(dat,data,dt,row,rowNA,i,j,money)

