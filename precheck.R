# Checks Soc Complx data from exportdat.csv for error
polities <- read.csv('polities.csv', header=TRUE)
polities <- polities[polities$Dupl=="n",]
Vars <- as.matrix(read.csv('variables.csv', header=TRUE))
SCdat <- matrix(nrow = 0, ncol = 0)
dat <- read.table('exportdat.csv', sep=",", header=TRUE, quote = "", colClasses = "character")
dat <- dat[dat$Section==Section1 | dat$Section==Section2 | dat$Section==Section3,] # Section is set in !MI.R

Vars[,1] <- paste(Vars[,2],Vars[,1]) #Creating unique variable/section combinations
dat[,5] <- paste(dat[,4],dat[,5]) #Creating unique variable/section combinations

for(i in 1:length(Vars[,1])){
   var <- Vars[i,1]
   dt <- dat[dat$Variable==var,]
   SCdat <- rbind(SCdat,dt)
}
dat <- SCdat
SCdat <- matrix(nrow = 0, ncol = 0)

for(i in 1:nrow(polities)){
   dt <- dat[dat$Polity==polities$PolID[i],]
   SCdat <- rbind(SCdat,dt)
}

SCdat <- SCdat[,c(1,2,5,6,7,8,9,10,11,12)]
row.names(SCdat) <- NULL

# Convert categorical values to numbers. Ignore warnings: they will be taken care off in the next step -- in errors
for(i in 1:nrow(SCdat)){
   for(j in 4:5){
       if(SCdat[i,j] == "present"){SCdat[i,j] <- "1"}   
      if(SCdat[i,j] == "inferred present"){SCdat[i,j] <- "0.9"}   
      if(SCdat[i,j] == "inferred absent"){SCdat[i,j] <- "0.1"}   
      if(SCdat[i,j] == "absent"){SCdat[i,j] <- "0"}   
      if(SCdat[i,j] == "none"){SCdat[i,j] <- "0"}
      if(SCdat[i,j] == "daily"){SCdat[i,j] <- "6"}   
      if(SCdat[i,j] == "weekly"){SCdat[i,j] <- "5"}
      if(SCdat[i,j] == "monthly"){SCdat[i,j] <- "4"}  
      if(SCdat[i,j] == "seasonally"){SCdat[i,j] <- "3"}  
      if(SCdat[i,j] == "yearly"){SCdat[i,j] <- "2"}   
      if(SCdat[i,j] == "once per generation"){SCdat[i,j] <- "1"}   
      if(SCdat[i,j] == "once in a lifetime"){SCdat[i,j] <- "0"}    
      if(SCdat[i,j] == "whole polity"){SCdat[i,j] <- "3"}      
      if(SCdat[i,j] == "majority"){SCdat[i,j] <- "2"}      
      if(SCdat[i,j] == "substantial minority"){SCdat[i,j] <- "1"}      
      if(SCdat[i,j] == "elites"){SCdat[i,j] <- "0"}          
      if(SCdat[i,j] == "inactive"){SCdat[i,j] <- "1"}
      if(SCdat[i,j] == "active"){SCdat[i,j] <- "2"}          
      if(SCdat[i,j] == "moralizing"){SCdat[i,j] <- "3"}          
      if(SCdat[i,j] == "monotheistic"){SCdat[i,j] <- "1"}     
      if(SCdat[i,j] == "polytheistic"){SCdat[i,j] <- "2"}                                                     
      if(SCdat[i,j] == "not applicable"){SCdat[i,j] <- "unknown"}      
      if(SCdat[i,j] == "suspected unknown"){SCdat[i,j] <- "unknown"}      
      if(SCdat[i,j] == "unknown"){SCdat[i,j] <- NA}     
   }}
SCdat <- SCdat[is.na(SCdat[,4])==FALSE,]
dat <- SCdat
for(i in 1:nrow(SCdat)){
   dat[i,4] <- as.numeric(SCdat[i,4])
}
errors <- SCdat[is.na(dat[,4]),]

SCdat <- SCdat[is.na(dat[,4])==FALSE,]
dat <- SCdat[SCdat[,5]!="",]
datNA <- dat
for(i in 1:nrow(dat)){
   datNA[i,5] <- as.numeric(dat[i,5])
}
errors <- rbind(errors,dat[is.na(datNA[,5]),])

# Change BCE to negative years and remove CE. 
for(i in 1:nrow(SCdat)){
   if(substr(SCdat[i,6], (nchar(SCdat[i,6]) - 2) , nchar(SCdat[i,6]) ) =="BCE" )
   {a <- -as.numeric(substr(SCdat[i,6], 1, (nchar(SCdat[i,6]) - 3)))
    SCdat[i,6] <- a}
   if(substr(SCdat[i,7], (nchar(SCdat[i,7]) - 2) , nchar(SCdat[i,7]) ) =="BCE" )
   {a <- -as.numeric(substr(SCdat[i,7], 1, (nchar(SCdat[i,7]) - 3)))
    SCdat[i,7] <- a}
}
for(i in 1:nrow(SCdat)){
   if(substr(SCdat[i,6], (nchar(SCdat[i,6]) - 1) , nchar(SCdat[i,6]) ) =="CE" )
   {a <- as.numeric(substr(SCdat[i,6], 1, (nchar(SCdat[i,6]) - 2)))
    SCdat[i,6] <- a}
   if(substr(SCdat[i,7], (nchar(SCdat[i,7]) - 1) , nchar(SCdat[i,7]) ) =="CE" )
   {a <- as.numeric(substr(SCdat[i,7], 1, (nchar(SCdat[i,7]) - 2)))
    SCdat[i,7] <- a}
}

dat <- SCdat[SCdat[,6]!="",]
for(i in 1:nrow(dat)){
   dat[i,6] <- as.numeric(dat[i,6])
}
errors <- rbind(errors,dat[is.na(dat[,6]),])

dat <- SCdat[SCdat[,7]!="",]
for(i in 1:nrow(dat)){
   dat[i,7] <- as.numeric(dat[i,7])
}
errors <- rbind(errors,dat[is.na(dat[,7]),])

write.csv(errors, file="errors.csv",  row.names=FALSE)
write.csv(SCdat, file="SCdat.csv",  row.names=FALSE)

rm(i,j,var,dat,dt,datNA,a,polities,Vars)

## Checks for polity inclusions
#dat <- read.table('exportdat.csv', sep=",", header=TRUE, quote = "", colClasses = "character")
#dat <- dat[dat$Section==Section,]
#polities_wiki <- levels(as.factor(dat$Polity))
#polities <- read.csv('polities.csv', header=TRUE)
#polities <- polities[polities$Dupl=="n",]
#polities <- levels(as.factor(polities$PolName))

#extras <- vector(length=0)
#for (i in 1:length(polities_wiki)){
#   if(any(polities_wiki[i] == polities)==FALSE){extras <- c(extras, polities_wiki[i])}
#}
#missing <- vector(length=0)
#for (i in 1:length(polities)){
#   if(any(polities[i] == polities_wiki)==FALSE){missing <- c(missing, polities[i])}
#}
#write.csv(extras, file="extras.csv",  row.names=FALSE)
#write.csv(missing, file="missing.csv",  row.names=FALSE)








