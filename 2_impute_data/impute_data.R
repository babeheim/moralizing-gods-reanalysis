
rm(list = ls())
source("../project_support.R")

dir_init("./temp")

Vars <- as.matrix(read.csv('./input/variables.csv', header=TRUE))
Vars <- Vars[Vars[,6]==Section1 | Vars[,6]==Section2 | Vars[,6]==Section3,] # Reduce the variables list to the Section set above

Vars[,1] <- paste(Vars[,2],Vars[,1]) #Creating unique variable/section combinations

polities <- read.csv('./input/polities.csv', header=TRUE)
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

write.csv(polities, file="./temp/polities.csv",  row.names=FALSE)

polities <- read.csv('./temp/polities.csv', header=TRUE)

NGAs <- levels(polities$NGA)

nrep <- 20
ImpDatRepl <- matrix(NA, nrow=0, ncol=0) 
for(irep in 1:nrep){
  print(irep)
  source("./code/ConstrMI.R") # this also writes "MIoutput.csv" which is passed on
  source("./code/AggrMI.R")
  source("./code/ImputeMI.R")
  ones <- matrix(data=1,nrow=length(AggrDat[,1]),ncol=1)
  colnames(ones) <- "irep"
  ImpDat <- cbind(AggrDat[,1:4],ImpDat,(ones*irep),AggrDat[,14:32])
  ImpDatRepl <- rbind(ImpDatRepl,ImpDat)
}

# Remove polity-dates that didn't yield 20 repl and post-colonial polities that
# couldn't be removed from multiple imputation due to bugs with only 1 polity/NGA

polities <- read.csv('./temp/polities.csv', header=TRUE)
polities <- polities[polities$PolID != "InGaroL",] #removing here because it caused bugs earlier
polities <- polities[polities$PolID != "CnHChin",] #removing here because it caused bugs earlier
polities <- polities[polities$PolID != "PgOrokL",] #removing here because it caused bugs earlier
write.csv(polities, file="./temp/polities.csv",  row.names=FALSE) 

#removing here because it seemed to create bugs when you have only 1 polity in an NGA, so couldn't remove earlier
ImpDatRepl <- ImpDatRepl[ImpDatRepl$PolID != "InGaroL",]
ImpDatRepl <- ImpDatRepl[ImpDatRepl$PolID != "CnHChin",]
ImpDatRepl <- ImpDatRepl[ImpDatRepl$PolID != "PgOrokL",]


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

write.csv(ImpDatRepl, file="./temp/ImpDatRepl.csv",  row.names=FALSE)

dir_init("./output")

file.copy(c("./temp/ImpDatRepl.csv", "./temp/MIoutput.csv"), "./output")

