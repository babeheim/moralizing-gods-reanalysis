
rm(list = ls())
source("../project_support.r")

dir_init("./temp")

Vars <- as.matrix(read.csv('./input/variables.csv', header=TRUE))
# Filter variables by sections of interest: "Social Complexity variables", "Ritual variables", "Religion and Normative Ideology
Vars <- Vars[Vars[,6]==Section1 | Vars[,6]==Section2 | Vars[,6]==Section3,] # Reduce the variables list to the Section set above

#Creating unique variable/section combinations
Vars[,1] <- paste(Vars[,2],Vars[,1])

# Filter polities to remove polities that are not being studied
polities <- read.csv('./input/polities.csv', header=TRUE)
#### Remove new NGAs
polities <- polities[polities$NGA != "Crete",]
polities <- polities[polities$NGA != "Galilee",]
polities <- polities[polities$NGA != "Middle Ganga",]
#### Remove low-coverage polities causing bugs in C. Java)
polities <- polities[polities$PolID != "IdCJBun",]
polities <- polities[polities$PolID != "IdKalin",] 
polities <- polities[polities$PolID != "EsHabsb",]
#### Remove post-colonial polities (also see below for removing post-colonial polities from NGAs with only 2 polities)
polities <- polities[polities$PolID != "USIllinL",]
polities <- polities[polities$PolID != "GbEmpir",]
polities <- polities[polities$PolID != "InEInCo",]
polities <- polities[polities$PolID != "InBritP",]
polities <- polities[polities$PolID != "RuYakuL",]
polities <- polities[polities$PolID != "UsIroqL",]

expect_equal(dim(polities), c(335, 11))

write.csv(polities, file="./temp/polities.csv",  row.names=FALSE)
polities <- read.csv('./temp/polities.csv', header=TRUE)
# Extract Natural Geographic Areas (NGAs)
NGAs <- levels(polities$NGA)

expect_equal(length(NGAs), 30)

# `nrep` (number of imputations) now defined in project_support.r

# create empty matrix
ImpDatRepl <- matrix(NA, nrow=0, ncol=0) 
for(irep in 1:nrep){
  print(paste("imputation", irep))
  # transpose data from long to wide and randomly select from disputed and range values
  source("./code/1ConstrMI.R")
 # Aggregate variables
  source("./code/2AggrMI.R")
  # impute missing aggregated variables 
  source("./code/3ImputeMI.R")
 # create matrix of ones for storing imputation number
  ones <- matrix(data=1,nrow=length(AggrDat[,1]),ncol=1)
  colnames(ones) <- "irep"
  # bind row of aggregated and imputed data for each NGA, polity, time
  ImpDat <- cbind(AggrDat[,1:4],ImpDat,(ones*irep),AggrDat[,14:32])
  # bind individual NGA, polity, times into one matrix
  ImpDatRepl <- rbind(ImpDatRepl,ImpDat)
}

expect_equal(dim(ImpDatRepl), c(417 * nrep, 33))

####### Remove polity-dates that didn't yield 20 repl #and post-colonial polities that couldn't be removed from multiple imputation due to bugs with only 1 polity/NGA

polities <- read.csv('./temp/polities.csv', header=TRUE)
polities <- polities[polities$PolID != "InGaroL",] #removing here because it caused bugs earlier
write.csv(polities, file="./temp/polities.csv",  row.names=FALSE) 
polities <- polities[polities$PolID != "CnHChin",] #removing here because it caused bugs earlier
write.csv(polities, file="./temp/polities.csv",  row.names=FALSE) 
polities <- polities[polities$PolID != "PgOrokL",] #removing here because it caused bugs earlier
write.csv(polities, file="./temp/polities.csv",  row.names=FALSE) 

expect_equal(dim(polities), c(332, 11))

# remove the same polities from the imputed data
ImpDatRepl <- ImpDatRepl[ImpDatRepl$PolID != "InGaroL",] #removing here because it seemed to create bugs when you have only 1 polity in an NGA, so couldn't remove earlier
ImpDatRepl <- ImpDatRepl[ImpDatRepl$PolID != "CnHChin",] #removing here because it seemed to create bugs when you have only 1 polity in an NGA, so couldn't remove earlier
ImpDatRepl <- ImpDatRepl[ImpDatRepl$PolID != "PgOrokL",] #removing here because it seemed to create bugs when you have only 1 polity in an NGA, so couldn't remove earlier

# create tempory data set
dat_temp <- ImpDatRepl
# remove duplicate rows for the same imputation repetition 
for(i in 1:nrow(polities)){
  # split data by polity
  dat <- ImpDatRepl[as.character(ImpDatRepl[,2])==as.character(polities[i,2]),]
  # if the polity has at least 1 row
  if(nrow(dat)!=0){
    # extract vector of unique times
    Time <- unique(dat$Time)
    for(j in 1:length(Time)){
      # split data by unique time
      dt <- dat[dat$Time==Time[j],]
      # if the number of imputation repetitions is not 20
      if(nrow(dt) != nrep){
        # print number of rows of dt
        print(nrow(dt))
        # print the NGA, Polity and Time
        print(dt[1,1:3])
        # mark duplicate rows with the same irep for deletion
        dat_temp[as.character(dat_temp$PolID)==as.character(dat$PolID[1]) & dat_temp$Time==Time[j],14] <- -99999
      }
    }
  }
}
# delete duplicate rows
ImpDatRepl <- dat_temp[dat_temp$irep!=-99999,]

write.csv(ImpDatRepl, file="./temp/ImpDatRepl.csv",  row.names=FALSE)

# 8280 rows in the original analysis, but it changes with nrep
# expect_equal(dim(ImpDatRepl), c(8280, 33))

dir_init("./output")

files <- c("./temp/polities.csv", "./temp/ImpDatRepl.csv", "./temp/MIoutput.csv")
file.copy(files, "./output")

print("multiple-imputation complete; ImpDatRepl.csv, MIoutput.csv created; polities.csv updated")
