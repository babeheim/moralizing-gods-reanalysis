
output <- matrix(nrow = 0, ncol = (4+nrow(Vars)))
   
for(iNGA in 1:(0+1*length(NGAs))){
   NGA <- NGAs[iNGA]
   dat <- read.table('./input/SCdat.csv', sep=",", header=TRUE, colClasses = "character")
   polities <- read.csv('./temp/polities.csv', header=TRUE)
   polities <- polities[polities[,1]==NGA,]      # Use only one NGA at a time
   polities <-polities[polities[,8]=="n",]       # Exclude duplicates
   row.names(polities) <- NULL
   dat <- dat[dat$NGA==NGA,]                     # Select data for the NGA
   row.names(dat) <- NULL
   dat_temp <- matrix(nrow = 0, ncol = 9)        # Make sure all data are for polities in polities.csv
   for(i in 1:nrow(polities)){
      dat_temp <- rbind(dat_temp,dat[dat$Polity == polities[i,2],])
      }
   dat <-dat_temp

# Group by NGA and polity, and randomly sample from variables with ranges assuming 90% intervals
# "Values that are coded with a confidence interval are sampled from a Gaussian distribution, with mean and variance that are estimated assuming that the interval covers 90% of the probability"
# Variables are: Social Scale Polity Population, Social Scale Polity territory, Social Scale Population of the largest settlement, Hierarchical Complexity Military levels, Hierarchical Complexity Religious levels, Hierarchical Complexity Administrative levels, Hierarchical Complexity Settlement hierarchy
# Randomly sample from ranges assuming 90% interval, Normal distribution
for(i in 1:nrow(dat)){
  # filter by range Value.Note
   if(dat$Value.Note[i]=="range") {
    # mean = 0.5 * (Value.From+Value.To)
      mean <- 0.5*(as.numeric(dat[i,4])+as.numeric(dat[i,5]))
    # sd = abs((Value.To-Value.From)/(2*1.645)
      sd <- abs((as.numeric(dat[i,5]) - as.numeric(dat[i,4]))/(2*1.645))
    # sample from Gaussian distrubition, using mean and sd, assuming that the range covers 90% of the probability
      dat[i,4] <- rnorm(1,mean = mean, sd = sd)
   }}
# Randomly sample from disputed and uncertain value notes and eliminate extra rows
# "For categorical or binary variables, we sample coded values in proportion to the number of categories that are presented as plausible."
# "In cases where experts disagree, each alternative coding has the same probability of being selected."
# Variables included are: 
# Disputed - Binary: Specialized Buildings: polity owned food storage site, Other Indigenous coins, Information Mnemonic devices, Professions Professional military officers, Bureaucracy characteristics Full-time bureaucrats, Bureaucracy characteristics Merit promotion, Law Judges, Professions Professional soldiers, Specialized Buildings: polity owned irrigation systems, Bureaucracy characteristics Examination system, Law Formal legal code, Information Script, Information Written records, Law Courts, Information Scientific literature, Specialized Buildings: polity owned drinking water supply systems, Information Philosophy, Specialized Buildings: polity owned markets
# Disputed - Continuous: Social Scale Polity Population, Social Scale Population of the largest settlement, Social Scale Polity territory
# Uncertain - all binary/ categorical
for(i in 1:(nrow(dat)-1)){
  # filter by disputed or uncertain Value.Note
   if(dat$Value.Note[i]=="disputed" | dat$Value.Note[i]=="uncertain"){value <- as.numeric(dat[i,4])
      for(j in (i+1):nrow(dat)){
        # group by NGA, polity, newvariable, Date.From and Date.To
         if(dat[i,1] == dat[j,1] & dat[i,2] == dat[j,2] & dat[i,3] == dat[j,3] & dat[i,6] == dat[j,6] & dat[i,7] == dat[j,7]){
           # create array of disputed/uncertain values
            value <- c(value,as.numeric(dat[j,4]))
            # mark duplicate rows for deletion
            dat[j,9] <- "delete"}}
            # sample from disputed/uncertain variables
            # this effectively works like 
            # sample(c(1000,2000), 1) will always produce an output or 1000 or 2000
            dat[i,4] <- sample(value, size=1)         
   }  }   
# remove duplicate rows to be deleted
dat <- dat[dat[,9] != "delete",]
# remove row names
row.names(dat) <- NULL
# extract only relevant variables
# Polity, Variable, Value.From, Date.From, Date.To
datSC <- dat[,c(2,3,4,6,7)]
colnames(datSC) <- c("PolID","Variable","Value", "Date", "DateTo")

# Construct output
# calculate the starting century
tmin <- ceiling(0.01*min(polities$Start[polities$NGA==NGA]))
# calculate the ending century
tmax <- floor(0.01*max(polities$End[polities$NGA==NGA]))
# combine seshat variables with short variable names (Vars$ShortName)
out <- matrix(nrow=c(length(100*tmin:tmax)),ncol=(3+nrow(Vars)))
# rename columns
# NGA, PolID, Date, PolPop, PolTerr, CapPop, etc. ... HeredStatus
# It must be noted that there are two variables with the short name SEAmoral, Normative Ideological precepts concerning morality/supernatural beings other supernatural enforcement of piety, Normative Ideological precepts concerning morality/supernatural beings other supernatural enforcement not related to piety
colnames(out) <- c("NGA","PolID", "Date", Vars[,3])  # Use short names for variables
# Populate the first column of out with NGA
out[,1] <- as.character(NGA)
# Populate Date with centuries, calculated by the range of each polity
out[,3] <- 100*tmin:tmax

# Polpulate the PolID variable of out with PolIDs from the polities dataset
for(i in 1:nrow(out)){ 
   for(j in 1:nrow(polities)){
  # if out$Date <= polities$End and out$Date >= polities$Start
      if( (as.numeric(out[i,3]) <= as.numeric(polities[j,5])) & 
             (as.numeric(out[i,3]) >= as.numeric(polities[j,4])) ){ 
     # replace out$PolID with polities$PolID
         out[i,2] <- as.character(polities[j,3]) 
      }}}

# Eliminate centuries for which a polity is lacking
out <- out[is.na(out[,2])==FALSE,]

# First populate 'out' with data tied to polities, not dates
for(ivar in 1:nrow(Vars)){
  # split data by variable per polity
   datV <- datSC[(datSC$Variable==Vars[ivar]) & (datSC$Date==""),]
   if(is.null(nrow(datV))){datV <- array(datV,c(1,5))}
   for(i in 1:nrow(datV)){
      for(j in 1:nrow(out)){
       # if split data has at least 1 row
         if(nrow(datV) != 0){
          # if out$PolID == datV$PolID, replace empty variable with value from datV
            if(out[j,2] == datV$PolID[i]){out[j,ivar+3] <- datV$Value[i]
            }}}}}

# Next populate 'out' with data tied to a single date
for(ivar in 1:nrow(Vars)){
  # split data by variable per polity
  datV <- datSC[((datSC$Variable==Vars[ivar]) & (datSC$Date!="") & (datSC$DateTo=="")),]
  if(is.null(nrow(datV))){datV <- array(datV,c(1,5))}
  for(i in 1:nrow(datV)){
    for(j in 1:nrow(out)){
      # if split data has at least 1 row
      if(nrow(datV) != 0){
        # calculate the century of the variable
        century <- 100*round(as.numeric(datV[i,4])/100)
        # if out$Date matches the century, replace empty variable with value from datV
        if(out[j,3] == as.character(century)){out[j,ivar+3] <- datV$Value[i]
        }}}}}


# Finally populate 'out' with data tied to a range of dates
for(ivar in 1:nrow(Vars)){
  # split data by variable per polity
  datV <- datSC[((datSC[,2]==Vars[ivar]) & (datSC[,4]!="") & (datSC[,5]!="")),]
  if(is.null(nrow(datV))){datV <- array(datV,c(1,5))}
  for(i in 1:nrow(datV)){
    for(j in 1:nrow(out)){
      # if split data has at least 1 row
      if(nrow(datV) != 0){
        # extract the Date (century) from out
        century <- as.numeric(out[j,3])
        # extract the start century from long data
        tmin <- as.numeric(datV[i,4])
        # extract the end century from long data
        tmax <- as.numeric(datV[i,5])
        # if the century is equal or larger than the start date and smaller than than or equal to the end date, replace empty variable with value from datV
        if(century >= tmin & century <= tmax){out[j,ivar+3] <- datV[i,3]
        }}}}}

# Calculate the proportion of data coded by century
#### NB: THIS CODES ONLY PROPORTION OF SOCIAL COMPLEXITY [ROWS 1:51] FOR CONSISTENCY s
# Extract column names of social complexity variables [1:51] which are used for calculating Propcoded
# Note: if the correct code was used to extract all social complexity variables
# variables$ShortName[variables$Section == "Social Complexity variables"]
# this would also include HeredStatus, however this code misses this variable
# create empty array
PropCoded <- array(0,c(nrow(out),2))
# the first column is out$Date
PropCoded[,1] <- out[,3]
colnames(PropCoded) <- c("Date1","PropCoded")
for(i in 1:nrow(out)){
  j <- 0  
  # lopp through only social complexity variables from wide data
  for(ivar in 1:nrow(Vars[1:51,])){
    if(is.na(out[i,ivar+3])){j <- j+1} 
  }
  # Calculate proportion coded
  PropCoded[i,2] <- 0.1*round((nrow(Vars[1:51,]) - j)/nrow(Vars[1:51,])*1000) # Keep 3 sign digit; NB: THIS CODES ONLY PROPORTION OF SOCIAL COMPLEXITY [ROWS 1:51] FOR CONSISTENCY s
}
# recombine PropCoded with wide data
out <- cbind(out[,1:3],PropCoded[,2],out[,4:(nrow(Vars)+3)])
colnames(out) <- c("NGA","PolID", "Date", "PropCoded", Vars[,3])  
output <- rbind(output,out)
}  # Closing the iNGA loop

# Remove rows of missing values
output <- output[output[,4] != 0,]
write.csv(output, file="./temp/MIoutput.csv",  row.names=FALSE)

rm(dat_temp,dat,datSC,datV,out,PropCoded,century,i,iNGA,ivar,j,mean,NGA,sd,tmax,tmin,value)
