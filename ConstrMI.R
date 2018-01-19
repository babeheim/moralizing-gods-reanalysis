# Multiple Imputation  
# Uses CSdat.csv from precheck.R
#  
# setwd("C:/Users/Peter Turchin/Google Drive/2.Seshat/1.R/PCA-MI")
output <- matrix(nrow = 0, ncol = (4+nrow(Vars)))
   
for(iNGA in 1:(0+1*length(NGAs))){
   NGA <- NGAs[iNGA]
   dat <- read.table('SCdat.csv', sep=",", header=TRUE, colClasses = "character")
   polities <- read.csv('polities.csv', header=TRUE)
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
# Randomly sample from ranges assuming 90% interval, Normal distribution
for(i in 1:nrow(dat)){
   if(dat$Value.Note[i]=="range") {
      mean <- 0.5*(as.numeric(dat[i,4])+as.numeric(dat[i,5]))
      sd <- abs((as.numeric(dat[i,5]) - as.numeric(dat[i,4]))/(2*1.645))
      dat[i,4] <- rnorm(1,mean = mean, sd = sd)
   }}
# Randomly sample from disputed and uncertain, eliminate extra rows
for(i in 1:(nrow(dat)-1)){
   if(dat$Value.Note[i]=="disputed" | dat$Value.Note[i]=="uncertain"){value <- as.numeric(dat[i,4])
      for(j in (i+1):nrow(dat)){
         if(dat[i,1] == dat[j,1] & dat[i,2] == dat[j,2] & dat[i,3] == dat[j,3] & dat[i,6] == dat[j,6] & dat[i,7] == dat[j,7]){
            value <- c(value,as.numeric(dat[j,4]))
            dat[j,9] <- "delete"}}
            dat[i,4] <- sample(value, size=1)         
   }  }   
dat <- dat[dat[,9] != "delete",]
row.names(dat) <- NULL
datSC <- dat[,c(2,3,4,6,7)]
colnames(datSC) <- c("PolID","Variable","Value", "Date", "DateTo")

# Construct output
tmin <- ceiling(0.01*min(polities$Start[polities$NGA==NGA]))
tmax <- floor(0.01*max(polities$End[polities$NGA==NGA]))
out <- matrix(nrow=c(length(100*tmin:tmax)),ncol=(3+nrow(Vars)))
colnames(out) <- c("NGA","PolID", "Date", Vars[,3])  # Use short names for variables
out[,1] <- as.character(NGA)
out[,3] <- 100*tmin:tmax

for(i in 1:nrow(out)){ 
   for(j in 1:nrow(polities)){
      if( (as.numeric(out[i,3]) <= as.numeric(polities[j,5])) & 
             (as.numeric(out[i,3]) >= as.numeric(polities[j,4])) ){ 
         out[i,2] <- as.character(polities[j,3]) 
      }}}
out <- out[is.na(out[,2])==FALSE,]   # Eliminate centuries for which a polity is lacking

# First populate 'out' with data tied to polities, not dates
for(ivar in 1:nrow(Vars)){
   datV <- datSC[(datSC$Variable==Vars[ivar]) & (datSC$Date==""),]
   if(is.null(nrow(datV))){datV <- array(datV,c(1,5))}
   for(i in 1:nrow(datV)){
      for(j in 1:nrow(out)){
         if(nrow(datV) != 0){
            if(out[j,2] == datV$PolID[i]){out[j,ivar+3] <- datV$Value[i]
            }}}}}

# Next populate 'out' with data tied to a single date
for(ivar in 1:nrow(Vars)){
   datV <- datSC[((datSC$Variable==Vars[ivar]) & (datSC$Date!="") & (datSC$DateTo=="")),]
   if(is.null(nrow(datV))){datV <- array(datV,c(1,5))}
   for(i in 1:nrow(datV)){
      for(j in 1:nrow(out)){
         if(nrow(datV) != 0){
            century <- 100*round(as.numeric(datV[i,4])/100)
            if(out[j,3] == as.character(century)){out[j,ivar+3] <- datV$Value[i]
            }}}}}


# Finally populate 'out' with data tied to a range of dates
for(ivar in 1:nrow(Vars)){
   datV <- datSC[((datSC[,2]==Vars[ivar]) & (datSC[,4]!="") & (datSC[,5]!="")),]
   if(is.null(nrow(datV))){datV <- array(datV,c(1,5))}
   for(i in 1:nrow(datV)){
      for(j in 1:nrow(out)){
         if(nrow(datV) != 0){
            century <- as.numeric(out[j,3])
            tmin <- as.numeric(datV[i,4])
            tmax <- as.numeric(datV[i,5])
            if(century >= tmin & century <= tmax){out[j,ivar+3] <- datV[i,3]
            }}}}}

# Calculate the proportion of data coded by century
PropCoded <- array(0,c(nrow(out),2))
PropCoded[,1] <- out[,3]
colnames(PropCoded) <- c("Date1","PropCoded")
for(i in 1:nrow(out)){
  j <- 0  
  for(ivar in 1:nrow(Vars[1:51,])){
       if(is.na(out[i,ivar+3])){j <- j+1} 
    }
  PropCoded[i,2] <- 0.1*round((nrow(Vars[1:51,]) - j)/nrow(Vars[1:51,])*1000) # Keep 3 sign digit; NB: THIS CODES ONLY PROPORTION OF SOCIAL COMPLEXITY [ROWS 1:51] FOR CONSISTENCY s
}
out <- cbind(out[,1:3],PropCoded[,2],out[,4:(nrow(Vars)+3)])
colnames(out) <- c("NGA","PolID", "Date", "PropCoded", Vars[,3])  
output <- rbind(output,out)
}  # Closing the iNGA loop

output <- output[output[,4] != 0,]  # Remove rows of missing values
write.csv(output, file="MIoutput.csv",  row.names=FALSE)

rm(dat_temp,dat,datSC,datV,out,PropCoded,century,i,iNGA,ivar,j,mean,NGA,sd,tmax,tmin,value)


