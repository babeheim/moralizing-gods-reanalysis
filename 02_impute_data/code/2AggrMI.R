# Aggregate data into few categories:
# scale, levels, government = (professions, bureaucracy, law), infrastr, writing, texts, money

data <- read.table('./temp/MIoutput.csv', sep=",", header=TRUE)
data <- data[data$PropCoded > 30,]
# Omit sparsely coded polities with PropCoded < 30%
row.names(data) <- NULL
# remove NGA, polity, Date and PropCoded
dat <- data[,5:96]
# convert to numeric 
for(i in 1:52){dat[,i] <- as.numeric(dat[,i])}

# Construct aggregated data
AggrDat <- matrix(NA, 0, 0)
for(i in 1:nrow(dat)){
  # scale variables, not-logged
   row <- dat[i,1:3]
  # levels + HeredStatus "levels"
  # Next come measures of hierarchical or vertical complexity (“levels of hierarchy” in Figure 1). These focus on the number of control/decision levels in the administrative, religious, and military hierarchies. Another measure of vertical complexity is the number of levels in the settlement hierarchy. The four hierarchical variables were averaged to yield the “levels of hierarchy” variable.
  # Administrative levels = AdmLev
  # Military levels = MilLev
  # Religious levels = ReligLev
  # Settlement hierarchy = SettlHier
  # However, c(4:7,52) actually subsets AdmLev MilLev ReligLev SettlHier SpecFreqLR (Frequency for the ritual specialist)
  # Is this meant to be HeredStatus rather than SpecFreqLR? In which case the column index is 98, not SpecFreqLR
   dt <- dat[i,c(4:7,52)]
   for(j in length(dt)){ row <- cbind(row,mean(dt[is.na(dt)==FALSE])) } # levels + HeredStatus
  # “Government” variables code for presence or absence of professional soldiers and officers, priests, bureaucrats, and judges. This class also includes characteristics of the bureaucracy and of the judicial system, and presence of specialized buildings (e.g., courts). Government variables were aggregated by adding the number of binary codes indicating “present” and dividing them by the total number of variables. The aggregated variable Government, thus is scaled between 0 and 1.
  # government "government"
  # ProfOfficer ProfSoldier ProfPriest FullTBur ExamSyst MeritProm GovtBldg Court LegCode Judge Lawyer
   dt <- dat[i,8:18]
   for(j in length(dt)){ row <- cbind(row,mean(dt[is.na(dt)==FALSE])) } # government
  # The variety of public goods and public works provided by the community is captured in “Infrastructure.”
  # infrastr "infrastr"
  # Irrigation WaterSuppl Market FoodStor Road Bridge Canal Port Mine Courier PostStation PostService
   dt <- dat[i,19:30]
   for(j in length(dt)){ row <- cbind(row,mean(dt[is.na(dt)==FALSE])) } # infrastr
  # Informational complexity is coded by the characteristics of the writing and record-keeping (more generally, informational) systems. We also record whether the society has developed specialized literature, including history, philosophy, and fiction. These binary codes were treated the same way as Government, yielding aggregated variables Infrastructure, Writing, and Texts
  # writing "writing"
  # Mnemonic NonWRecord WRecord Script NonPhWrit PhAlph Lists
   dt <- dat[i,31:37]
   for(j in length(dt)){ row <- cbind(row,mean(dt[is.na(dt)==FALSE])) } # writing
  # texts "texts"
  # Calendar SacrTxt ReligLit PractLit History Philosophy SciLit Fiction
   dt <- dat[i,38:45]
   for(j in length(dt)){ row <- cbind(row,mean(dt[is.na(dt)==FALSE])) } # texts
  # Finally, the sophistication of the cash economy is reflected in Monetary System, which can take values between 0 and 6, reflecting the “most sophisticated” monetary instrument present in the coded society (Figure 1 in the main article). For example, if precious metals were used as money, while foreign and indigenous coins and paper currency were absent, Money would take the value of 3. If on the other hand, paper currency was present, the value of the aggregated variable is 6. 
  # money "money"
  # Article Token PrecMetal ForCoin IndigCoin PaperCurr
   money <- as.numeric(dat[i,46:51])*1:6                                # money
   money <- money[is.na(money)==FALSE]
   row <- cbind(row,NA)
   if(length(money)!=0){row[9] <- max(money) }
   # Largest: freq "FreqLR"
   # SpecFreqLR PartFreqLR AudFreqLR
   dt <- dat[i,52:54]
   for(j in length(dt)){ row <- cbind(row,mean(dt[is.na(dt)==FALSE])) } # Largest: freq 
  # conformity "ConformityLR"
  # OrthopraxyLR OrthodoxyLR
   dt <- dat[i,55:56]
   for(j in length(dt)){ row <- cbind(row,mean(dt[is.na(dt)==FALSE])) } # conformity
  # Widespread: freq "FreqWR"
  # SpecFreqWR PartFreqWR AudFreqWR
   dt <- dat[i,57:59]
   for(j in length(dt)){ row <- cbind(row,mean(dt[is.na(dt)==FALSE])) } # Widespread: freq
  # conformity "ConformityWR"
  # OrthopraxyWR OrthodoxyWR
   dt <- dat[i,60:61]
   for(j in length(dt)){ row <- cbind(row,mean(dt[is.na(dt)==FALSE])) } # conformity
  # Frequent: freq "FreqFR"
  # SpecFreqFR PartFreqFR AudFreqFR
   dt <- dat[i,62:64]
   for(j in length(dt)){ row <- cbind(row,mean(dt[is.na(dt)==FALSE])) } # Frequent: freq
  # conformity "ConformityFR"
  # OrthopraxyFR OrthodoxyFR
   dt <- dat[i,65:66]
   for(j in length(dt)){ row <- cbind(row,mean(dt[is.na(dt)==FALSE])) } # conformity
  # Euphoric: freq "FreqER"
  # SpecFreqER PartFreqER AudFreqER
   dt <- dat[i,67:69]
   for(j in length(dt)){ row <- cbind(row,mean(dt[is.na(dt)==FALSE])) } # Euphoric: freq
  # conformity "ConformityER"
  # OrthopraxyER OrthodoxyER
   dt <- dat[i,70:71]
   for(j in length(dt)){ row <- cbind(row,mean(dt[is.na(dt)==FALSE])) } # conformity
  # Dysphoric: freq "FreqDR"
  # SpecFreqDR PartFreqDR AudFreqDR
   dt <- dat[i,72:74]
   for(j in length(dt)){ row <- cbind(row,mean(dt[is.na(dt)==FALSE])) } # Dysphoric: freq
  # conformity "ConformityDR"
  # OrthopraxyDR OrthodoxyDR
   dt <- dat[i,75:76]
   for(j in length(dt)){ row <- cbind(row,mean(dt[is.na(dt)==FALSE])) } # conformity
  # Big Gods: MHG "MHG"
  # HighGods
   dt <- dat[i,81]
   for(j in length(dt)){ row <- cbind(row,mean(dt[is.na(dt)==FALSE])) } # Big Gods: MHG
  # Watts et al. developed and coded a new variable they call “Broad Supernatural Punishment” which may more closely match this relaxed definition of “Big Gods” than does the traditional “Moralizing High Gods” variable. Watts et al. define “Broad Supernatural Punishment” as follows: “For BSP to be coded as present in a culture there must be the concept of a supernatural agent or process that reliably monitors and punishes selfish actions, and this concept must (i) be widely advocated within the community, (ii) involve punishment of a broad range of selfish behaviours and (iii) apply to a wide range of community members.”
  # Because “selfish actions” can occur in a variety of domains, Seshat subdivides the types of supernatural enforcement of morality possible based on nine proposed categories of morality. For this study, we focused on three domains that are relevant to the establishment of large-scale cooperation:
  # 1) fairness (sharing of resources; e.g. dividing disputed resources, bargaining, redistribution of wealth),
  # 2) reciprocity (e.g. fulfilling contracts, returning gifts, repaying debts, upholding trust), and
  # 3) in-group loyalty (the need to remain loyal to unrelated members of the same group; e.g. helping coreligionists, going to war for one's group).
  # BSP was coded as present if any one of these three sub-types of selfish actions was supernaturally enforced.
  # Big Gods: BSP "BSP"
  # SEFair SERecip SEIngroup
   dt <- dat[i,c(83:84,87)]
   for(j in length(dt)){ row <- cbind(row,mean(dt[is.na(dt)==FALSE])) } # Big Gods: BSP
  # written records "WrittenRecords"
  # WRecord
   dt <- dat[i,33]
   for(j in length(dt)){ row <- cbind(row,mean(dt[is.na(dt)==FALSE])) } # written records
  # religious hierarchy "ReligiousHier"
  # ReligLev
   dt <- dat[i,6]
   for(j in length(dt)){ row <- cbind(row,mean(dt[is.na(dt)==FALSE])) } # religious hierarchy
   rowNA <- row
  # replace missing values with -99999
   rowNA[is.na(row)] <- -99999
   if(length(AggrDat)==0){AggrDat <- rbind(AggrDat,rowNA)}
   if(length(AggrDat)>0)
   {if(all(AggrDat[length(AggrDat[,1]),]==rowNA)==FALSE){AggrDat <- rbind(AggrDat,rowNA)} } 
}
# replace -99999 with NA
for(i in 1:length(AggrDat[,1])){
   for(j in 1:length(AggrDat[1,])){
      if(AggrDat[i,j]==-99999){AggrDat[i,j] <- NA}
   }}
# recombine aggregated data with NGA, PolID, Date and PropCoded
AggrDat <- cbind( data[row.names(AggrDat),c(1,2,3,4)] ,AggrDat)
# convert NGA from factor to character
AggrDat[,1] <- as.character(AggrDat[,1])
# convert PolID from factor to character
AggrDat[,2] <- as.character(AggrDat[,2])
# rename columns (why is Date renamed to Time?)
colnames(AggrDat) <- c("NGA", "PolID","Time", "PropCoded", "PolPop","PolTerr","CapPop","levels", "government", "infrastr", "writing", "texts", "money",   "FreqLR", "ConformityLR", "FreqWR",  "ConformityWR", "FreqFR", "ConformityFR", "FreqER", "ConformityER", "FreqDR", "ConformityDR", "MHG", "BSP", "WrittenRecords","ReligiousHier")
row.names(AggrDat) <- NULL

# Add General Moralistic Punishment variable
AggrDat$GeneralMoralisticPunishment<-ifelse(AggrDat$BSP>0.1,1,0)
# Add Moralising High Gods variable
# For consistency with previous studies that have generally used the “Moralizing High Gods” variable from the Ethnographic Atlas, the presence of MHG was coded as a binary variable based on this original definition:
# "As outlined by Murdock33 (1967:52), a high god follows the definition of Guy Swanson34 (1960: chapter III and appendix 1) as "a spiritual being who is believed to have created all reality and/or to be its ultimate governor, even though his sole act was to create other spirits who, in turn, created or control the natural world"… (1) "Absent or not reported," (2) "Present but not active in human affairs," (3) "Present and active in human affairs but not supportive of human morality" and (4) "Present, active, and specifically supportive of human morality"
# Thus, a coding of high gods “present, active, and specifically supportive of human morality” was coded as MHG being present, while all other types were coded as absent. 
AggrDat$MoralisingHighGods<-ifelse(AggrDat$MHG==3,1,0)
# Add Moralising Gods variable
# if either General Moralistic Punishment or Moralising High Gods are present, Moralising Gods are considered to be present
AggrDat$MoralisingGods<-ifelse(AggrDat$GeneralMoralisticPunishment>0.1 | AggrDat$MoralisingHighGods>0.1 ,1,0)
# Add Doctrinal Mode Variable
# The “Modes of Religiosity” hypothesis focuses on two factors that facilitate standardization of a body of beliefs and practices. First, high frequency (e.g. daily or weekly) collective rituals facilitate easy detection of deviations from the orthodox canon. Second, religious hierarchy enables enforcement of authorized belief and practice. Seshat codes five different types of rituals: the most frequent, most widespread, largest scale, most euphoric, and most dysphoric rituals. For each ritual, frequency is coded as daily, weekly, monthly, seasonally, yearly, generationally, or once-in-a-lifetime. Seshat also codes for levels of religious hierarchy. One represents no levels of religious hierarchy beyond the local priest/shaman, while higher numbers represent multiple levels of hierarchy (e.g. senior priests, High Druids).
# Making inferences about prehistoric rituals requires using various measurable archaeological proxies. Previous research has established that both frequent rituals and multi-level religious hierarchies tend to co-occur with other features of doctrinal rituals (e.g., low arousal). Not all of these features can always be found in the archaeological record, so in this paper we use the appearance of either religious hierarchy or frequent rituals as proxies for the appearance of doctrinal rituals. Doctrinal rituals were thus coded as present if the most frequent ritual occurred weekly or daily, or if there was evidence of multiple levels of religious hierarchy.
AggrDat$DoctrinalMode<-ifelse(AggrDat$FreqLR>4.5 | AggrDat$FreqWR>4.5 | AggrDat$FreqFR>4.5 | AggrDat$FreqER>4.5 | AggrDat$FreqDR>4.5 | AggrDat$ReligiousHier>=2,1,0)
# Add writing variable (why have this when there is the Written Records variable - just to convert into binary?)
AggrDat$Writing<-ifelse(AggrDat$WrittenRecords>0.1,1,0)

write.csv(AggrDat, file="./temp/MIAggrDat.csv",  row.names=FALSE)

rm(dat,data,dt,row,rowNA,i,j,money)

