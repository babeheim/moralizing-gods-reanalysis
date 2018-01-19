# moralizing-gods-test
Code used for analyses in Whitehouse, Francois, Savage, et al., "Complex societies and doctrinal rituals precede moralizing gods throughout world history"

Uploading code to share with collaborators, then reviewers, then readers

#To run the primary analyses:
1) Download all files
2) Download exportdat.csv (sent privately)
3) Open "!MI.R"
4) Modify the following code to set the working directory to the appropriate directory where you have saved these files:
setwd("/Users/pesavage/Documents/Research/Oxford Seshat/Data/SCBigGodsOct2017")
5) Copy and paste all code in "!MI.R" into R or R Studio, and then you're done (takes ~30 minutes, mainly to perform pre-check and run 20 imputations of the data
6) Once you've run it once, you can save much time by starting at the following point in the !MI.R code:
######### end of the new scrape section
(if you do need to redo the multiple imputation process from the beginning, make to sure to replace the "polities.csv" file with the original one, as during the analysis process the polities.csv file is changed in ways that won't allow it to be used again from the beginning)

#To perform confirmatory analyses described in the Methods:
1) Removing hierarchy: 
Re-run from "# end of the new scrape section" with following changes:
-Replace "5:13" with "c(5:7,9:13)" in !MI.R line 105 ("ImpDat <- ImpDatRepl[ImpDatRepl$irep==irep,5:13]")
-Replace "9" with "8" in !MI.R line 101 ("Rotations <- matrix(NA,0,9)")
-Replace "9" with "8" in !MI.R line 102 ("PropVar <- matrix(NA,0,9)")

2a) Scale variables only:
Re-run from "# end of the new scrape section" with following changes:
-Replace "5:13" with "5:8" in !MI.R line 105 ("ImpDat <- ImpDatRepl[ImpDatRepl$irep==irep,5:13]")
-Replace "9" with "4" in !MI.R line 101 ("Rotations <- matrix(NA,0,9)")
-Replace "9" with "4" in !MI.R line 102 ("PropVar <- matrix(NA,0,9)")

2b) Non-scale variables only:
Re-run from "# end of the new scrape section" with following changes:
-Replace "5:13" with "9:13" in !MI.R line 105 ("ImpDat <- ImpDatRepl[ImpDatRepl$irep==irep,5:13]")
-Replace "9" with "5" in !MI.R line 101 ("Rotations <- matrix(NA,0,9)")
-Replace "9" with "5" in !MI.R line 102 ("PropVar <- matrix(NA,0,9)")

3) Big gods:
Re-run from "# end of the new scrape section" after replacing all instances in sub-scripts of "MoralisingGods" with "MoralisingHighGods". 

4) Dating uncertainty:
add following line between "MG<-..." and "MGAppear<-..." lines in BigGodAnalysesEditedV2.R: 
MG<-as.data.frame(MG %>% group_by(PolID) %>% sample_n(size = 1)) 
(NB: add "library(dplyr)" before the loop, but note that R bugs when you load both dplyr and plyr)

5a) 700-year max time-window:
Re-run from "# end of the new scrape section" with following changes:
-Replace "2050" with "750" in BigGodAnalysesEditedV2.R line 67 ("out <-subset(out, out[,5]<2050)")

5b) All possible time-windows:
Re-run from "# end of the new scrape section" after deleting/hashing out BigGodAnalysesEditedV2.R line 67 ("out <-subset(out, out[,5]<2050)")
