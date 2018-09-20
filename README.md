# moralizing-gods

#The following code was used for analyses in Whitehouse, Francois, Savage, et al., "Complex societies precede moralizing gods throughout world history", with an "exportdat.csv" file scraped from the Seshat database on 10 Jan 2018.

######
CC By-NC SA License

Copyright (c) 2018 Peter Turchin and Patrick E. Savage

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software under the conditions of Creative Commons Attribution Non-Commercial (CC By-NC SA) licensing (https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode), subject to the following conditions:

Please include the following text in any publication using these data:

This research employed data from the Seshat Databank (seshatdatabank.info) under Creative Commons Attribution Non-Commercial (CC By-NC SA) licensing (https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode).

and cite:

1) Whitehouse, Francois, Savage, et al., "Complex societies precede moralizing gods throughout world history"

2) Turchin, P. et al. Quantitative historical analysis uncovers a single dimension of complexity that structures global variation in human social organization. Proc. Natl. Acad. Sci. U. S. A. 115, E144-E151 (2018).

3) Turchin, P. Fitting dynamical regression models to Seshat data. Cliodynamics 9(1):25-58 (2018).

4) Turchin P. et al. 2015. Seshat: The Global History Databank. Cliodynamics 6(1):77-107. 

The views and conclusions contained in this document are those of the authors and should not be interpreted as representing the official positions, either expressed or implied, of the Seshat Databank, its collaborative scholarly community, or the Evolution Institute.

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
######

#To run the primary analyses:
1) Download all files from https://github.com/pesavage/moralizing-gods
2) Download exportdat.csv (sent privately, or downloaded as Supplementary Table 1 [in which case it may need to be renamed as "exportdat.csv"])
3) Open "!MoralizingGods.R"
4) Modify the following code to set the working directory to the appropriate directory where you have saved these files:

setwd("/Users/pesavage/Documents/Research/Papers/Unpublished/Whitehouse Francois Savage et al Moralizing Gods/Nature resubmission/Reanalysis R code")

5) Copy and paste all code in "!MoralizingGods.R" into R or R Studio, and then you're done (takes ~30 minutes, mainly to perform pre-check and run 20 imputations of the data
6) Once you've run it once, you can save much time by starting at the following point in the !MI.R code:

######### end of the new scrape section
(if you do need to redo the multiple imputation process from the beginning, make to sure to replace the "polities.csv" file with the original one, as during the analysis process the polities.csv file is changed in ways that won't allow it to be used again from the beginning)

#To perform confirmatory analyses described in the Methods:
1) Removing hierarchy: 
Re-run from "# end of the new scrape section" with following changes:

-Replace "5:13" with "c(5:7,9:13)" in !MoralizingGods.R line 138 ("ImpDat <- ImpDatRepl[ImpDatRepl$irep==irep,5:13]")

-Replace "9" with "8" in !MoralizingGods.R line 134 ("Rotations <- matrix(NA,0,9)")

-Replace "9" with "8" in !MoralizingGods.R line 135 ("PropVar <- matrix(NA,0,9)")

2a) Scale variables only:
Re-run from "# end of the new scrape section" with following changes:

-Replace "5:13" with "5:8" in !MoralizingGods.R line 138 ("ImpDat <- ImpDatRepl[ImpDatRepl$irep==irep,5:13]")

-Replace "9" with "4" in !MoralizingGods.R line 134 ("Rotations <- matrix(NA,0,9)")

-Replace "9" with "4" in !MoralizingGods.R line 135 ("PropVar <- matrix(NA,0,9)")

2b) Non-scale variables only:
Re-run from "# end of the new scrape section" with following changes:

-Replace "5:13" with "9:13" in !MoralizingGods.R line 138 ("ImpDat <- ImpDatRepl[ImpDatRepl$irep==irep,5:13]")

-Replace "9" with "5" in !MoralizingGods.R line 134 ("Rotations <- matrix(NA,0,9)")

-Replace "9" with "5" in !MoralizingGods.R line 135 ("PropVar <- matrix(NA,0,9)")

3) Moralizing High gods:
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

#To perform additional analyses described in the Methods focusing on role of doctrinal rituals rather than moralizing gods:

1) Doctrinal ritual defined via religious hierarchy or ritual frequency:
-Re-run from "# end of the new scrape section" after replacing all instances in sub-scripts of "MoralisingGods" with "DoctrinalMode"

-Then, using this new "BigGodAnalysesEditedV2.R" file, re-run from "# end of the new scrape section" again with following changes:

2) Doctrinal ritual defined via religious hierarchy only:
-Add the following code to BigGodAnalysesEditedV2.R line 13:
dat$DoctrinalMode<-ifelse(dat$ReligiousHier>=2,1,0)
-and the following code to RegrDat.R line 18:
data$DoctrinalMode<-ifelse(dat$ReligiousHier>=2,1,0) 


3) Doctrinal ritual defined via ritual frequency only:
-Add the following code to BigGodAnalysesEditedV2.R line 13:
dat$DoctrinalMode<-ifelse(dat$FreqLR>4.5 | dat$FreqWR>4.5 | dat$FreqFR>4.5 | dat$FreqER>4.5 | dat$FreqDR>4.5,1,0)
-and the following code to RegrDat.R line 18:
data$DoctrinalMode<-ifelse(dat$FreqLR>4.5 | dat$FreqWR>4.5 | dat$FreqFR>4.5 | dat$FreqER>4.5 | dat$FreqDR>4.5,1,0) 
