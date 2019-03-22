# The following code was used for analyses in Whitehouse, Francois, Savage, et al.,
# "Complex societies and doctrinal rituals precede moralizing gods throughout world history",
# with an "exportdat.csv" file scraped from the Seshat database on 19 Jan 2018.

# This software was developed by Peter Turchin and Patrick Savage (Copyright 29 Jan 2018).

# For the general methodology underlying these analyses please refer to: 
# Turchin, P. et al. Quantitative historical analysis uncovers a single dimension of complexity
#   that structures global variation in human social organization. Proc. Natl. Acad. Sci. U. S. A. 115, E144-E151 (2018).
# Turchin, P. Fitting dynamical regression models to Seshat data. Cliodynamics 9, (2018).  

######

# Make sure to use a full polities.csv file (if you re-run using the output the multiple
# imputation loop will run into bugs when it comes across NGAs with only one polity)

### NB: Using shortened variables.csv file that only uses variables of interest to save computational time.
# If using full file, may need to renumber variables in Aggr.MI

rm(list = ls())

source("./project_support.R")

dir_init("./1_precheck/input/")
files <- c("./data/polities.csv", "./data/variables.csv", "./data/exportdat.csv")
file.copy(files, "./1_precheck/input/")
setwd("./1_precheck")
source("precheck.R")
setwd("..")

dir_init("./2_impute_data/input")
files <- c("./data/polities.csv", "./data/variables.csv")
files <- c(files, "./1_precheck/output/SCdat.csv")
file.copy(files, "./2_impute_data/input/")
setwd("./2_impute_data")
source("./impute_data.R")
setwd("..")
# takes a few mins

dir_init("./3_run_pca/input")
files <- c("./data/polities.csv", "./data/NGAcoords.csv")
files <- c(files, "./2_impute_data/output/MIoutput.csv", "./2_impute_data/output/ImpDatRepl.csv")
file.copy(files, "./3_run_pca/input/")
setwd("./3_run_pca")
source("./run_pca.R")
setwd("..")

dir_init("./4_run_tests/input")
files <- c("./data/polities.csv")
files <- c(files, "./3_run_pca/output/PC1_traj_merged.csv")
file.copy(files, "./4_run_tests/input/")
setwd("./4_run_tests")
source("BigGodAnalysesEditedV2.R")
setwd("..")

dir_init("./5_prep_regrdat/input")
files <- c("./data/polities.csv", "./data/DistMatrix.csv")
files <- c(files, "./3_run_pca/output/PC1_traj_merged.csv")
file.copy(files, "./5_prep_regrdat/input/")
setwd("./5_prep_regrdat")
source("RegrDat.R")
setwd("..")

dir_init("./6_fit_regressions/input")
file.copy("./5_prep_regrdat/output/RegrDat.csv", "./6_fit_regressions/input/")
setwd("./6_fit_regressions")
source("LogistRegr.R")
setwd("..")

dir_init("./7_create_map/input")
file.copy("./data/map.csv", "./7_create_map/input/")
setwd("./7_create_map")
source("Seshat_BG_map.r")
setwd("..")
