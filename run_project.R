
rm(list = ls())

source("./project_support.R")

dir_init("./1_precheck/input/")
files <- c("./data/polities.csv", "./data/variables.csv", "./data/exportdat.csv")
file.copy(files, "./1_precheck/input/")
setwd("./1_precheck")
source("precheck.R")
setwd("..")

dir_init("./2_impute_data/input/")
files <- c("./data/polities.csv", "./data/variables.csv")
files <- c(files, "./1_precheck/output/SCdat.csv")
file.copy(files, "./2_impute_data/input/")
setwd("./2_impute_data")
source("./impute_data.R")
setwd("..")
# takes a few mins

dir_init("./3_run_pca/input")
files <- c("./data/NGAcoords.csv")
files <- c(files, "./2_impute_data/output/polities.csv",
  "./2_impute_data/output/MIoutput.csv", "./2_impute_data/output/ImpDatRepl.csv")
file.copy(files, "./3_run_pca/input/")
setwd("./3_run_pca")
source("./run_pca.R")
setwd("..")

dir_init("./4_run_tests/input")
files <- c("./2_impute_data/output/polities.csv")
files <- c(files, "./3_run_pca/output/PC1_traj_merged.csv")
file.copy(files, "./4_run_tests/input/")
setwd("./4_run_tests")
source("BigGodAnalysesEditedV2.R")
setwd("..")

dir_init("./5_prep_regrdat/input")
files <- c("./data/DistMatrix.csv")
files <- c(files, "./2_impute_data/output/polities.csv")
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
