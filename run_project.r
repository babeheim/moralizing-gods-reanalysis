
rm(list = ls())

source("./project_support.r")

dir_init("./01_precheck/input/")
files <- c("./data/polities.csv", "./data/variables.csv",
  "./data/exportdat.csv")
file.copy(files, "./01_precheck/input/")
setwd("./01_precheck")
source("precheck.r")
setwd("..")

dir_init("./02_impute_data/input/")
files <- c("./data/polities.csv", "./data/variables.csv")
files <- c(files, "./01_precheck/output/SCdat.csv")
file.copy(files, "./02_impute_data/input/")
setwd("./02_impute_data")
source("./impute_data.r")
setwd("..")

dir_init("./03_run_pca/input")
files <- c("./data/NGAcoords.csv")
files <- c(files, "./02_impute_data/output/polities.csv",
  "./02_impute_data/output/MIoutput.csv",
  "./02_impute_data/output/ImpDatRepl.csv")
file.copy(files, "./03_run_pca/input/")
setwd("./03_run_pca")
source("./run_pca.r")
setwd("..")

dir_init("./04_prep_comparisons/input")
files <- c("./02_impute_data/output/polities.csv")
files <- c(files, "./03_run_pca/output/PC1_traj_merged.csv")
file.copy(files, "./04_prep_comparisons/input/")
setwd("./04_prep_comparisons")
source("prep_comparisons.r")
setwd("..")

dir_init("./05_draw_figures/input")
files <- c("./03_run_pca/output/PC1_traj_merged.csv")
files <- c(files, "./04_prep_comparisons/output/PrePostComparison.csv")
file.copy(files, "./05_draw_figures/input/")
setwd("./05_draw_figures")
source("draw_figures.r")
setwd("..")

dir_init("./06_prep_regression_data/input")
files <- c("./data/DistMatrix.csv")
files <- c(files, "./02_impute_data/output/polities.csv")
files <- c(files, "./03_run_pca/output/PC1_traj_merged.csv")
file.copy(files, "./06_prep_regression_data/input/")
setwd("./06_prep_regression_data")
source("prep_regression_data.r")
setwd("..")

dir_init("./07_fit_regressions/input")
file.copy("./06_prep_regression_data/output/RegrDat.csv",
  "./07_fit_regressions/input/")
setwd("./07_fit_regressions")
source("fit_regressions.r")
setwd("..")

dir_init("./08_create_map/input")
file.copy("./data/map.csv", "./08_create_map/input/")
setwd("./08_create_map")
source("create_map.r")
setwd("..")

# publication analysis above; reanalysis begins here

dir_init("./09_draw_missingness/input")
file.copy("./06_prep_regression_data/output/RegrDat.csv",
  "./09_draw_missingness/input/")
setwd("./09_draw_missingness")
source("draw_missingness.r")
setwd("..")

dir_init("./10_fit_new_regressions/input")
file.copy("./06_prep_regression_data/output/RegrDat.csv",
  "./10_fit_new_regressions/input/")
setwd("./10_fit_new_regressions")
source("fit_new_regressions.r")
setwd("..")

dir_init("./11_draw_new_figures/input")
file.copy("./04_prep_comparisons/output/PrePostComparison.csv",
  "./11_draw_new_figures/input/")
file.copy("./05_draw_figures/output/SCNorm.csv",
  "./11_draw_new_figures/input/")
file.copy("./06_prep_regression_data/output/RegrDat.csv",
  "./11_draw_new_figures/input/")
files <- list.files("./10_fit_new_regressions/output", full.names = TRUE)
file.copy(files, "./11_draw_new_figures/input/")
setwd("./11_draw_new_figures")
source("draw_new_figures.r")
setwd("..")

source("./project_support.r")
files <- c("./02_impute_data/output/polities.csv")
files <- c(files, "./03_run_pca/output/PC1_traj_merged.csv")
file.copy(files, "./12_causal_analysis/", overwrite = TRUE)
setwd("./12_causal_analysis")
render("Reanalysis_of_BigGodAnalyses.Rmd")
setwd("..")

source("./project_support.r")
files <- c("./02_impute_data/output/ImpDatRepl.csv")
file.copy(files, "./13_MG_writing_analysis/", overwrite = TRUE)
setwd("./13_MG_writing_analysis")
source("MG_writing_analysis.r")
setwd("..")

dir_init("./output")
files <- list.files("./09_draw_missingness/output", pattern = ".", full.names = TRUE)
files <- c(files, list.files("./11_draw_new_figures/output", pattern = ".", full.names = TRUE))
file.copy(files, "./output")
