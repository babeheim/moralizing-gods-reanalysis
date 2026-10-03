
rm(list = ls())


restore_environment <- function(project = ".") {

  project <- normalizePath(
    project,
    winslash = "/",
    mustWork = TRUE
  )

  lockfile <- file.path(
    project,
    "renv.lock"
  )

  activate_file <- file.path(
    project,
    "renv",
    "activate.R"
  )

  if (!file.exists(lockfile)) {
    stop(
      "renv.lock was not found: ",
      lockfile
    )
  }

  if (!file.exists(activate_file)) {
    stop(
      "renv activation script was not found: ",
      activate_file
    )
  }

  old_wd <- getwd()

  on.exit(
    setwd(old_wd),
    add = TRUE
  )

  setwd(project)

  # Bootstrap and activate the project's own renv environment.
  #
  # renv/activate.R can bootstrap renv itself on a fresh machine,
  # so renv does not need to be installed globally beforehand.
  source(
    activate_file,
    local = .GlobalEnv
  )

  if (!requireNamespace(
    "renv",
    quietly = TRUE
  )) {
    stop(
      "renv could not be bootstrapped from renv/activate.R."
    )
  }

  # Restore exactly the package versions recorded in renv.lock.
  renv::restore(
    project = project,
    lockfile = lockfile,
    prompt = FALSE,
    retry = FALSE
  )

  # Explicitly ensure that this running R process is using the
  # restored project library.
  renv::load(
    project = project,
    quiet = TRUE
  )

  invisible(TRUE)
}


project_root <- normalizePath(
  getwd(),
  winslash = "/",
  mustWork = TRUE
)

restore_environment(project_root)






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

# original analysis above; reanalysis begins here

dir_init("./09_explore_data/input")
file.copy("./06_prep_regression_data/output/RegrDat.csv",
  "./09_explore_data/input/")
setwd("./09_explore_data")
source("explore_data.r")
setwd("..")

dir_init("./10_test_forward_bias/input")
files <- "./02_impute_data/output/polities.csv"
files <- c(files, "./03_run_pca/output/PC1_traj_merged.csv")
file.copy(files, "./10_test_forward_bias/input/")
setwd("./10_test_forward_bias")
source("test_forward_bias.r")
setwd("..")

dir_init("./11_fit_hierarchical_models/input")
files <- "./02_impute_data/output/polities.csv"
files <- c(files, "./03_run_pca/output/PC1_traj_merged.csv")
file.copy(files, "./11_fit_hierarchical_models/input/")
setwd("./11_fit_hierarchical_models")
source("fit_hierarchical_models.r")
setwd("..")

dir_init("./12_fit_revised_binomials/input")
file.copy("./06_prep_regression_data/output/RegrDat.csv",
  "./12_fit_revised_binomials/input/")
setwd("./12_fit_revised_binomials")
source("fit_revised_binomials.r")
setwd("..")

dir_init("./13_explore_models/input")
files <- c("./04_prep_comparisons/output/PrePostComparison.csv",
  "./05_draw_figures/output/SCNorm.csv",
  "./06_prep_regression_data/output/RegrDat.csv")
files <- c(files, list.files("./12_fit_revised_binomials/output", full.names = TRUE))
file.copy(files, "./13_explore_models/input/")
setwd("./13_explore_models")
source("explore_models.r")
setwd("..")

dir_init("./output")
files <- list.files("./09_explore_data/output", pattern = ".png", full.names = TRUE)
files <- c(files, list.files("./10_test_forward_bias/output", pattern = ".png", full.names = TRUE))
files <- c(files, list.files("./11_fit_hierarchical_models/output", pattern = ".png", full.names = TRUE))
files <- c(files, list.files("./13_explore_models/output", pattern = ".png", full.names = TRUE))
file.copy(files, "./output")
