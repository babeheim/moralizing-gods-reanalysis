library(maps)
library(plotrix)
library(dplyr)
library(plyr) # manually call plyr::rename to avoid problems

library(testthat)
library(viridis)
library(rethinking) # github.com/rmcelreath/rethinking

set.seed(1234)

# number of imputations in `2_impute_data`
nrep <- 5
# this was 20 in the publication version

Section1 <- "Social Complexity variables"
Section2 <- "Ritual variables"
Section3 <- "Religion and Normative Ideology"

dir_init <- function(path, verbose=FALSE, overwrite=TRUE){
  if(substr(path, 1, 2)!='./') stop('path argument must be formatted
    with "./" at beginning')
  contents <- dir(path, recursive=TRUE)
  if(dir.exists(path)){
    if(overwrite){
      if(verbose){
        if(length(contents)==0) print(paste('folder ', path, ' created.', sep=""))
        if(length(contents)>0) print(paste('folder ', path, ' wiped of ', length(contents), ' files/folders.', sep=""))
      }
      if(dir.exists(path)) unlink(path, recursive=TRUE)
      dir.create(path)
    }
  } else {
    if(verbose){
      print(paste('folder ', path, ' created.', sep=""))
    }
    dir.create(path)
  }
}

col_alpha <- function (acol, alpha = 0.2){
  acol <- col2rgb(acol)
  acol.red <- acol["red",]/255
  acol.green <- acol["green",]/255
  acol.blue <- acol["blue",]/255
  acol <- mapply(function(red, green, blue, alphas) rgb(red, green, blue, alphas), acol.red, acol.green, acol.blue, alpha)
  return(as.character(acol))
}