library(maps)
library(plotrix)
library(plyr) # this has to be first, also manually call plyr::rename to avoid problems
library(dplyr)

library(testthat)
library(viridis)
library(rethinking) # github.com/rmcelreath/rethinking

library(dplyr)
library(glmmTMB)
library(glmmADMB)
library(lme4)
library(DHARMa)
library(bbmle)
library(reshape)
library(yarrr)

library(rmarkdown)

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

density_offset <- function(y, scale = 1) {
  dens <- density(y)
  y_dens <- sapply(y, function(z) dens$y[which.min(abs(z - dens$x))])
  y_dens <- y_dens / max(y_dens)
  y_offset <- rnorm(length(y), 0, sd = scale * sqrt(y_dens))
  return(y_offset)
}

texttab <- function(dataframe, hlines = NA){
  inmat <- as.matrix(dataframe)
  inmat <- rbind(colnames(inmat), inmat)
  if(!is.null(rownames(inmat))) inmat <- cbind(rownames(inmat), inmat)
  output <- character(nrow(inmat))
  for(i in 1:nrow(inmat)){
    add.amps <- paste(inmat[i,], collapse=" & ")
    output[i] <- paste(add.amps, "\\\\", sep=" ")
  }
  if(all(!is.na(hlines))){
    for(i in 1:length(hlines)) output <- append(output, "\\hline", hlines[i]+(i-1))
  }
  return(output)
}

psign <- function(samples){
  if(mean(samples)>0) output <- round(mean(samples < 0), 3)
  if(mean(samples)<0) output <- round(mean(samples > 0), 3)
  output <- sprintf("%.2f", output)
  output[output=="0.00"] <- "$<$0.01"
  return(output)
}
