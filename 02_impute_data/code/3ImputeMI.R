# Single Imputation function using stochastic regression

# For missing data, we impute values as follows. Suppose that, for some polity, we have a missing value for variable A and coded values for variables B–H. We select a subset of cases from the full dataset, in which all values of A–H variables have values and build a regression model for A. Not all predictors B–H may be relevant to predicting A, and thus, the first step is selecting which of the predictors should enter the model (information on model selection is given below). After the optimal model is identified, we estimate its parameters. Then, we go back to the polity (where variable A is missing) and use the known values of predictor variables for this polity to calculate the expected value of A using the estimated regression coefficients. However, we do not simply substitute the missing value with the expected one (because as explained above, this is known to result in biased estimates). Instead, we sample from the posterior distribution characterizing the prediction of the regression model (in practice, we randomly sample the regression residual and add it to the expected value). We applied the same approach to each missing value in the dataset, yielding an imputed dataset without gaps. The overall imputation procedure was repeated 20 times, yielding 20 imputed sets that were used in the analyses.
AggrDat <- read.table('./temp/MIAggrDat.csv', sep=",", header=TRUE)
# subset PolPop, PolTerr, CapPop, levels, government, infrastr, writing, texts, money
dat <- AggrDat[,5:13]
# log base 10 PolPop, PolTerr and CapPop
dat[,1:3] <- log10(dat[,1:3])
ImpDat <- dat

# 1:ncol(dat)
for(j in 1:length(dat[1,])){ 
  # 1:nrow(dat)
  for(i in 1:length(dat[,1])){ 
    # if cell is NA
    if(is.na(dat[i,j])==TRUE){
      # create index of columns in dat (1:ncol(dat) or 1:9)
      index <- c(1:length(dat[1,]))
      # index which columns have NA values per row and extract an index of only columns with non-NA values
      index <- index[is.na(dat[i,])==FALSE] 
      # extract variables with column of interest and columns without missing values
      RegrDat <- dat[,c(j,index)]
      # remove rows with missing values in column of interest
      RegrDat <- RegrDat[is.na(RegrDat[,1])==FALSE,] 
      # fit linear model
      fit <- lm(RegrDat)
      # extract probability Pr(>|t|)
      Pval <- summary(fit)$coefficients[,4]  
      # remove model intercept
      Pval <- Pval[-1]
      # if all p values are not less than 0.05 then set the smallest p value to 0.01
      if (all((Pval < 0.05)==FALSE)){Pval[Pval==min(Pval)] <- 0.01}
      # index all p values < 0.05
      index <- index[Pval < 0.05]  
      # extract column of interest and all columns with p < 0.05
      RegrDat <- dat[,c(j,index)] 
      # remove rows in column of interest with missing values
      RegrDat <- RegrDat[is.na(RegrDat[,1])==FALSE,]
      # fit linear model 
      fit <- lm(RegrDat)      
      # extract predictors from data
      predictors <- dat[i,c(j,index)]  
      # replace missing value with 1
      predictors[1] <- 1
      # extract model coefficients
      coeff <- coefficients(fit) 
      # randomly sample from model residuals and add the sum of coefficients * predictors
      ImpDat[i,j] <- sum(coeff*predictors) + sample(fit$residuals,1) 
    }
  }
}

rm(dat,predictors,RegrDat,coeff,fit,i,index,j,Pval)
