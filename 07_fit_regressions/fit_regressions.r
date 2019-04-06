

rm(list = ls())
source("../project_support.r")

dir_init("./temp")

#########

RegrDat <- read.csv("./input/RegrDat.csv", stringsAsFactors = FALSE)

RD1 <- RegrDat[is.na(RegrDat$Lag1) == FALSE,]
RD2 <- RD1[is.na(RD1$Lag2) == FALSE,]

# LogistRegrDat = as.data.frame(cbind(RD2$MG,RD2[,c(36,51:54)]))
## replaced this line with an explicit call
LogistRegrDat <- RD2[,c("MG", "Mean", "Lag1", "Lag2", "Phylogeny", "Space")]

##### Logistic Regression 
## dat = 1st column response variable (binary), subsequent columns predictors
LogistRegrDat <- LogistRegrDat[is.na(LogistRegrDat[,1])==FALSE,]

expect_equal(dim(LogistRegrDat), c(801, 6))
expect_equal(sum(is.na(LogistRegrDat[,1])), 0)

Y <- LogistRegrDat[,1]
X <- LogistRegrDat[,2:length(LogistRegrDat[1,])]
if(length(LogistRegrDat[1,])==2){reslt <- glm(Y ~ X, family=binomial(link='logit'))}
if(length(LogistRegrDat[1,])>2){reslt <- glm(Y ~ ., data = X, family=binomial(link='logit'))}

print(summary(reslt))

expect_true(abs(as.numeric(logLik(reslt)) - (-79)) < 1)
expect_true(abs(coef(reslt)[1] - (-7.9)) < 0.15)
expect_true(abs(coef(reslt)[2] - (9.9)) < 0.1)



ObsPred <- cbind(LogistRegrDat[,1], predict(reslt, type="response"))
colnames(ObsPred) <- c("Obs","Pred")

R_tjur <- mean(ObsPred[ObsPred[,1]==1,2]) - mean(ObsPred[ObsPred[,1]==0,2])
R_tjur <- 0.001*round(R_tjur*1000)
print(R_tjur)

par(mfrow=c(2,1))
hist(ObsPred[ObsPred[,1]==0,2], main=paste("Response = 0       Tjur's coeff. discr. = ",R_tjur), xlab="Predicted Value")
hist(ObsPred[ObsPred[,1]==1,2], main=paste("Response = 1       Tjur's coeff. discr. = ",R_tjur), xlab="Predicted Value")
par(mfrow=c(1,1))

rm(ObsPred, R_tjur, X, Y)
