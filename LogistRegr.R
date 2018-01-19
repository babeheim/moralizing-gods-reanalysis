##### Logistic Regression 
## dat = 1st column response variable (binary), subsequent columns predictors
LogistRegrDat <- LogistRegrDat[is.na(LogistRegrDat[,1])==FALSE,]
Y <- LogistRegrDat[,1]
X <- LogistRegrDat[,2:length(LogistRegrDat[1,])]
if(length(LogistRegrDat[1,])==2){reslt <- glm(Y ~ X, family=binomial(link='logit'))}
if(length(LogistRegrDat[1,])>2){reslt <- glm(Y ~ ., data = X, family=binomial(link='logit'))}

print(summary(reslt))

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
