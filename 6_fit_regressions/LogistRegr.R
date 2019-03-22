

rm(list = ls())
source("../project_support.R")

RegrDat <- read.csv("./input/RegrDat.csv", stringsAsFactors = FALSE)

RD1 <- RegrDat[is.na(RegrDat$Lag1) == FALSE,]
RD2 <- RD1[is.na(RD1$Lag2) == FALSE,]
LogistRegrDat = as.data.frame(cbind(RD2$MG,RD2[,c(36,51:54)]))

##### Logistic Regression 
## dat = 1st column response variable (binary), subsequent columns predictors
LogistRegrDat <- LogistRegrDat[is.na(LogistRegrDat[,1])==FALSE,]
Y <- LogistRegrDat[,1]
X <- LogistRegrDat[,2:length(LogistRegrDat[1,])]
if(length(LogistRegrDat[1,])==2){reslt <- glm(Y ~ X, family=binomial(link='logit'))}
if(length(LogistRegrDat[1,])>2){reslt <- glm(Y ~ ., data = X, family=binomial(link='logit'))}

summary(reslt)

# new analysis block

post <- extract.samples(reslt)

d <- LogistRegrDat

scs <- seq(0, 1, by = 0.01)

has_mg_mean <- rep(NA, length(scs))
has_mg_sd <- rep(NA, length(scs))
has_mg_lb <- rep(NA, length(scs))
has_mg_ub <- rep(NA, length(scs))

for (i in 1:length(scs)) {
  has_mg <- logistic(post$Intercept + post$Lag1 * 0 + post$Mean * scs[i])
  has_mg_mean[i] <- mean(has_mg)
  has_mg_sd[i] <- sd(has_mg)
  has_mg_lb[i] <- HPDI(has_mg)[1]
  has_mg_ub[i] <- HPDI(has_mg)[2]
}


has_mg <- logistic(post$Intercept + post$Lag1 * 0 + post$Mean * 0.1)
has_mg <- logistic(post$Intercept + post$Lag1 * 0 + post$Mean * 0.4)

# d$Pop <- RD2$CapPop
# d$NGA <- as.character(RD2$NGA)
# d$Time <- RD2$Time

d$pr_mg_mean <- NA
d$pr_mg_sd <- NA
d$pr_mg_lb <- NA
d$pr_mg_ub <- NA

for (i in 1:nrow(d)) {
  has_mg <- logistic(
    post$Intercept +
    post$Lag1 * d$Lag1[i] +
    post$Lag2 * d$Lag2[i] +
    post$Mean * d$Mean[i] +
    post$Space * d$Space[i] +
    post$Phylogeny * d$Phylogeny[i]
  )
  d$pr_mg_mean[i] <- mean(has_mg)
  d$pr_mg_sd[i] <- sd(has_mg)
  d$pr_mg_lb[i] <- HPDI(has_mg)[1]
  d$pr_mg_ub[i] <- HPDI(has_mg)[2]
}

plot(scs, has_mg_mean, ylim = c(0, 1), type = "l", ylab = "pr(moralizing gods)", xlab = "social complexity")
polygon(c(scs, rev(scs)), c(has_mg_ub, rev(has_mg_lb)), border = NA, col = col.alpha("dodgerblue", 0.1))

tar <- which(d$Lag1 == 0)
points(d$Mean[tar], d$pr_mg_mean[tar], col = col.alpha("black", 0.8))

abline(v = 0.4, lty = 2)

####

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
