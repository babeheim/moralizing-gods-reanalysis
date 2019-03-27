

rm(list = ls())
source("../project_support.r")

dir_init("./temp")

RegrDat <- read.csv("./input/RegrDat.csv", stringsAsFactors = FALSE)

RD1 <- RegrDat[is.na(RegrDat$Lag1) == FALSE,]
RD2 <- RD1[is.na(RD1$Lag2) == FALSE,]
# LogistRegrDat = as.data.frame(cbind(RD2$MG,RD2[,c(36,51:54)]))
## replaced this line with an explicit call
LogistRegrDat <- RD2[,c("NGA", "Time", "MG", "Mean", "Lag1", "Lag2", "Phylogeny", "Space", "MG_missing")]

##### Logistic Regression 
## dat = 1st column response variable (binary), subsequent columns predictors

# publication results here

d <- LogistRegrDat

d$key <- paste(d$NGA, d$Time)

appearances <- d$key[which(d$Lag1 == 0)] # preserve these cases for comparison

reslt <- glm(MG ~ Mean + Lag1 + Lag2 + Phylogeny + Space, data = d, family=binomial(link='logit'))

summary(reslt)

# make a plot to show its the same as what's about to come next

post <- extract.samples(reslt)

scs <- seq(0, 1, by = 0.01)

has_mg_mean <- rep(NA, length(scs))
has_mg_sd <- rep(NA, length(scs))
has_mg_lb <- rep(NA, length(scs))
has_mg_ub <- rep(NA, length(scs))

for (i in 1:length(scs)) {
  has_mg <- logistic(post$Intercept + post$Mean * scs[i])
  has_mg_mean[i] <- mean(has_mg)
  has_mg_sd[i] <- sd(has_mg)
  has_mg_lb[i] <- HPDI(has_mg)[1]
  has_mg_ub[i] <- HPDI(has_mg)[2]
}

has_mg <- logistic(post$Intercept + post$Mean * (-0.5))
mean(has_mg)
sd(has_mg)

has_mg <- logistic(post$Intercept + post$Lag1 * 0 + post$Mean * 0.42)
mean(has_mg)
sd(has_mg)

# calculation posterior distributions on p for each case
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

png("./temp/pred_comparison.png", res = 300, height = 5, width = 10, units = "in")

par(mfrow = c(1, 2))

plot(scs, has_mg_mean, ylim = c(0, 1), type = "l",
  ylab = "pr(moralizing gods)", xlab = "social complexity",
  main = "published analysis")

polygon(c(scs, rev(scs)), c(has_mg_ub, rev(has_mg_lb)),
  border = NA, col = col.alpha("dodgerblue", 0.2))

tar <- which(d$key %in% appearances)
print(paste(length(tar), "cases after lag1 = 0"))
points(d$Mean[tar], d$pr_mg_mean[tar], pch = 16,
  col = col_alpha("black", 0.5))

# for (i in 1:length(tar)) lines(c(d$Mean[tar[i]], d$Mean[tar[i]]),
#   c(d$pr_mg_lb[tar[i]], d$pr_mg_ub[tar[i]]))



## bounding box of interior figure
xmin <- 0.15
xmax <- 0.62
ymin <- 0.45
ymax <- 0.95

# prep variables on the new scale

d$logPop <- RD2$PolPop
d$Pop <- 10^d$logPop

d$Mean2 <- (d$Mean - 0) / 1 * (xmax - xmin) + xmin
d$logPop2 <- (d$logPop - min(d$logPop[tar])) / max(d$logPop[tar]) * (ymax - ymin) + ymin

d$na_col <- ifelse(d$MG_missing == 1, "gray", "black")

# plot

mg <- 0.02 # don't mix up with outcome variable `MG`!

tar <- which(d$key %in% appearances & d$MG_missing == 1)
points(d$Mean2[tar], d$logPop2[tar], pch = 16, col = col_alpha(d$na_col[tar], 0.8))
tar <- which(d$key %in% appearances & d$MG_missing == 0)
points(d$Mean2[tar], d$logPop2[tar], pch = 16, col = col_alpha(d$na_col[tar], 0.8))

lines(c(xmin, xmin) - mg, c(ymin, ymax))
lines(c(xmin, xmax), c(ymin, ymin) - mg)

xtck <- seq(0, 1, by = 0.25) * (xmax - xmin) + xmin
ytck <- rep(ymin, times = length(xtck))
for(i in 1:length(xtck)) lines(c(xtck[i], xtck[i]), c(ytck[i] - mg, ytck[i] - 2 * mg))

text(xtck, ytck - 4 * mg, labels = c("0", "", "0.5", "", "1"))

ytck <- (seq(1, 8, by = 1) - 1)/8 * (ymax - ymin) + ymin
xtck <- rep(xmin, times = length(ytck))
for(i in 1:length(ytck)) lines(c(xtck[i] - mg, xtck[i] - 2 * mg), c(ytck[i], ytck[i]))

text(xtck - 4 * mg, ytck, labels = c("", "2", "", "4", "", "6", "", "8"))

text(xmin + (xmax - xmin) / 2, ymin - 6 * mg,
  label = "social complexity", cex = 0.8)
text(xmin - 6 * mg, ymin + (ymax - ymin) / 2,
  label = "log10 population", srt = 90, cex = 0.8)

text(xmin + (xmax - xmin) / 3, ymax - 2 * mg,
  label = "observed", col = "black")
text(xmin + (xmax - xmin) / 3, ymax - 5 * mg,
  label = "missing", col = gray(0.4))

m0 <- lm(logPop2 ~ Mean2, data = d)

curve(coef(m0)[1] + coef(m0)[2] * x, from = xmin,
  to = xmax, lty = 2, add = TRUE)





## second analysis, this time dropping the 0's that were NA's

RegrDat <- read.csv("./input/RegrDat.csv", stringsAsFactors = FALSE)

# put the NA's back
RegrDat$MG[RegrDat$MG_missing == 1] <- NA

# kill lagged variables that were missing values

NGAs <- sort(unique(RegrDat$NGA))

for (i in 1:length(NGAs)) {
  nga_rows <- which(RegrDat$NGA == NGAs[i])
# if the MG_missing is 1, the Lag1 one time step forward must be NA, if exists
  if (length(nga_rows) > 1) {
    for (j in 1:(length(nga_rows) - 1)) {
      if (RegrDat$MG_missing[nga_rows[j]] == 1) RegrDat$Lag1[nga_rows[j + 1]] <- NA
    }
  }
# if the MG_missing is 1, the Lag2 two time steps forward must be NA, if exists
  if (length(nga_rows) > 2) {
    for (j in 1:(length(nga_rows) - 2)) {
      if (RegrDat$MG_missing[nga_rows[j]] == 1) RegrDat$Lag2[nga_rows[j + 2]] <- NA
    }
  }
}

# RegrDat[which(RegrDat$NGA == NGAs[i]),c("NGA", "Time", "MG", "Lag1", "Lag2", "MG_missing")]; i <- i + 1

## replaced this line with an explicit call
LogistRegrDat <- RegrDat[,c("NGA", "Time", "MG", "Mean", "Lag1", "Lag2", "Phylogeny", "Space", "MG_missing")]

d <- LogistRegrDat

d$key <- paste(d$NGA, d$Time)


mg_na <- which(is.na(d$MG))
dm <- d[-mg_na, ]

reslt <- glm(MG ~ Mean + Space, data = d, family=binomial(link='logit'))

summary(reslt)


# make a plot to show its the same as what's about to come next

post <- extract.samples(reslt)

scs <- seq(0, 1, by = 0.01)

has_mg_mean <- rep(NA, length(scs))
has_mg_sd <- rep(NA, length(scs))
has_mg_lb <- rep(NA, length(scs))
has_mg_ub <- rep(NA, length(scs))

for (i in 1:length(scs)) {
  has_mg <- logistic(post$Intercept + post$Mean * scs[i])
  has_mg_mean[i] <- mean(has_mg)
  has_mg_sd[i] <- sd(has_mg)
  has_mg_lb[i] <- HPDI(has_mg)[1]
  has_mg_ub[i] <- HPDI(has_mg)[2]
}

# calculation posterior distributions on p for each case
d$pr_mg_mean <- NA
d$pr_mg_sd <- NA
d$pr_mg_lb <- NA
d$pr_mg_ub <- NA

for (i in 1:nrow(d)) {
  has_mg <- logistic(
    post$Intercept +
    post$Mean * d$Mean[i] +
    post$Space * d$Space[i]
  )
  d$pr_mg_mean[i] <- mean(has_mg)
  d$pr_mg_sd[i] <- sd(has_mg)
  d$pr_mg_lb[i] <- HPDI(has_mg)[1]
  d$pr_mg_ub[i] <- HPDI(has_mg)[2]
}

plot(scs, has_mg_mean, ylim = c(0, 1), type = "l",
  ylab = "pr(moralizing gods)", xlab = "social complexity",
  main = "analysis without imputed outcomes")

polygon(c(scs, rev(scs)), c(has_mg_ub, rev(has_mg_lb)),
  border = NA, col = col.alpha("dodgerblue", 0.2))

all(appearances %in% d$key) # use `appearances` above to id cases we want to predict
tar <- which(d$key %in% appearances)
print(paste(length(tar), "cases with same IDs in new version"))
points(d$Mean[tar], d$pr_mg_mean[tar], pch = 16,
  col = col_alpha("black", 0.5))

# for (i in 1:length(tar)) lines(c(d$Mean[tar[i]], d$Mean[tar[i]]),
#   c(d$pr_mg_lb[tar[i]], d$pr_mg_ub[tar[i]]))

dev.off()


print("regressions fit")
