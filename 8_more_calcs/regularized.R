
rm(list = ls())
source("../project_support.R")

# visualize the missingness problem

RegrDat <- read.csv("./input/RegrDat.csv", stringsAsFactors = FALSE)

d <- RegrDat

d$MG[d$MG_missing == 1] <- NA

NGAs <- sort(unique(RegrDat$NGA))
nga_col <- viridis(length(NGAs))

d$nga_col <- nga_col[match(d$NGA, NGAs)]

plot(d$Time, d$Mean, type = "n")

for (i in 1:length(NGAs)) {
  my_rows <- which(d$NGA == NGAs[i])
  points(d$Time[my_rows], d$Mean[my_rows],
    col = d$nga_col[my_rows], type = "l")
}

d$ybp <- 2000 - d$Time

plot(d$ybp, d$Mean, type = "n",
  col = d$nga_col, log = "x",
  xlim = c(12000, 100)) # cool trick!

for (i in 1:length(NGAs)) {
  my_rows <- which(d$NGA == NGAs[i])
  points(d$ybp[my_rows], d$Mean[my_rows],
    col = d$nga_col[my_rows], type = "l")
}

library(dplyr)

d$mg_col <- case_when(
  d$MG == 1 ~ "black",
  d$MG == 0 ~ "white",
  is.na(d$MG) ~ "white"
)

d$mg_outline <- ifelse(is.na(d$MG), "gray", "black")

# 5 x 6 matrix with each pop to look at

d$time_nga <- NA

for (i in 1:length(NGAs)) {
  my_rows <- which(d$NGA == NGAs[i])
  min_year <- min(d$Time[my_rows])
  d$time_nga[my_rows] <- d$Time[my_rows] - min_year
}

png("test.png", res = 300, units = "in", height = 12, width = 16)

par(mfrow = c(5, 6))

for (i in 1:length(NGAs)) {

  my_rows <- which(d$NGA == NGAs[i])
  plot(d$Time[my_rows], d$PolPop[my_rows], ylim = c(1, 8.5),
    type = "n", main = NGAs[i], frame.plot = FALSE,
    ylab = "log10 population", xlab = "years common era",
    tck = 0.05, las = 1)

  points(d$Time[my_rows], d$PolPop[my_rows],
    col = nga_col[i], type = "l")
  if(length(my_rows) == 1) points(d$Time[my_rows], d$PolPop[my_rows],
    pch = 20, cex = 0.2, col = nga_col[i])

  axis(1, tck = 0.02, at = seq(min(d$Time[my_rows]), max(d$Time[my_rows]), by = 100),
    labels = FALSE)

  my_points <- which(d$NGA == NGAs[i] & d$MG == 1)
  points(d$Time[my_points], d$PolPop[my_points],
    bg = d$mg_col[my_points], col = d$mg_outline[my_points], pch = 21)
  my_points <- which(d$NGA == NGAs[i] & d$MG == 0)
  points(d$Time[my_points], d$PolPop[my_points],
    bg = d$mg_col[my_points], col = d$mg_outline[my_points], pch = 21)

}

dev.off()


keep <- rep(NA, length(NGAs))

for (i in 1:length(NGAs)) {
  my_rows <- which(d$NGA == NGAs[i])
  keep[i] <- any(d$MG[my_rows] == 1) & any(is.na(d$MG[my_rows] | d$MG[my_rows] == 0))
}

NGA_keep <- NGAs[which(keep)]



png("test2.png", res = 300, units = "in", height = 12, width = 8)

par(mfrow = c(5, 3))

plot(c(0, 10), c(0, 10), type = "n", frame.plot = FALSE, axes = FALSE,
 xlab = "", ylab = "", main = "Legend")
# ann = all axis labels and titles

symbol_line <- 0.5
symbol_anchor <- 8

lines(symbol_line + c(-0.3, 0.3), c(symbol_anchor, symbol_anchor))
points(symbol_line, symbol_anchor, pch = 21, cex = 1, bg = "black")
text(symbol_line + 0.5, symbol_anchor, "moralizing gods present", pos = 4)

lines(symbol_line + c(-0.3, 0.3), c(symbol_anchor, symbol_anchor) - 2)
points(symbol_line, symbol_anchor - 2, pch = 21, cex = 1, bg = "white")
text(symbol_line + 0.5, symbol_anchor - 2, "moralizing gods absent", pos = 4)

lines(symbol_line + c(-0.3, 0.3), c(symbol_anchor, symbol_anchor) - 4)
text(symbol_line + 0.5, symbol_anchor - 4, "no data on moralizing gods", pos = 4)
points(symbol_line, symbol_anchor - 4, col = col.alpha("gray", 0.8), pch = 1)


for (i in 1:length(NGA_keep)) {

  my_rows <- which(d$NGA == NGA_keep[i])
  plot(d$Time[my_rows], d$PolPop[my_rows], ylim = c(1, 8.5),
    type = "n", main = NGA_keep[i], frame.plot = FALSE,
    ylab = "log10 population", xlab = "years common era",
    tck = 0.05, las = 1)

  points(d$Time[my_rows], d$PolPop[my_rows],
    col = nga_col[i], type = "l")
  if(length(my_rows) == 1) points(d$Time[my_rows], d$PolPop[my_rows],
    pch = 20, cex = 0.2, col = nga_col[i])

  axis(1, tck = 0.02, at = seq(min(d$Time[my_rows]), max(d$Time[my_rows]), by = 100),
    labels = FALSE)

  my_points <- which(d$NGA == NGA_keep[i] & is.na(d$MG))
  points(d$Time[my_points], d$PolPop[my_points],
    col = col.alpha("gray", 0.8), pch = 1)
  my_points <- which(d$NGA == NGA_keep[i] & d$MG == 1)
  points(d$Time[my_points], d$PolPop[my_points],
    bg = d$mg_col[my_points], col = d$mg_outline[my_points], pch = 21)
  my_points <- which(d$NGA == NGA_keep[i] & d$MG == 0)
  points(d$Time[my_points], d$PolPop[my_points],
    bg = d$mg_col[my_points], col = d$mg_outline[my_points], pch = 21)

}

dev.off()

# subset to a representative **8** for discussion

NGA_keep <- c("North Colombia", "Big Island Hawaii", "Valley of Oaxaca", "Orkhon Valley", "Latium")

png("test3.png", res = 300, units = "in", height = 4.8, width = 8)

par(mfrow = c(2, 3))

plot(c(0, 10), c(0, 10), type = "n", frame.plot = FALSE, axes = FALSE,
 xlab = "", ylab = "", main = "Legend")
# ann = all axis labels and titles

symbol_line <- 0.5
symbol_anchor <- 8

lines(symbol_line + c(-0.3, 0.3), c(symbol_anchor, symbol_anchor))
points(symbol_line, symbol_anchor, pch = 21, cex = 1, bg = "black")
text(symbol_line + 0.5, symbol_anchor, "moralizing gods present", pos = 4)

lines(symbol_line + c(-0.3, 0.3), c(symbol_anchor, symbol_anchor) - 2)
points(symbol_line, symbol_anchor - 2, pch = 21, cex = 1, bg = "white")
text(symbol_line + 0.5, symbol_anchor - 2, "moralizing gods absent", pos = 4)

lines(symbol_line + c(-0.3, 0.3), c(symbol_anchor, symbol_anchor) - 4)
text(symbol_line + 0.5, symbol_anchor - 4, "no data on moralizing gods", pos = 4)
points(symbol_line, symbol_anchor - 4, col = col.alpha("gray", 0.8), pch = 1)



for (i in 1:length(NGA_keep)) {

  my_rows <- which(d$NGA == NGA_keep[i])
  plot(d$Time[my_rows], d$PolPop[my_rows], ylim = c(1, 8.5),
    type = "n", main = NGA_keep[i], frame.plot = FALSE,
    ylab = "log10 population", xlab = "years common era",
    tck = 0.05, las = 1)

  points(d$Time[my_rows], d$PolPop[my_rows],
    col = nga_col[i], type = "l")
  if(length(my_rows) == 1) points(d$Time[my_rows], d$PolPop[my_rows],
    pch = 20, cex = 0.2, col = nga_col[i])

  axis(1, tck = 0.02, at = seq(min(d$Time[my_rows]), max(d$Time[my_rows]), by = 100),
    labels = FALSE)

  my_points <- which(d$NGA == NGA_keep[i] & is.na(d$MG))
  points(d$Time[my_points], d$PolPop[my_points],
    col = col.alpha("gray", 0.8), pch = 1)
  my_points <- which(d$NGA == NGA_keep[i] & d$MG == 1)
  points(d$Time[my_points], d$PolPop[my_points],
    bg = d$mg_col[my_points], col = d$mg_outline[my_points], pch = 21)
  my_points <- which(d$NGA == NGA_keep[i] & d$MG == 0)
  points(d$Time[my_points], d$PolPop[my_points],
    bg = d$mg_col[my_points], col = d$mg_outline[my_points], pch = 21)

}

dev.off()







# now redo with centering, regularization, pooling, etc.

rm(list = ls())
source("../project_support.R")

d <- read.csv("./input/RegrDat.csv", stringsAsFactors = FALSE)

drop <- which(is.na(d$Lag1) | is.na(d$Lag2))
d <- d[-drop, ]

d$key <- paste(d$NGA, d$Time)
appearances <- d$key[which(d$Lag1 == 0)]

original_model <- alist(
  MG <- dbinom(1, p),
  logit(p) <- a +
    b_sc * Mean +
    b_l1 * Lag1 + 
    b_l2 * Lag2 +
    b_ph * Phylogeny +
    b_sp * Space,
  a ~ dnorm(0, 7),
  b_sc ~ dnorm(0, 7),
  b_l1 ~ dnorm(0, 7),
  b_l2 ~ dnorm(0, 7),
  b_ph ~ dnorm(0, 7),
  b_sp ~ dnorm(0, 7)
)

dm <- d[, c("MG", "Mean", "Lag1", "Lag2", "Phylogeny", "Space")]
m1 <- map2stan(original_model, data = dm)






d <- read.csv("./input/RegrDat.csv", stringsAsFactors = FALSE)

d$key <- paste(d$NGA, d$Time)

# center social complexity at 0.5
d$Mean_c <- d$Mean - 0.5

dm <- d[, c("NGA", "MG", "Mean_c", "Phylogeny", "Space", "MG_missing")]

drop <- which(d$MG_missing == 1)
dm <- dm[-drop, ]

NGAs <- sort(unique(dm$NGA))
dm$nga <- match(dm$NGA, NGAs)

revised_model <- alist(
  MG <- dbinom(1, p),
  logit(p) <- a +
    a_nga[nga] +
    b_sc * Mean_c +
    b_ph * Phylogeny +
    b_sp * Space,
  a ~ dnorm(0, 1),
  a_nga[nga] ~ dnorm(0, a_sigma),
  a_sigma ~ dexp(1),
  b_sc ~ dnorm(0, 4),
  b_ph ~ dnorm(0, 4),
  b_sp ~ dnorm(0, 4)
)

m2 <- map2stan(revised_model, data = dm)




png("./replyfig_stan.png", res = 300, height = 5, width = 10, units = "in")

par(mfrow = c(1, 2))

d <- read.csv("./input/RegrDat.csv", stringsAsFactors = FALSE)

drop <- which(is.na(d$Lag1) | is.na(d$Lag2))
d <- d[-drop, ]

d$key <- paste(d$NGA, d$Time)
appearances <- d$key[which(d$Lag1 == 0)]

post <- extract.samples(m1)

scs <- seq(0, 1, by = 0.01)

has_mg_mean <- rep(NA, length(scs))
has_mg_sd <- rep(NA, length(scs))
has_mg_lb <- rep(NA, length(scs))
has_mg_ub <- rep(NA, length(scs))

for (i in 1:length(scs)) {
  has_mg <- logistic(post$a + post$b_sc * scs[i])
  has_mg_mean[i] <- mean(has_mg)
  has_mg_sd[i] <- sd(has_mg)
  has_mg_lb[i] <- HPDI(has_mg)[1]
  has_mg_ub[i] <- HPDI(has_mg)[2]
}

# identify cases to evaluate!!

d$pr_mg_mean <- NA
d$pr_mg_sd <- NA
d$pr_mg_lb <- NA
d$pr_mg_ub <- NA

for (i in 1:nrow(d)) {
  has_mg <- logistic(
    post$a +
    post$b_l1 * d$Lag1[i] +
    post$b_l2 * d$Lag2[i] +
    post$b_sc * d$Mean[i] +
    post$b_sp * d$Space[i] +
    post$b_ph * d$Phylogeny[i]
  )
  d$pr_mg_mean[i] <- mean(has_mg)
  d$pr_mg_sd[i] <- sd(has_mg)
  d$pr_mg_lb[i] <- HPDI(has_mg)[1]
  d$pr_mg_ub[i] <- HPDI(has_mg)[2]
}


# plot left figure

plot(scs, has_mg_mean, ylim = c(0, 1), type = "l",
  ylab = "pr(moralizing gods)", xlab = "social complexity",
  main = "original analysis of MG innovation", xaxt = "n")
polygon(c(scs, rev(scs)), c(has_mg_ub, rev(has_mg_lb)),
  border = NA, col = col.alpha("dodgerblue", 0.2))

axis(1, at = seq(0, 1, by = 0.2),
  labels = seq(0, 1, by = 0.2))

tar <- which(d$key %in% appearances)
points(d$Mean[tar], d$pr_mg_mean[tar], pch = 16, col = "black")

# for (i in 1:length(tar)) lines(c(d$Mean[tar[i]], d$Mean[tar[i]]),
#   c(d$pr_mg_lb[tar[i]], d$pr_mg_ub[tar[i]]), col = col.alpha("black", 0.2))

# same as publication figure

## add interior plot

## bounding box of interior figure
xmin <- 0.15
xmax <- 0.62
ymin <- 0.45
ymax <- 0.95

# prep variables on the new scale

d$logPop <- d$PolPop
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











d <- read.csv("./input/RegrDat.csv", stringsAsFactors = FALSE)

d$key <- paste(d$NGA, d$Time)

# center social complexity at 0.5
d$Mean_c <- d$Mean - 0.5

post <- extract.samples(m2)

scs <- seq(0, 1, by = 0.01) - 0.5

has_mg_mean <- rep(NA, length(scs))
has_mg_sd <- rep(NA, length(scs))
has_mg_lb <- rep(NA, length(scs))
has_mg_ub <- rep(NA, length(scs))

for (i in 1:length(scs)) {
  has_mg <- logistic(post$a + post$b_sc * scs[i])
  has_mg_mean[i] <- mean(has_mg)
  has_mg_sd[i] <- sd(has_mg)
  has_mg_lb[i] <- HPDI(has_mg)[1]
  has_mg_ub[i] <- HPDI(has_mg)[2]
}

d$pr_mg_mean <- NA
d$pr_mg_sd <- NA
d$pr_mg_lb <- NA
d$pr_mg_ub <- NA

for (i in 1:nrow(d)) {
  if(d$NGA[i] %in% NGAs) {
    nga <- match(d$NGA[i], NGAs)
    nga_offset <- post$a[nga]
  } else {
    nga_offset <- 0
  }
  has_mg <- logistic(
    post$a +
    nga_offset +
    post$b_sc * d$Mean_c[i] +
    post$b_sp * d$Space[i] +
    post$b_ph * d$Phylogeny[i]
  )
  d$pr_mg_mean[i] <- mean(has_mg)
  d$pr_mg_sd[i] <- sd(has_mg)
  d$pr_mg_lb[i] <- HPDI(has_mg)[1]
  d$pr_mg_ub[i] <- HPDI(has_mg)[2]
}

plot(scs, has_mg_mean, ylim = c(0, 1), type = "l",
  ylab = "pr(moralizing gods)", xlab = "social complexity",
  main = "analysis without missing values", xaxt = "n")
polygon(c(scs, rev(scs)), c(has_mg_ub, rev(has_mg_lb)),
  border = NA, col = col.alpha("dodgerblue", 0.2))

axis(1, at = seq(0, 1, by = 0.2) - 0.5, labels = seq(0, 1, by = 0.2))

tar <- which(d$key %in% appearances)
points(d$Mean_c[tar], d$pr_mg_mean[tar], pch = 16, col = "black")

# for (i in 1:length(tar)) lines(c(d$Mean_c[tar[i]], d$Mean_c[tar[i]]),
#   c(d$pr_mg_lb[tar[i]], d$pr_mg_ub[tar[i]]), col = col.alpha("black", 0.2))


dev.off()




# now predict on the famous 12 over their existences

NGAs <- c("Upper Egypt", "Susiana", "Konya Plain",
  "Middle Yellow River Valley", "Kachi Plain", "Sogdiana",
  "Latium", "Deccan", "Paris Basin", "Orkhon Valley", 
  "Kansai", "Niger Inland Delta")

# calculate "time of first appearance" and subtract off

d$time_fa <- NA

for (i in 1:length(NGAs)) {
  my_rows <- which(d$NGA == NGAs[i])
  my_first_mg_row <- my_rows[min(which(d$MG[my_rows] == 1))]
  d$time_fa[my_rows] <- d$Time[my_rows] - d$Time[my_first_mg_row]
}


# we need to calculate the 

pr_threshold <- 0.5
density_threshold <- 0.97

d$hit <- NA

for (i in 1:nrow(d)) {
  if(d$NGA[i] %in% NGAs) {
    nga <- match(d$NGA[i], NGAs)
    nga_offset <- post$a[nga]
  } else {
    nga_offset <- 0
  }
  has_mg <- logistic(
    post$a +
    nga_offset +
    post$b_sc * d$Mean_c[i] +
    post$b_sp * d$Space[i] +
    post$b_ph * d$Phylogeny[i]
  )
  d$hit[i] <- mean(has_mg > pr_threshold) > density_threshold
}

min_year_50 <- rep(NA, length(NGAs))

for(i in 1:length(NGAs)) {
  dm <- d[which(d$NGA == NGAs[i]),]
  min_year_50[i] <- min(dm$time_fa[dm$hit == 1])
}



pr_threshold <- 0.9
density_threshold <- 0.97

d$hit <- NA

for (i in 1:nrow(d)) {
  if(d$NGA[i] %in% NGAs) {
    nga <- match(d$NGA[i], NGAs)
    nga_offset <- post$a[nga]
  } else {
    nga_offset <- 0
  }
  has_mg <- logistic(
    post$a +
    nga_offset +
    post$b_sc * d$Mean_c[i] +
    post$b_sp * d$Space[i] +
    post$b_ph * d$Phylogeny[i]
  )
  d$hit[i] <- mean(has_mg > pr_threshold) > density_threshold
}

min_year_90 <- rep(NA, length(NGAs))

for(i in 1:length(NGAs)) {
  dm <- d[which(d$NGA == NGAs[i]),]
  min_year_90[i] <- min(dm$time_fa[dm$hit == 1])
}

nga_dat <- data.frame(NGA = NGAs, min_year_50, min_year_90)

write.csv(nga_dat, "earliest_mg_estimates.csv", row.names = FALSE)


png("revised_EDfit1.png", res = 300, height = 8, width = 10, units = "in")

par(mfrow = c(3, 4))

for(i in 1:length(NGAs)) {
  dm <- d[which(d$NGA == NGAs[i]),]
  plot(dm$time_fa, dm$pr_mg_mean, ylim = c(0, 1),
    xlim = c(-4000, 100), type = "l",
    xlab ="years before first apperance", 
    ylab = "pr(moralizing gods present)",
    main = NGAs[i])

  polygon(c(dm$time_fa, rev(dm$time_fa)), c(dm$pr_mg_lb, rev(dm$pr_mg_ub)),
    border = NA, col = col.alpha("firebrick", 0.2))
  # tar <- which(dm$hit == 1)
  # points(dm$time_fa[tar], dm$pr_mg_lb[tar], pch = 20, col = "red")
  abline(h = 0.5, col = "red")
}

dev.off()


# more than 90% sure more likely than not

# now make that fukcing fig again



rm(list = ls())
source("../project_support.R")

dir_init("./temp")

polities <- read.csv('./input/polities.csv', header=TRUE)

#New scripts for automated analysis of rates of change in social complexity pre/post moralising gods/doctrinal mode/writing

dat <- read.table("./input/PC1_traj_merged.csv", sep=",", header=TRUE) # from 3_run_pca
dat$NGA<-as.character(dat$NGA)
NGAs <- levels(polities$NGA)
NGAs <- NGAs[NGAs != "Crete"]    #### Remove new NGAs
NGAs <- NGAs[NGAs != "Galilee"]

#Overall rate (beginning to end of polity)
out <- matrix(NA, nrow=0, ncol=4)
for(i in 1:length(NGAs)){
  dt <- dat[dat$NGA == NGAs[i],]
  Earliest<-subset(dt,Time==min(Time))
  Latest<-subset(dt,Time==max(Time))
  MG<-subset(dt,MoralisingGods=="1") #Replace "MoralisingGods" with "DoctrinalMode" or "Writing" to do these analyses
  MGAppear<-subset(MG, Time==min(Time))
  rates<-cbind(MGAppear$NGA,
    (MGAppear$Mean-Earliest$Mean)/(MGAppear$Time-Earliest$Time),
    ((Latest$Mean-MGAppear$Mean)/(Latest$Time-MGAppear$Time)),
    (MGAppear$End-MGAppear$Start))
  out <- rbind(out,rates)
}
colnames(out)<-c("NGA","PreRate","PostRate","MGUncertainty")
#mean(out$MGUncertainty) #Use this while replacing "MoralisingGods" as above to get uncertainty values for time-series

# bret: wat happens when you have no moralizing god observations? they just drop exclude the whole NGA
# 10 NGAs drop out as a result, since they never have MG
# 7 have no pre-rate (so they start with MG)
# 1 has no post-rate (hawaii it has MG only in the last obs)
# so from 30 NGAs we are down to 12 with a pre and post rate
# and of those

# Of 30 NGAs, 10 are excluded from the t-test because no MG ever. Then, 7 more are gone because MG = 1 in all observations, then Hawaii is out because MG=1 only in the last observation (no "post-MG rate"). So the surviving 12 NGAs are in your table Joe, and in every case but Middle Yellow River the MG goes from NA to 1.

write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE) #Exporting/importing to force it to read as numeric (there is probably a more elegant way to do this)
out<-read.table("./temp/FullRates.csv", sep=",", header=TRUE)

out$Difference<-out[,2]-out[,3]
write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE)


#Full time windows (up to 10,000 years before and after moralizing gods)
NGAs <- levels(out$NGA)
out <- matrix(NA, nrow=0, ncol=5)
for(i in 1:length(NGAs)){
  dt <- dat[dat$NGA == NGAs[i],]
  MG<-subset(dt,MoralisingGods=="1") #Replace "MoralisingGods" with "DoctrinalMode" or "Writing" to do these analyses
  #library(dplyr)
  #MG<-as.data.frame(MG %>% group_by(PolID) %>% sample_n(size = 1)) #randomly samples so there is only one century per polity
  MGAppear<-subset(MG, Time==min(Time))
  for(j in 1: 100){
    Earliest<-subset(dt, Time==MGAppear$Time-j*100) 
    Latest<-subset(dt, Time==MGAppear$Time+j*100)
    rates<-cbind(MGAppear$NGA,ifelse(class(Earliest$Time)=="NULL","NA",(MGAppear$Mean-Earliest$Mean)/(MGAppear$Time-Earliest$Time)),ifelse(class(Latest$Time)=="NULL","NA",((Latest$Mean-MGAppear$Mean)/(Latest$Time-MGAppear$Time))),(MGAppear$End-MGAppear$Start),j*100)
    out <- rbind(out,rates)
  }
  out <- rbind(out,rates)
}
colnames(out)<-c("NGA","PreRate","PostRate","MGUncertainty","TimeWindow")

#mean(out$MGUncertainty) #Use this while replacing "MoralisingGods" as above to get uncertainty values for time-series

write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE) #Exporting/importing to force it to read as numeric (there is probably a more elegant way to do this)
out<-read.table("./temp/FullRates.csv", sep=",", header=TRUE)

out$Difference<-out[,3]-out[,2]

for(i in 2:length(out[,5])){
  out[i,7]<-out[i,5]-out[i-1,5]
}
out <-subset(out, out[,7]!=0) #getting rid of bug when the final row repeats in each NGA

write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE)

out <-subset(out, out[,5]<2050) #Change this to modify time-window restriction from 700 years pre/post moralizing gods (<750) or # out to use full time-window

write.csv(out, file="./temp/EqualRates.csv",  row.names=FALSE)

#bar chart paired
my.values<-mean(1000*out[,6],na.rm=TRUE)
err1<-1.96*std.error(1000*out[,6],na.rm=TRUE)
x <- barplot(my.values,names.arg="After - before moralizing gods",ylim = c(-1.3,1),ylab="Rate of increase in social complexity (SC/ky)")
arrows(x,my.values-err1 ,x,my.values+err1, code=3, angle=90, length=.1)
abline(h=0)

#significance test
print(t.test(out[,3], out[,2],paired=TRUE))

#histogram of differences
hist(1000*out[,6],col="gray",breaks="FD",xlim=c(-15,5),ylim=c(0,80),xlab="Rates of change in social complexity before vs. after moralizing gods (SC/ky)", ylab="Frequency")
abline(v = 0,col="black")

#Mean Pre:post ratio
print(mean(out[,2],na.rm=TRUE)/mean(out[,3],na.rm=TRUE))

###Calculate p-values etc.


data <- matrix(NA, nrow=0, ncol=5)

for(i in 1:length(NGAs)){
  dt <- out[out$NGA == NGAs[i],]
  ot <- dat[dat$NGA == NGAs[i],]
  MG<-subset(ot,MoralisingGods=="1") #Replace "MoralisingGods" with "DoctrinalMode" or "Writing" to do these analyses
  #library(dplyr)
  #MG<-as.data.frame(MG %>% group_by(PolID) %>% sample_n(size = 1)) #randomly samples so there is only one century per polity
  MGAppear<-subset(MG, Time==min(Time))
  DM<-subset(ot,DoctrinalMode=="1") #Replace "MoralisingGods" with "DoctrinalMode" or "Writing" to do these analyses
  DMAppear<-subset(DM, Time==min(Time))
  WR<-subset(ot,Writing=="1") #Replace "MoralisingGods" with "DoctrinalMode" or "Writing" to do these analyses
  WRAppear<-subset(WR, Time==min(Time))
  my.values<-c(NGAs[i],ifelse(class(mean(dt[,3],na.rm=TRUE))=="NULL","NA",mean(dt[,3],na.rm=TRUE)),ifelse(class(mean(dt[,2],na.rm=TRUE))=="NULL","NA",mean(dt[,2],na.rm=TRUE)),ifelse(class(std.error(dt[,3],na.rm=TRUE))=="NULL","NA",1.96*std.error(dt[,3],na.rm=TRUE)),ifelse(class(std.error(dt[,2],na.rm=TRUE))=="NULL","NA",1.96*std.error(dt[,2],na.rm=TRUE)))
  data <- rbind(data,my.values)
}
colnames(data)<-c("NGA","PostRate","PreRate","PostConfInt","PreConfInt")
write.csv(data, file="./temp/PrePostComparisonFull.csv",  row.names=FALSE)
data<-read.table("./temp/PrePostComparisonFull.csv", sep=",", header=TRUE)
data<-as.data.frame(data)
data$Difference<-data[,2]-data[,3]
data[,2:6]<-data[,2:6]*1000
write.csv(data, file="./temp/PrePostComparisonFull.csv",  row.names=FALSE)



#Full values for matching pre-/post-NGAs
out <-subset(out, out[,6]<1000) #Removing windows without matching pre-/post-MG rates
write.csv(out, file="./temp/FullRates.csv",  row.names=FALSE)

# oh snap, they don't include all the NGAs in this figure...

out <- read.table("./temp/FullRates.csv", sep=",", header=TRUE)
NGAs <- levels(out$NGA)

data <- matrix(NA, nrow=0, ncol=15)

for(i in 1:length(NGAs)){
  dt <- out[out$NGA == NGAs[i],]
  ot <- dat[dat$NGA == NGAs[i],]
  MG<-subset(ot,MoralisingGods=="1") #Replace "MoralisingGods" with "DoctrinalMode" or "Writing" to do these analyses
  #library(dplyr)
  #MG<-as.data.frame(MG %>% group_by(PolID) %>% sample_n(size = 1)) #randomly samples so there is only one century per polity
  MGAppear<-subset(MG, Time==min(Time))
  DM<-subset(ot,DoctrinalMode=="1") #Replace "MoralisingGods" with "DoctrinalMode" or "Writing" to do these analyses
  DMAppear<-subset(DM, Time==min(Time))
  WR<-subset(ot,Writing=="1") #Replace "MoralisingGods" with "DoctrinalMode" or "Writing" to do these analyses
  WRAppear<-subset(WR, Time==min(Time))
  my.values<-c(NGAs[i], mean(dt[,3]), mean(dt[,2]),
    1.96 * std.error(dt[,3]), 1.96 * std.error(dt[,2]),
    t.test(dt[,3],dt[,2])$p.value,
    t.test(dt[,3],dt[,2])$parameter, length(dt[,2]),
    t.test(dt[,3],dt[,2])$statistic, DMAppear$Time - MGAppear$Time,
    WRAppear$Time - MGAppear$Time,
    MGAppear$End - MGAppear$Start, # the Seshat historical range w/in which the first appearance lies
    DMAppear$End - DMAppear$Start,
    WRAppear$End - WRAppear$Start,
    MGAppear$Time
  )
  data <- rbind(data,my.values)
}

colnames(data)<-c("NGA","PostRate","PreRate","PostConfInt","PreConfInt","p","df","n","t","PreMGRitual","PreMGWriting","RangeMGAppear","RangeDMAppear","RangeWRAppear","MGAppear")
write.csv(data, file="./temp/PrePostComparison.csv",  row.names=FALSE)
data<-read.table("./temp/PrePostComparison.csv", sep=",", header=TRUE)
data<-as.data.frame(data)
data$Difference<-data[,2]-data[,3]
data[,c(2:5,16)]<-data[,c(2:5,16)]*1000
write.csv(data, file="./temp/PrePostComparison.csv",  row.names=FALSE)

#Test signifiance of doctrinal ritual preceding moralizing gods
print(t.test(data$PreMGRitual))

######Normalize time-series centered around moralising god appearance
out <- matrix(NA, nrow=0, ncol=0)
for(i in 1:length(NGAs)){
  dt <- dat[dat$NGA == NGAs[i],]
  Earliest<-subset(dt,Time==min(Time))  
  Latest<-subset(dt,Time==max(Time))  
  MG<-subset(dt,MoralisingGods=="1") #Replace "MoralisingGods" with "DoctrinalMode" or "Writing" to do these analyses
  MGAppear<-subset(MG, Time==min(Time))
  dt$Time.norm<-dt$Time-MGAppear$Time
  out <- rbind(out,dt)
}
out$MGUncertainty<-out$End-out$Start
write.csv(out, file="./temp/TimeNorm.csv",  row.names=FALSE) 

#Merge Normalized times
dat <- read.table("./temp/TimeNorm.csv", sep=",", header=TRUE)
out<-unique(dat$Time.norm)

#library(plyr)
#library(dplyr) ##Bugs when I load dplyr and plyr (even when only loading dplyr after plyr)!

for(i in 1:length(NGAs)){
  dt <- dat[dat$NGA == NGAs[i],]
  out<-merge(out,dt[,c("Time.norm","Mean")],by.x="x",by.y="Time.norm",all=TRUE)
  out<-rename(out, c("Mean" = NGAs[i]))
}

MeanPCs<-out[,2:(1+length(NGAs))]

out$Mean <- apply(MeanPCs,1,mean,na.rm=TRUE)
out$Lower <- out$Mean - 1.96*apply(MeanPCs,1, std.error,na.rm=TRUE)
out$Upper <- out$Mean + 1.96*apply(MeanPCs,1, std.error,na.rm=TRUE)
write.csv(out, file="./temp/SCNorm.csv",  row.names=FALSE) 








#plot normalized times
FullImpDat<-read.csv('./temp/SCNorm.csv', header=TRUE)

####Plot 12 NGA time-series with pre- and post-moralizing god SC data
#First average normalized to moralizing gods time
out<-read.csv('./temp/SCNorm.csv', header=TRUE)
data<-read.csv('./temp/PrePostComparison.csv', header=TRUE)

#data <- FullImpDat
y <- out$Mean
x <- out$x
l <- out$Lower
u <- out$Upper
#earliest<-read.csv('EarliestAppearance.csv', header=TRUE)
#a<-earliest[,5]
#b<-earliest[,6]
#c<-earliest[,2]
#d<-earliest[,4]

#ylim <- c(min(na.rm=TRUE,y), max(na.rm=TRUE,y)) #This sets y-axis from the minimum to maximum values throughout the whole global dataset
ylim <- c(0,1) #This sets x-axis to min (0) to max (1) social complexity
#xlim <- c(min(na.rm=TRUE,x), max(na.rm=TRUE,x)) #Use this instead to set x-axis from earliest to latest time in dataset
xlim <- c(-2000,2000) #This sets x-axis to 1,000 years before and after the appearance of moralising gods
pch=19
cex=.5
xaxt='s'
yaxt='s'
xaxs="i"
yaxs="i"
ann=FALSE
#v1=min(subset(x,a==1))
#v2=min(subset(x,c==1))
linecol1<-"red"
linecol2<-"green"
linecol3<-"blue"
linecol4<-"orange"
lty1<-2
lty2<-4
lty3<-3
lty4<-2
lwd<-1.5
type="l"
h=.6

col1<-rgb(0,0,0,max=255,alpha=50)
#col2<-rgb(255,0,0,max=255,alpha=125)
#col3<-rgb(0,255,0,max=255,alpha=125)

png("revised_fig2.png", res = 300, height = 5, width = 5, units = "in")

plot(x, y, ylim=ylim, xlim=xlim, pch=pch, cex=cex, xaxt=xaxt, yaxt=yaxt,type=type,xaxs=xaxs,yaxs=yaxs,
  main = "avg. emergence of MG by SC, no NA to 0", xlab = "years until first record of MG", ylab = "social complexity")

panel.first = rect(
  0,
  -1e6,
  0 + mean(data$RangeMGAppear),
  1e6,
  col=col1, border=NA
)

#panel.first = rect(c(mean(data$PreMGRitual)-(1.96*std.error(data$PreMGRitual)), mean(data$PreMGWriting)-(1.96*std.error(data$PreMGWriting)), 0), -1e6, c(mean(data$PreMGRitual)+ (1.96*std.error(data$PreMGRitual)), mean(data$PreMGWriting)+(1.96*std.error(data$PreMGWriting)), 0+mean(data$RangeMGAppear)), 1e6, col=c(col1,col2,col3), border=NA)
#This was earlier version that included writing
#polygon(c(x, rev(x)) , c(u, rev(l)) , col = 'grey' , border = NA) #out for now because of bugs
lines(x, y,type="l") 
lines(x, u,type="l",lty="dotted") 
lines(x, l,type="l",lty="dotted")
#abline(h=0.6,lty="dashed")



nga_dat <- read.csv("./earliest_mg_estimates.csv", stringsAsFactors = FALSE)

min_year_50_mean <- mean(nga_dat$min_year_50)
min_year_50_sd <- sd(nga_dat$min_year_50)/sqrt(nrow(nga_dat))

polygon(
  c(min_year_50_mean + c(-1.96, 1.96) * min_year_50_sd, 
  min_year_50_mean + c(1.96, -1.96) * min_year_50_sd),
  c(0, 0, 1, 1),
  border = NA,
  col = col.alpha("firebrick", 0.3)
)

abline(v = min_year_50_mean, col = "firebrick")

text(-1000, 0.7, "predicted MG\nemergence", srt = 90)
text(60, 0.3, "MG first recorded", srt = 90)

dev.off()




src <- read.csv("./input/HighGodSource.csv", stringsAsFactors = FALSE)

