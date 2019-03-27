
rm(list = ls())

source("../project_support.r")

dir_init("./temp")


# tapply(d$PolPop, d$MG_missing, mean)
# tapply(d$Writing, d$MG_missing, mean)
# tapply(d$Money, d$MG_missing, mean)

# 92% of present had writing
# (the remaining were documented by visitors)
# 84% of the missingness were pre-literate

# states w/o data avg was 375,000
# states with data avg pop of 1.7 million

# 

# or were these values made up too??

# the NAs are a problem here too, aren't they?

# now plot

load("./input/m1.rdata")

d <- read.csv("./input/RegrDat.csv", stringsAsFactors = FALSE)

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
  if(!is.na(d$Lag1[i]) & !is.na(d$Lag2[i])) {
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
}




png("./temp/m1_predictions_missingness.png", res = 300, height = 5, width = 5, units = "in")

plot(scs, has_mg_mean, ylim = c(0, 1), type = "l",
  ylab = "pr(moralizing gods)", xlab = "social complexity",
  main = "original analysis of MG innovation", xaxt = "n")
polygon(c(scs, rev(scs)), c(has_mg_ub, rev(has_mg_lb)),
  border = NA, col = col.alpha("dodgerblue", 0.2))

axis(1, at = seq(0, 1, by = 0.2),
  labels = seq(0, 1, by = 0.2))

tar <- which(d$MG_missing == 1)
for (i in 1:length(tar)) lines(c(d$Mean[tar[i]], d$Mean[tar[i]]),
  c(d$pr_mg_lb[tar[i]], d$pr_mg_ub[tar[i]]), col = col.alpha("dodgerblue", 0.05))

points(d$Mean[tar], d$pr_mg_mean[tar], pch = 16, col = "dodgerblue", cex = 0.6)


# now add missingness plot picture-in-picture

## add interior plot

## bounding box of interior figure
xmin <- 0.15
xmax <- 0.62
ymin <- 0.45
ymax <- 0.95

rect(xmin, xmax, ymin, ymax, border = NA, col = "white")

# prep variables on the new scale

d$logPop <- d$PolPop
d$Pop <- 10^d$logPop

d$Mean2 <- (d$Mean - 0) / 1 * (xmax - xmin) + xmin
d$logPop2 <- (d$logPop - min(d$logPop[tar])) / max(d$logPop[tar]) * (ymax - ymin) + ymin

d$na_col <- ifelse(d$MG_missing == 1, "gray", "black")

# plot

mg <- 0.02 # don't mix up with outcome variable `MG`!

tar <- which(d$MG_missing == 1)
points(d$Mean2[tar], d$logPop2[tar], pch = 16, col = col_alpha(d$na_col[tar], 0.8), cex = 0.6)
tar <- which(d$MG_missing == 0)
points(d$Mean2[tar], d$logPop2[tar], pch = 16, col = col_alpha(d$na_col[tar], 0.8), cex = 0.6)

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

dev.off()



# plot predictions from m2, which could not see the MG that are NA 

load("./input/m2.rdata")

d <- read.csv("./input/RegrDat.csv", stringsAsFactors = FALSE)

# 17 "0s", 319 "1s"

NGAs <- c(
  "Big Island Hawaii",         "Cambodian Basin",          
  "Central Java",              "Chuuk Islands",            
  "Cuzco",                     "Deccan",                   
  "Garo Hills",                "Ghanaian Coast",           
  "Iceland",                   "Kachi Plain",              
  "Kansai",                    "Kapuasi Basin",            
  "Konya Plain",               "Latium",                   
  "Lena River Valley",         "Lowland Andes",            
  "Middle Yellow River Valley","Niger Inland Delta",       
  "North Colombia",            "Orkhon Valley",            
  "Oro PNG",                   "Paris Basin",              
  "Sogdiana",                  "Susiana",                  
  "Upper Egypt",               "Yemeni Coastal Plain"    
)

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


# now predict on the focal 12 NGAs over their existences using m2

# calculate "time of first appearance" and subtract off

NGAs <- c("Upper Egypt", "Susiana", "Konya Plain",
  "Middle Yellow River Valley", "Kachi Plain", "Sogdiana",
  "Latium", "Deccan", "Paris Basin", "Orkhon Valley", 
  "Kansai", "Niger Inland Delta")

d$time_fa <- NA

for (i in 1:length(NGAs)) {
  my_rows <- which(d$NGA == NGAs[i])
  my_first_mg_row <- my_rows[min(which(d$MG[my_rows] == 1))]
  d$time_fa[my_rows] <- d$Time[my_rows] - d$Time[my_first_mg_row]
}

png("./temp/revised_EDfit1.png", res = 300, height = 8, width = 10, units = "in")

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
  abline(h = 0.5, col = "red")
}

dev.off()


# now show the predictions for the missing values

drop <- which(dm$MG_missing == 0)
dm <- d[-drop, ]

png("./temp/m2_missingness_predictions.png", res = 300, height = 5, width = 5, units = "in")

plot(scs, has_mg_mean, ylim = c(0, 1), type = "l",
  ylab = "pr(moralizing gods)", xlab = "social complexity",
  main = "predicted probability of moral gods", xaxt = "n")
polygon(c(scs, rev(scs)), c(has_mg_ub, rev(has_mg_lb)),
  border = NA, col = col.alpha("dodgerblue", 0.2))

axis(1, at = seq(0, 1, by = 0.2) - 0.5, labels = seq(0, 1, by = 0.2))

for (i in 1:nrow(dm)) lines(c(dm$Mean_c[i], dm$Mean_c[i]),
  c(dm$pr_mg_lb[i], dm$pr_mg_ub[i]), col = col.alpha("black", 0.05))
points(dm$Mean_c, dm$pr_mg_mean, pch = 16, col = "black", cex = 0.6)

dev.off()



# define the evidence threshold for "first appearance" analysis

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

write.csv(nga_dat, "./temp/earliest_mg_estimates.csv", row.names = FALSE)

# its weird that min_year_90 for upper egypt is 800 years after first appearance








# aggregate first appearance estimates for revised fig2 from paper

SCNorm <- read.csv('./input/SCNorm.csv', stringsAsFactors = FALSE)
data <- read.csv('./input/PrePostComparison.csv', stringsAsFactors = FALSE)
nga_dat <- read.csv("./temp/earliest_mg_estimates.csv", stringsAsFactors = FALSE)

min_year_50_mean <- mean(nga_dat$min_year_50)
min_year_50_sd <- sd(nga_dat$min_year_50)/sqrt(nrow(nga_dat))


png("./temp/revised_fig2.png", res = 300, height = 5, width = 5, units = "in")

col1 <- rgb(0, 0, 0, max = 255, alpha = 50)

plot(SCNorm$x, SCNorm$Mean, type = "n", ylim = c(0, 1), xlim = c(-2000, 2000), ann = FALSE,
  xaxs = "i", yaxs = "i")

# MG appear period
rect(0, 0, 0 + mean(data$RangeMGAppear), 1,
  border = NA, col = col1)
# why would you do that? why not center the range at 0?

lines(SCNorm$x, SCNorm$Mean, type="l") 
lines(SCNorm$x, SCNorm$Upper, type="l",lty="dotted") 
lines(SCNorm$x, SCNorm$Lower, type="l",lty="dotted")

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




