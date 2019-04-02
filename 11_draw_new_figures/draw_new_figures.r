
rm(list = ls())

source("../project_support.r")

dir_init("./temp")

# now plot

load("./input/m1.rdata")

d <- read.csv("./input/RegrDat.csv", stringsAsFactors = FALSE)

# add centering
d$Mean_c <- d$Mean - 0.5

d$na_col <- ifelse(d$MG_missing == 1, "firebrick", "dodgerblue")

NGAs_full <- sort(unique(d$NGA))
NGA_cols <- sample(viridis(length(NGAs_full)))

d$nga_col <- NGA_cols[match(d$NGA, NGAs_full)]

post <- extract.samples(m1)

# plot counterfactual estimates for sc/mg relationship
# centered on average levels for missing values

ph_na_mean <- mean(d$Phylogeny[d$MG_missing == 1])
sp_na_mean <- mean(d$Space[d$MG_missing == 1])

scs <- seq(0, 1, by = 0.01) - 0.5
has_mg_mean <- rep(NA, length(scs))
has_mg_sd <- rep(NA, length(scs))
has_mg_lb <- rep(NA, length(scs))
has_mg_ub <- rep(NA, length(scs))

for (i in 1:length(scs)) {
  logit_p <- post$a +
    post$b_l1 * 0 +
    post$b_l2 * 0 +
    post$b_sc * scs[i] +
    post$b_sp * sp_na_mean +
    post$b_ph * ph_na_mean
  has_mg <- logistic(logit_p)
  has_mg_mean[i] <- mean(has_mg)
  has_mg_sd[i] <- sd(has_mg)
  has_mg_lb[i] <- HPDI(has_mg)[1]
  has_mg_ub[i] <- HPDI(has_mg)[2]
}


# summarize m1's posterior predictions for each case

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
      post$b_sc * d$Mean_c[i] +
      post$b_sp * d$Space[i] +
      post$b_ph * d$Phylogeny[i]
    )
    d$pr_mg_mean[i] <- mean(has_mg)
    d$pr_mg_sd[i] <- sd(has_mg)
    d$pr_mg_lb[i] <- HPDI(has_mg)[1]
    d$pr_mg_ub[i] <- HPDI(has_mg)[2]
  }
}

has_mg_mean_m1 <- has_mg_mean

png("./temp/m1_predictions_missingness.png", res = 300, height = 5, width = 5, units = "in")

plot(scs, has_mg_mean, ylim = c(0, 1), type = "l",
  ylab = "pr(moralizing gods)", xlab = "social complexity",
  main = "retrodicted probability of MG appearance", xaxt = "n")
polygon(c(scs, rev(scs)), c(has_mg_ub, rev(has_mg_lb)),
  border = NA, col = col_alpha("black", 0.2))

axis(1, at = seq(0, 1, by = 0.2) - 0.5,
  labels = seq(0, 1, by = 0.2))

# tar <- which(d$MG_missing == 1)
tar <- 1:nrow(d)
for (i in 1:length(tar)) lines(c(d$Mean_c[tar[i]], d$Mean_c[tar[i]]),
  c(d$pr_mg_lb[tar[i]], d$pr_mg_ub[tar[i]]), col = col_alpha(d$na_col[tar[i]], 0.05))

points(d$Mean_c[tar], d$pr_mg_mean[tar], pch = 16, col = d$na_col, cex = 0.6)

dev.off()


# plot predictions from m2, which could not see the MG that are NA 

load("./input/m2.rdata")

post <- extract.samples(m2)

scs <- seq(0, 1, by = 0.01) - 0.5

has_mg_mean <- rep(NA, length(scs))
has_mg_sd <- rep(NA, length(scs))
has_mg_lb <- rep(NA, length(scs))
has_mg_ub <- rep(NA, length(scs))

for (i in 1:length(scs)) {
  logit_p <- post$a +
    post$b_sc * scs[i] +
    post$b_sp * sp_na_mean +
    post$b_ph * ph_na_mean
  has_mg <- logistic(logit_p)
  has_mg_mean[i] <- mean(has_mg)
  has_mg_sd[i] <- sd(has_mg)
  has_mg_lb[i] <- HPDI(has_mg)[1]
  has_mg_ub[i] <- HPDI(has_mg)[2]
}

# these are the NGAs used in the m2 varying effects terms
NGAs_m2 <- c(
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

d$pr_mg_mean <- NA
d$pr_mg_sd <- NA
d$pr_mg_lb <- NA
d$pr_mg_ub <- NA

for (i in 1:nrow(d)) {
  if(d$NGA[i] %in% NGAs_m2) {
    nga <- match(d$NGA[i], NGAs_m2)
    nga_offset <- post$a_nga[, nga]
  } else {
    nga_offset <- rnorm(length(post$a_sigma), 0, post$a_sigma)
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


# now show the predictions for all missing values

png("./temp/m2_missingness_predictions.png", res = 300, height = 5, width = 5, units = "in")

plot(scs, has_mg_mean, ylim = c(0, 1), type = "l",
  ylab = "pr(moralizing gods)", xlab = "social complexity",
  main = "predicted probability of moral gods", xaxt = "n")
polygon(c(scs, rev(scs)), c(has_mg_ub, rev(has_mg_lb)),
  border = NA, col = col.alpha("dodgerblue", 0.2))

axis(1, at = seq(0, 1, by = 0.2) - 0.5, labels = seq(0, 1, by = 0.2))

# for (i in 1:nrow(d)) lines(c(d$Mean_c[tar[i]], d$Mean_c[tar[i]]),
#   c(d$pr_mg_lb[tar[i]], d$pr_mg_ub[tar[i]]), col = col.alpha("black", 0.05))

tar <- which(d$MG_missing == 1)
points(d$Mean_c[tar], d$pr_mg_mean[tar], pch = 16, col = col_alpha(d$nga_col[tar], 0.8), cex = 0.6)

points(scs, has_mg_mean_m1, col = "gray", type = "l", lty = 2)
text(0.3, 0.3, "original model", col = "gray", srt = 60)

dev.off()


# now predict on the focal 12 NGAs over their existences using m2

NGAs_short <- c("Upper Egypt", "Susiana", "Konya Plain",
  "Middle Yellow River Valley", "Kachi Plain", "Sogdiana",
  "Latium", "Deccan", "Paris Basin", "Orkhon Valley", 
  "Kansai", "Niger Inland Delta")

# calculate "time of first appearance" and subtract off

d$time_fa <- NA

for (i in 1:length(NGAs_short)) {
  my_rows <- which(d$NGA == NGAs_short[i])
  my_first_mg_row <- my_rows[min(which(d$MG[my_rows] == 1))]
  d$time_fa[my_rows] <- d$Time[my_rows] - d$Time[my_first_mg_row]
}

# define the evidence threshold for "first appearance" analysis

# `hit` = if 80% of the posterior density is over the p = 0.5 mark

pr_threshold <- 0.5
density_threshold <- 0.8

d$hit <- NA

for (i in 1:nrow(d)) {
  if(d$NGA[i] %in% NGAs_m2) {
    nga <- match(d$NGA[i], NGAs_m2)
    nga_offset <- post$a_nga[, nga]
  } else {
    nga_offset <- rnorm(length(post$a_sigma), 0, post$a_sigma)
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

min_year_50 <- rep(NA, length(NGAs_short))

for(i in 1:length(NGAs_short)) {
  dn <- d[which(d$NGA == NGAs_short[i] & d$time_fa < 0),]
  if (any(dn$hit == 1)) min_year_50[i] <- min(dn$time_fa[dn$hit == 1])
}

# `hit` = if 80% of the posterior density is over the p = 0.9 mark

pr_threshold <- 0.9
density_threshold <- 0.8

d$hit <- NA

for (i in 1:nrow(d)) {
  if(d$NGA[i] %in% NGAs_m2) {
    nga <- match(d$NGA[i], NGAs_m2)
    nga_offset <- post$a_nga[, nga]
  } else {
    nga_offset <- rnorm(length(post$a_sigma), 0, post$a_sigma)
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

min_year_90 <- rep(NA, length(NGAs_short))

for(i in 1:length(NGAs_short)) {
  dn <- d[which(d$NGA == NGAs_short[i] & d$time_fa < 0),]
  if (any(dn$hit == 1)) min_year_90[i] <- min(dn$time_fa[dn$hit == 1])
}

nga_dat <- data.frame(NGA = NGAs_short, min_year_50, min_year_90)

# visualize each NGA seperately of the 12

png("./temp/revised_EDfit1.png", res = 300, height = 8, width = 10, units = "in")

par(mfrow = c(3, 4))

for(i in 1:length(NGAs_short)) {
  dm <- d[which(d$NGA == NGAs_short[i]),]
  plot(dm$time_fa, dm$pr_mg_mean, ylim = c(0, 1),
    xlim = c(-4000, 100), type = "l",
    xlab ="years before first apperance", 
    ylab = "pr(moralizing gods present)",
    main = NGAs_short[i])

  abline(v = nga_dat$min_year_50[i], lwd = 2)

  polygon(c(dm$time_fa, rev(dm$time_fa)), c(dm$pr_mg_lb, rev(dm$pr_mg_ub)),
    border = NA, col = col.alpha("firebrick", 0.2))
  abline(h = 0.5, col = "red")
}

dev.off()

# aggregate first appearance estimates for revised fig2 from paper

SCNorm <- read.csv('./input/SCNorm.csv', stringsAsFactors = FALSE)
data <- read.csv('./input/PrePostComparison.csv', stringsAsFactors = FALSE)

min_year_50_mean <- mean(nga_dat$min_year_50, na.rm = TRUE)
min_year_50_sd <- sd(nga_dat$min_year_50, na.rm = TRUE)/sqrt(sum(!is.na(nga_dat$min_year_50)))


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




