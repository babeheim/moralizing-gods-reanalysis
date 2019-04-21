
rm(list = ls())

source("../project_support.r")

dir_init("./temp")

#####


print("load models and extract their posterior samples")

load("./input/m1.rdata")
m1_post <- extract.samples(m1)

load("./input/m1_alt1.rdata")
m1_alt1_post <- extract.samples(m1_alt1)

load("./input/m1_alt2.rdata")
m1_alt2_post <- extract.samples(m1_alt2)

load("./input/m2.rdata")
m2_post <- extract.samples(m2)



print("make m1 tables")

varnames <- c("Intercept", "Social Complexity", "Lag1", "Lag2",
  "Phylogeny", "Space", "N", "Deviance")
parnames <- c("a", "b_sc", "b_l1", "b_l2", "b_ph", "b_sp")

n <- ncol(m1_post$p)
dev <- sprintf("%.1f", -2 * mean(m1_post$lp__))
est <- sprintf("%.2f",  unlist(lapply(m1_post, mean))[parnames])
se <- sprintf("%.2f", sqrt(unlist(lapply(m1_post, var))[parnames]))
estse <- paste0(est, " (", se, ")")
psigns <- c("-", psign(m1_post$b_sc), psign(m1_post$b_l1),
  psign(m1_post$b_l2), psign(m1_post$b_ph), psign(m1_post$b_sp))
x <- data.frame(m1_est = c(estse, n, dev), m1_psigns = c(psigns, "", ""))

n <- ncol(m1_alt1_post$p)
dev <- sprintf("%.1f", -2 * mean(m1_alt1_post$lp__))
est <- sprintf("%.2f",  unlist(lapply(m1_alt1_post, mean))[parnames])
se <- sprintf("%.2f", sqrt(unlist(lapply(m1_alt1_post, var))[parnames]))
estse <- paste0(est, " (", se, ")")
psigns <- c("-", psign(m1_alt1_post$b_sc), psign(m1_alt1_post$b_l1),
  psign(m1_alt1_post$b_l2), psign(m1_alt1_post$b_ph), psign(m1_alt1_post$b_sp))
x$alt1_est <- c(estse, n, dev)
x$alt1_psigns <- c(psigns, "", "")

n <- ncol(m1_alt2_post$p)
dev <- sprintf("%.1f", -2 * mean(m1_alt2_post$lp__))
est <- sprintf("%.2f",  unlist(lapply(m1_alt2_post, mean))[parnames])
se <- sprintf("%.2f", sqrt(unlist(lapply(m1_alt2_post, var))[parnames]))
estse <- paste0(est, " (", se, ")")
psigns <- c("-", psign(m1_alt2_post$b_sc), psign(m1_alt2_post$b_l1),
  psign(m1_alt2_post$b_l2), psign(m1_alt2_post$b_ph), psign(m1_alt2_post$b_sp))
x$alt2_est <- c(estse, n, dev)
x$alt2_psigns <- c(psigns, "", "")

cell_cols <- rep(c("Est. (SE)", "P(sign)"), 3)
colnames(x) <- cell_cols
rownames(x) <- varnames

writeLines(texttab(x, hlines = c(1, 7)), "./temp/model1_variations.txt")



print("make m2 tables")

varnames <- c("Intercept", "Social Complexity", "Phylogeny",
  "Space", "NGA Varying Effect", "N", "Deviance")
parnames <- c("a", "b_sc", "b_ph", "b_sp", "a_sigma")

n <- ncol(m2_post$p)
dev <- sprintf("%.1f", -2 * mean(m2_post$lp__))
est <- sprintf("%.2f",  unlist(lapply(m2_post, mean))[parnames])
se <- sprintf("%.2f", sqrt(unlist(lapply(m2_post, var))[parnames]))
estse <- paste0(est, " (", se, ")")
psigns <- c("-", psign(m2_post$b_sc), psign(m2_post$b_ph),
  psign(m2_post$b_sp), "-")
x <- data.frame(m1_est = c(estse, n, dev), m1_psigns = c(psigns, "", ""))

colnames(x) <- c("Est. (SE)", "P(sign)")
rownames(x) <- varnames

writeLines(texttab(x, hlines = c(1, 6)), "./temp/model2.txt")



print("calculate posterior predictions and counterfactual predictions")

density_threshold <- 0.8

d <- read.csv("./input/RegrDat.csv", stringsAsFactors = FALSE)

# add centering
d$Mean_c <- d$Mean - 0.5

d$na_col <- ifelse(d$MG_missing == 1, "firebrick", "dodgerblue")

NGAs_full <- sort(unique(d$NGA))
NGA_cols <- viridis(length(NGAs_full))

d$nga_col <- NGA_cols[match(d$NGA, NGAs_full)]

# plot counterfactual estimates for sc/mg relationship
# centered on average levels for missing values

Mean_c <- seq(0, 1, by = 0.01) - 0.5
Lag1 <- 0
Lag2 <- 0
Phylogeny <- mean(d$Phylogeny[d$MG_missing == 1])
Space <- mean(d$Space[d$MG_missing == 1])

counter <- expand.grid(Mean_c, Lag1, Lag2, Phylogeny, Space)
colnames(counter) <- c("Mean_c", "Lag1", "Lag2", "Phylogeny", "Space")

# first do m1

counter$m1_pr_mg_mean <- NA
counter$m1_pr_mg_sd <- NA
counter$m1_pr_mg_lb <- NA
counter$m1_pr_mg_ub <- NA

for (i in 1:nrow(counter)) {
  logit_p <- m1_post$a +
    m1_post$b_l1 * counter$Lag1[i] +
    m1_post$b_l2 * counter$Lag2[i] +
    m1_post$b_sc * counter$Mean_c[i] +
    m1_post$b_sp * counter$Space[i] +
    m1_post$b_ph * counter$Phylogeny[i]
  pr_mg <- logistic(logit_p)
  counter$m1_pr_mg_mean[i] <- mean(pr_mg)
  counter$m1_pr_mg_sd[i] <- sd(pr_mg)
  counter$m1_pr_mg_lb[i] <- HPDI(pr_mg)[1]
  counter$m1_pr_mg_ub[i] <- HPDI(pr_mg)[2]
}

# summarize m1's posterior predictions for each case

d$m1_pr_mg_mean <- NA
d$m1_pr_mg_sd <- NA
d$m1_pr_mg_lb <- NA
d$m1_pr_mg_ub <- NA
d$m1_hit50 <- NA

for (i in 1:nrow(d)) {
  if (!is.na(d$Lag1[i]) & !is.na(d$Lag2[i])) {
    pr_mg <- logistic(
      m1_post$a +
      m1_post$b_l1 * d$Lag1[i] +
      m1_post$b_l2 * d$Lag2[i] +
      m1_post$b_sc * d$Mean_c[i] +
      m1_post$b_sp * d$Space[i] +
      m1_post$b_ph * d$Phylogeny[i]
    )
    d$m1_pr_mg_mean[i] <- mean(pr_mg)
    d$m1_pr_mg_sd[i] <- sd(pr_mg)
    d$m1_pr_mg_lb[i] <- HPDI(pr_mg)[1]
    d$m1_pr_mg_ub[i] <- HPDI(pr_mg)[2]
    d$m1_hit50[i] <- mean(pr_mg > 0.5) > density_threshold
  }
}

# now do this for m1_alt1

counter$m1_alt1_pr_mg_mean <- NA
counter$m1_alt1_pr_mg_sd <- NA
counter$m1_alt1_pr_mg_lb <- NA
counter$m1_alt1_pr_mg_ub <- NA

for (i in 1:nrow(counter)) {
  logit_p <- m1_alt1_post$a +
    m1_alt1_post$b_l1 * counter$Lag1[i] +
    m1_alt1_post$b_l2 * counter$Lag2[i] +
    m1_alt1_post$b_sc * counter$Mean_c[i] +
    m1_alt1_post$b_sp * counter$Space[i] +
    m1_alt1_post$b_ph * counter$Phylogeny[i]
  pr_mg <- logistic(logit_p)
  counter$m1_alt1_pr_mg_mean[i] <- mean(pr_mg)
  counter$m1_alt1_pr_mg_sd[i] <- sd(pr_mg)
  counter$m1_alt1_pr_mg_lb[i] <- HPDI(pr_mg)[1]
  counter$m1_alt1_pr_mg_ub[i] <- HPDI(pr_mg)[2]
}

# summarize m1_alt1's posterior predictions for each case

density_threshold <- 0.8

d$m1_alt1_pr_mg_mean <- NA
d$m1_alt1_pr_mg_sd <- NA
d$m1_alt1_pr_mg_lb <- NA
d$m1_alt1_pr_mg_ub <- NA
d$m1_alt1_hit50 <- NA

for (i in 1:nrow(d)) {
  if (!is.na(d$Lag1[i]) & !is.na(d$Lag2[i])) {
    pr_mg <- logistic(
      m1_alt1_post$a +
      m1_alt1_post$b_l1 * d$Lag1[i] +
      m1_alt1_post$b_l2 * d$Lag2[i] +
      m1_alt1_post$b_sc * d$Mean_c[i] +
      m1_alt1_post$b_sp * d$Space[i] +
      m1_alt1_post$b_ph * d$Phylogeny[i]
    )
    d$m1_alt1_pr_mg_mean[i] <- mean(pr_mg)
    d$m1_alt1_pr_mg_sd[i] <- sd(pr_mg)
    d$m1_alt1_pr_mg_lb[i] <- HPDI(pr_mg)[1]
    d$m1_alt1_pr_mg_ub[i] <- HPDI(pr_mg)[2]
    d$m1_alt1_hit50[i] <- mean(pr_mg > 0.5) > density_threshold
  }
}

# and for m1_alt2

counter$m1_alt2_pr_mg_mean <- NA
counter$m1_alt2_pr_mg_sd <- NA
counter$m1_alt2_pr_mg_lb <- NA
counter$m1_alt2_pr_mg_ub <- NA

for (i in 1:nrow(counter)) {
  logit_p <- m1_alt2_post$a +
    m1_alt2_post$b_l1 * counter$Lag1[i] +
    m1_alt2_post$b_l2 * counter$Lag2[i] +
    m1_alt2_post$b_sc * counter$Mean_c[i] +
    m1_alt2_post$b_sp * counter$Space[i] +
    m1_alt2_post$b_ph * counter$Phylogeny[i]
  pr_mg <- logistic(logit_p)
  counter$m1_alt2_pr_mg_mean[i] <- mean(pr_mg)
  counter$m1_alt2_pr_mg_sd[i] <- sd(pr_mg)
  counter$m1_alt2_pr_mg_lb[i] <- HPDI(pr_mg)[1]
  counter$m1_alt2_pr_mg_ub[i] <- HPDI(pr_mg)[2]
}

# summarize m1_alt2's posterior predictions for each case

d$m1_alt2_pr_mg_mean <- NA
d$m1_alt2_pr_mg_sd <- NA
d$m1_alt2_pr_mg_lb <- NA
d$m1_alt2_pr_mg_ub <- NA
d$m1_alt2_hit50 <- NA

for (i in 1:nrow(d)) {
  if (!is.na(d$Lag1[i]) & !is.na(d$Lag2[i])) {
    pr_mg <- logistic(
      m1_alt2_post$a +
      m1_alt2_post$b_l1 * d$Lag1[i] +
      m1_alt2_post$b_l2 * d$Lag2[i] +
      m1_alt2_post$b_sc * d$Mean_c[i] +
      m1_alt2_post$b_sp * d$Space[i] +
      m1_alt2_post$b_ph * d$Phylogeny[i]
    )
    d$m1_alt2_pr_mg_mean[i] <- mean(pr_mg)
    d$m1_alt2_pr_mg_sd[i] <- sd(pr_mg)
    d$m1_alt2_pr_mg_lb[i] <- HPDI(pr_mg)[1]
    d$m1_alt2_pr_mg_ub[i] <- HPDI(pr_mg)[2]
    d$m1_alt2_hit50[i] <- mean(pr_mg > 0.5) > density_threshold
  }
}



print("create a three-panel plot including each alternative")

png("./temp/alternative_missingness.png", res = 300,
  height = 3, width = 8, units = "in")

par(mfrow = c(1, 3))

plot(counter$Mean_c, counter$m1_pr_mg_mean, ylim = c(0, 1), type = "l",
  ylab = "pr(moralizing gods)", xlab = "social complexity",
  main = "original model: all 'NA' to '0'", xaxt = "n")
polygon(c(counter$Mean_c, rev(counter$Mean_c)),
  c(counter$m1_pr_mg_ub, rev(counter$m1_pr_mg_lb)),
  border = NA, col = col_alpha("black", 0.2))

axis(1, at = seq(0, 1, by = 0.2) - 0.5,
  labels = seq(0, 1, by = 0.2))

for (i in 1:nrow(d)) lines(c(d$Mean_c[i], d$Mean_c[i]),
  c(d$m1_pr_mg_lb[i], d$m1_pr_mg_ub[i]), col = col_alpha(d$na_col[i], 0.05))

points(d$Mean_c, d$m1_pr_mg_mean, pch = 16, col = d$na_col, cex = 0.6)

points(-0.45, 0.8, col = "firebrick", pch = 20)
text(-0.45, 0.8, "no MG data", col = "firebrick", pos = 4)

points(-0.45, 0.7, col = "dodgerblue", pch = 20)
text(-0.45, 0.7, "has MG data", col = "dodgerblue", pos = 4)

plot(counter$Mean_c, counter$m1_alt1_pr_mg_mean, ylim = c(0, 1), type = "l",
  ylab = "pr(moralizing gods)", xlab = "social complexity",
  main = "alternative 1: 96% 'NA' to '1'", xaxt = "n")
polygon(c(counter$Mean_c, rev(counter$Mean_c)),
  c(counter$m1_alt1_pr_mg_ub, rev(counter$m1_alt1_pr_mg_lb)),
  border = NA, col = col_alpha("black", 0.2))

axis(1, at = seq(0, 1, by = 0.2) - 0.5,
  labels = seq(0, 1, by = 0.2))

for (i in 1:nrow(d)) {
  lines(c(d$Mean_c[i], d$Mean_c[i]), c(d$m1_alt1_pr_mg_lb[i],
    d$m1_alt1_pr_mg_ub[i]), col = col_alpha(d$na_col[i], 0.05))
}

points(d$Mean_c, d$m1_alt1_pr_mg_mean, pch = 16, col = d$na_col, cex = 0.6)

points(0, 0.3, col = "firebrick", pch = 20)
text(0, 0.3, "no MG data", col = "firebrick", pos = 4)

points(0, 0.2, col = "dodgerblue", pch = 20)
text(0, 0.2, "has MG data", col = "dodgerblue", pos = 4)

plot(counter$Mean_c, counter$m1_alt2_pr_mg_mean, ylim = c(0, 1), type = "l",
  ylab = "pr(moralizing gods)", xlab = "social complexity",
  main = "alternative 2: 50% NA each to 0/1", xaxt = "n")
polygon(c(counter$Mean_c, rev(counter$Mean_c)),
  c(counter$m1_alt2_pr_mg_ub, rev(counter$m1_alt2_pr_mg_lb)),
  border = NA, col = col_alpha("black", 0.2))

axis(1, at = seq(0, 1, by = 0.2) - 0.5,
  labels = seq(0, 1, by = 0.2))

for (i in 1:nrow(d)) {
  lines(c(d$Mean_c[i], d$Mean_c[i]), c(d$m1_alt2_pr_mg_lb[i],
    d$m1_alt2_pr_mg_ub[i]), col = col_alpha(d$na_col[i], 0.05))
}

points(d$Mean_c, d$m1_alt2_pr_mg_mean, pch = 16, col = d$na_col, cex = 0.6)

points(0, 0.3, col = "firebrick", pch = 20)
text(0, 0.3, "no MG data", col = "firebrick", pos = 4)

points(0, 0.2, col = "dodgerblue", pch = 20)
text(0, 0.2, "has MG data", col = "dodgerblue", pos = 4)

dev.off()


print("plot predictions from m2")

counter$m2_pr_mg_mean <- NA
counter$m2_pr_mg_sd <- NA
counter$m2_pr_mg_lb <- NA
counter$m2_pr_mg_ub <- NA

for (i in 1:nrow(counter)) {
  logit_p <- m2_post$a +
    m2_post$b_sc * counter$Mean_c[i] +
    m2_post$b_sp * counter$Space[i] +
    m2_post$b_ph * counter$Phylogeny[i]
  pr_mg <- logistic(logit_p)
  counter$m2_pr_mg_mean[i] <- mean(pr_mg)
  counter$m2_pr_mg_sd[i] <- sd(pr_mg)
  counter$m2_pr_mg_lb[i] <- HPDI(pr_mg)[1]
  counter$m2_pr_mg_ub[i] <- HPDI(pr_mg)[2]
}

# these are the NGAs used in the m2 varying effects terms
NGAs_m2 <- c(
  "Big Island Hawaii",          "Cambodian Basin",
  "Central Java",               "Chuuk Islands",
  "Cuzco",                      "Deccan",
  "Garo Hills",                 "Ghanaian Coast",
  "Iceland",                    "Kachi Plain",
  "Kansai",                     "Kapuasi Basin",
  "Konya Plain",                "Latium",
  "Lena River Valley",          "Lowland Andes",
  "Middle Yellow River Valley", "Niger Inland Delta",
  "North Colombia",             "Orkhon Valley",
  "Oro PNG",                    "Paris Basin",
  "Sogdiana",                   "Susiana",
  "Upper Egypt",                "Yemeni Coastal Plain"
)

d$m2_pr_mg_mean <- NA
d$m2_pr_mg_sd <- NA
d$m2_pr_mg_lb <- NA
d$m2_pr_mg_ub <- NA
d$m2_hit50 <- NA

for (i in 1:nrow(d)) {
  if (d$NGA[i] %in% NGAs_m2) {
    nga <- match(d$NGA[i], NGAs_m2)
    nga_offset <- m2_post$a_nga[, nga]
  } else {
    nga_offset <- rnorm(length(m2_post$a_sigma), 0, m2_post$a_sigma)
  }
  pr_mg <- logistic(
    m2_post$a +
    nga_offset +
    m2_post$b_sc * d$Mean_c[i] +
    m2_post$b_sp * d$Space[i] +
    m2_post$b_ph * d$Phylogeny[i]
  )
  d$m2_pr_mg_mean[i] <- mean(pr_mg)
  d$m2_pr_mg_sd[i] <- sd(pr_mg)
  d$m2_pr_mg_lb[i] <- HPDI(pr_mg)[1]
  d$m2_pr_mg_ub[i] <- HPDI(pr_mg)[2]
  d$m2_hit50[i] <- mean(pr_mg > 0.5) > density_threshold
}



print("estimate MG first emergence from model posteriors")

NGAs_short <- c("Upper Egypt", "Susiana", "Konya Plain",
  "Middle Yellow River Valley", "Kachi Plain", "Sogdiana",
  "Latium", "Deccan", "Paris Basin", "Orkhon Valley",
  "Kansai", "Niger Inland Delta")

d$time_to_first_obs <- NA

for (i in 1:length(NGAs_short)) {
  nga_rows <- which(d$NGA == NGAs_short[i])
  nga_first_mg_row <- nga_rows[min(which(d$MG[nga_rows] == 1))]
  d$time_to_first_obs[nga_rows] <- d$Time[nga_rows] - d$Time[nga_first_mg_row]
}

# define the evidence threshold for "first appearance" analysis

nga_dat <- data.frame(NGA = NGAs_short)
nga_dat$m1_year_appear_50 <- NA
nga_dat$m2_year_appear_50 <- NA

for (i in 1:nrow(nga_dat)) {
  dn <- d[which(d$NGA == nga_dat$NGA[i] & d$time_to_first_obs < 0), ]
  if (any(dn$m1_hit50 == 1, na.rm = TRUE)) {
    hit_rows <- which(dn$m1_hit50 == 1)
    nga_dat$m1_year_appear_50[i] <- min(dn$time_to_first_obs[hit_rows])
  }
  if (any(dn$m2_hit50 == 1)) {
    nga_dat$m2_year_appear_50[i] <- min(dn$time_to_first_obs[dn$m2_hit50 == 1])
  }
}



print("visualize each NGA individually")

png("./temp/revised_EDfig1_m1.png", res = 300,
  height = 8, width = 10, units = "in")

par(mfrow = c(3, 4))

for (i in 1:nrow(nga_dat)) {
  dm <- d[which(d$NGA == nga_dat$NGA[i]), ]
  plot(dm$time_to_first_obs, dm$m1_pr_mg_mean, ylim = c(0, 1),
    xlim = c(-4000, 100), type = "l",
    xlab = "years before first apperance",
    ylab = "pr(moralizing gods present)",
    main = nga_dat$NGA[i])

  abline(v = nga_dat$m1_year_appear_50[i], lwd = 2)

  tar <- which(!is.na(dm$m1_pr_mg_lb))
  polygon(c(dm$time_to_first_obs[tar], rev(dm$time_to_first_obs[tar])),
    c(dm$m1_pr_mg_lb[tar], rev(dm$m1_pr_mg_ub[tar])),
    border = NA, col = col.alpha("firebrick", 0.2))
  abline(h = 0.5, col = "red")
  abline(v = 0, lty = 2)
}

dev.off()

png("./temp/revised_EDfig1_m2.png", res = 300,
  height = 8, width = 10, units = "in")

par(mfrow = c(3, 4))

for (i in 1:nrow(nga_dat)) {
  dm <- d[which(d$NGA == nga_dat$NGA[i]), ]
  plot(dm$time_to_first_obs, dm$m2_pr_mg_mean, ylim = c(0, 1),
    xlim = c(-4000, 100), type = "l", lwd = 2,
    xlab = "years before first apperance",
    ylab = "pr(moralizing gods present)",
    main = nga_dat$NGA[i])

  abline(v = nga_dat$m2_year_appear_50[i], lwd = 1, lty = 2)

  polygon(c(dm$time_to_first_obs, rev(dm$time_to_first_obs)),
    c(dm$m2_pr_mg_lb, rev(dm$m2_pr_mg_ub)),
    border = NA, col = col.alpha("firebrick", 0.2))
  abline(h = 0.5, col = "red")
}

dev.off()



print("aggregate first appearance estimates for revised fig2 from paper")

SCNorm <- read.csv("./input/SCNorm.csv", stringsAsFactors = FALSE)
data <- read.csv("./input/PrePostComparison.csv", stringsAsFactors = FALSE)

year_appear_50_mean <- mean(nga_dat$m2_year_appear_50, na.rm = TRUE)
year_appear_50_se <- sd(nga_dat$m2_year_appear_50, na.rm = TRUE) /
  sqrt(sum(!is.na(nga_dat$m2_year_appear_50)))

expect_true(abs(year_appear_50_mean - (-997)) < 50)
expect_true(abs(year_appear_50_se - (200)) < 50)

png("./temp/m2_predictions_fig2_combined.png", res = 300,
  height = 5, width = 10, units = "in")

par(mfrow = c(1, 2))

plot(counter$Mean_c, counter$m2_pr_mg_mean, ylim = c(0, 1), type = "l",
  ylab = "pr(moralizing gods)", xlab = "social complexity",
  main = "", xaxt = "n")
polygon(c(counter$Mean_c, rev(counter$Mean_c)),
  c(counter$m2_pr_mg_ub, rev(counter$m2_pr_mg_lb)),
  border = NA, col = col.alpha("dodgerblue", 0.2))

axis(1, at = seq(0, 1, by = 0.2) - 0.5, labels = seq(0, 1, by = 0.2))

tar <- which(d$MG_missing == 1)

points(d$Mean_c[tar], d$m2_pr_mg_mean[tar], pch = 16,
  col = col_alpha(d$nga_col[tar], 0.8), cex = 0.6)

points(d$Mean_c[tar], d$m1_pr_mg_mean[tar], pch = 16,
  col = col_alpha("gray", 0.8), cex = 0.6)

points(counter$Mean_c, counter$m1_pr_mg_mean,
  col = gray(0.45), type = "l", lty = 2)
text(0.4, 0.5, "original model", col = gray(0.45), srt = 63)

text(-0.48, 0.93, "A", cex = 2)

col1 <- rgb(0, 0, 0, max = 255, alpha = 50)

plot(SCNorm$x, SCNorm$Mean, type = "n", ylim = c(0, 1),
  xlim = c(-2000, 2000), xaxs = "i", yaxs = "i",
  xlab = "years from first MG observation", ylab = "social complexity")

# 'MG observed' period
rect(0, 0, 0 + mean(data$RangeMGAppear), 1,
  border = NA, col = col1)

lines(SCNorm$x, SCNorm$Mean, type = "l")
lines(SCNorm$x, SCNorm$Upper, type = "l", lty = "dotted")
lines(SCNorm$x, SCNorm$Lower, type = "l", lty = "dotted")

polygon(
  c(year_appear_50_mean + c(-1.96, 1.96) * year_appear_50_se,
  year_appear_50_mean + c(1.96, -1.96) * year_appear_50_se),
  c(0, 0, 1, 1),
  border = NA,
  col = col.alpha("firebrick", 0.3)
)

abline(v = year_appear_50_mean, col = "firebrick")

text(year_appear_50_mean, 0.7, "predicted MG\nemergence", srt = 90)
text(80, 0.3, "MG first recorded", srt = 90)

text(-1800, 0.9, "B", cex = 2)

dev.off()

#########

dir_init("./output")

files <- list.files("./temp", full.names = TRUE)
file.copy(files, "./output")
