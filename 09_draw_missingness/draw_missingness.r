
rm(list = ls())
source("../project_support.r")

dir_init("./temp")

# explore the missingness patterns

d <- read.csv("./input/RegrDat.csv", stringsAsFactors = FALSE)

d$MG_known <- 1 - d$MG_missing

d$MG_code <- case_when(
  d$MG_known == 1 & d$MG == 0 ~ 1,
  d$MG_known == 1 & d$MG == 1 ~ 2,
  d$MG_known == 0 ~ 3
)

# a barbell plot

dm <- d[, c("MG", "Mean", "Lag1", "Lag2", "MG_known")]

# drop cases to create publication dataset
drop <- which(is.na(dm$Lag1) | is.na(dm$Lag2))
dm <- dm[-drop, ]

n_known_present <- sum(dm$MG == 1)
n_known_absent <- sum(dm$MG == 0 & dm$MG_known == 1)
n_unknown <- sum(dm$MG_known == 0)

# 801 in the original analysis
expect_equal(nrow(dm), 801)
expect_equal(n_known_present, 299)
expect_equal(n_known_absent, 12)
expect_equal(n_unknown, 490)

expect_true(abs(cor(dm$MG_known, dm$MG) - 0.97) < 0.01)

png("./temp/barbell.png", res = 300, units = "in", height = 5, width = 6)

plot(1, 1, type = "n", xlim = c(0.5, 3.5), ylim = c(0, 1),
  frame.plot = FALSE, ylab = "social complexity", xaxt = "n",
  xlab = "'moralizing gods' status")

abline(h = seq(0, 1, 0.2), col = "gray")

axis(1, at = c(1, 2, 3), col = NA, col.ticks = NA,
  labels = c(paste0("\nabsent\nn=", n_known_absent),
  paste0("\npresent\nn=", n_known_present), paste0("\nunknown\nn=", n_unknown)))

d$MG_col <- ifelse(d$MG == 1, "firebrick", "dodgerblue")

boxplot(Mean ~ MG_code, data = d, add = TRUE, ann = FALSE,
  frame.plot = FALSE, axes = FALSE, outline = FALSE, col = "white")

tar <- which(d$MG_known == 1 & d$MG == 0)
points(1 + density_offset(d$Mean[tar], scale = 0.09), d$Mean[tar], pch = 16,
  col = col_alpha(d$MG_col[tar], 0.8))
tar <- which(d$MG_known == 1 & d$MG == 1)
points(2 + density_offset(d$Mean[tar], scale = 0.15), d$Mean[tar], pch = 16,
  col = col_alpha(d$MG_col[tar], 0.8))
tar <- which(d$MG_known == 0)
points(3 + density_offset(d$Mean[tar], scale = 0.15), d$Mean[tar], pch = 16,
  col = col.alpha("gray", 0.5))

dev.off()


# calc missingness patterns

m1 <- lm(PolPop ~ MG_known, data = d)

pop1 <- 10^(coef(m1)[1] + coef(m1)[2] * 1) # 2.9 million people
pop0 <- 10^(coef(m1)[1] + coef(m1)[2] * 0) # 7000 ppl

expect_true(abs(pop1 - 2921630) < 10000)
expect_true(abs(pop0 - 6954) < 100)

m2 <- glm(Writing ~ MG_known, data = d, family = "binomial")

pr_read1 <- logistic(coef(m2)[1] + coef(m2)[2]) # 0.92
pr_read0 <- logistic(coef(m2)[1]) # 0.16

expect_true(abs(pr_read1 - 0.9256) < 0.01)
expect_true(abs(pr_read0 - 0.1629) < 0.01)

exp(coef(m2)[2]) # OR: 64

m3 <- lm(d$PolPop ~ d$Mean, data = d)

10^(coef(m1)[1] + coef(m1)[2] * 0.2) # 23k
10^(coef(m1)[1] + coef(m1)[2] * 0.4) # 78k
10^(coef(m1)[1] + coef(m1)[2] * 0.6) # 260k

# make a population vs social complexity figure with missigness shading

expect_equal(nrow(d), 864)

d$logPop <- d$PolPop
d$Pop <- 10^d$logPop

d$na_col <- ifelse(d$MG_missing == 1, "gray", "black")

cor(d$Mean, d$logPop) # 0.94

png("./temp/sc_pop.png", res = 300, units = "in", height = 5, width = 5)

# plot

plot(d$Mean, d$logPop, col = col_alpha(d$na_col, 0.8), pch = 16,
  xlim = c(0, 1), xlab = "social complexity",
  ylab = "log10 population", ylim = c(1.5, 8.5), frame.plot = FALSE)

abline(lm(logPop ~ Mean, data = d), lty = 2)

rect(0.05, 6, 0.35, 7.5)
points(0.1, 7, col = "black", pch = 20)
text(0.1, 7, label = "observed", col = "black", pos = 4)
points(0.1, 6.5, col = "gray", pch = 20)
text(0.1, 6.5, label = "missing", col = gray(0.4), pos = 4)

dev.off()

#########

dir_init("./output")

files <- list.files("./temp", full.names = TRUE)
file.copy(files, "./output")
