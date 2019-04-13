
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

10^(coef(m1)[1] + coef(m1)[2] * 1) # 2.9 million people
10^(coef(m1)[1] + coef(m1)[2] * 0) # 7000 ppl

m2 <- glm(Writing ~ MG_known, data = d, family = "binomial")

logistic(coef(m2)[1] + coef(m2)[2]) # 0.92
logistic(coef(m2)[1]) # 0.16
exp(coef(m2)[2]) # OR: 64

m3 <- lm(d$PolPop ~ d$Mean, data = d)

10^(coef(m1)[1] + coef(m1)[2] * 0.2) # 23k
10^(coef(m1)[1] + coef(m1)[2] * 0.4) # 78k
10^(coef(m1)[1] + coef(m1)[2] * 0.6) # 260k


# visualize the missingness problem

d$MG[d$MG_missing == 1] <- NA

NGAs <- sort(unique(d$NGA))
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

png("./temp/Thirty.png", res = 300, units = "in", height = 12, width = 16)

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

  axis(1, tck = 0.02, at = seq(min(d$Time[my_rows]),
    max(d$Time[my_rows]), by = 100), labels = FALSE)

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
  keep[i] <- any(d$MG[my_rows] == 1) &
    any(is.na(d$MG[my_rows] | d$MG[my_rows] == 0))
}

NGA_keep <- NGAs[which(keep)]



png("./temp/Thirteen.png", res = 300, units = "in", height = 12, width = 8)

par(mfrow = c(5, 3))

plot(c(0, 10), c(0, 10), type = "n", frame.plot = FALSE, axes = FALSE,
 xlab = "", ylab = "", main = "Legend")

symbol_line <- 0.5
symbol_anchor <- 8

lines(symbol_line + c(-0.3, 0.3), c(symbol_anchor, symbol_anchor))
points(symbol_line, symbol_anchor, pch = 21, cex = 1, bg = "black")
text(symbol_line + 0.5, symbol_anchor, "moralizing gods present", pos = 4)

lines(symbol_line + c(-0.3, 0.3), c(symbol_anchor, symbol_anchor) - 2)
points(symbol_line, symbol_anchor - 2, pch = 21, cex = 1, bg = "white")
text(symbol_line + 0.5, symbol_anchor - 2, "moralizing gods absent", pos = 4)

lines(symbol_line + c(-0.3, 0.3), c(symbol_anchor, symbol_anchor) - 4)
text(symbol_line + 0.5, symbol_anchor - 4,
  "no data on moralizing gods", pos = 4)

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

  axis(1, tck = 0.02, at = seq(min(d$Time[my_rows]),
    max(d$Time[my_rows]), by = 100), labels = FALSE)

  my_points <- which(d$NGA == NGA_keep[i] & d$MG == 1)
  points(d$Time[my_points], d$PolPop[my_points],
    bg = d$mg_col[my_points], col = d$mg_outline[my_points], pch = 21)
  my_points <- which(d$NGA == NGA_keep[i] & d$MG == 0)
  points(d$Time[my_points], d$PolPop[my_points],
    bg = d$mg_col[my_points], col = d$mg_outline[my_points], pch = 21)

}

dev.off()

# subset to a representative 5 for discussion

NGA_keep <- c("North Colombia", "Big Island Hawaii",
  "Valley of Oaxaca", "Orkhon Valley", "Latium")

png("./temp/Five.png", res = 300, units = "in", height = 4.8, width = 8)

par(mfrow = c(2, 3))

plot(c(0, 10), c(0, 10), type = "n", frame.plot = FALSE, axes = FALSE,
 xlab = "", ylab = "", main = "Legend")

symbol_line <- 0.5
symbol_anchor <- 8

lines(symbol_line + c(-0.3, 0.3), c(symbol_anchor, symbol_anchor))
points(symbol_line, symbol_anchor, pch = 21, cex = 1, bg = "black")
text(symbol_line + 0.5, symbol_anchor, "moralizing gods present", pos = 4)

lines(symbol_line + c(-0.3, 0.3), c(symbol_anchor, symbol_anchor) - 2)
points(symbol_line, symbol_anchor - 2, pch = 21, cex = 1, bg = "white")
text(symbol_line + 0.5, symbol_anchor - 2, "moralizing gods absent", pos = 4)

lines(symbol_line + c(-0.3, 0.3), c(symbol_anchor, symbol_anchor) - 4)
text(symbol_line + 0.5, symbol_anchor - 4,
  "no data on moralizing gods", pos = 4)
points(symbol_line, symbol_anchor - 4, col = col.alpha("gray", 0.8), pch = 1)

for (i in 1:length(NGA_keep)) {

  my_rows <- which(d$NGA == NGA_keep[i])
  plot(d$Time[my_rows], d$PolPop[my_rows], ylim = c(1, 8.5),
    type = "n", main = NGA_keep[i], frame.plot = FALSE,
    ylab = "log10 population", xlab = "years common era",
    tck = 0.05, las = 1)

  points(d$Time[my_rows], d$PolPop[my_rows],
    col = nga_col[i], type = "l")
  if (length(my_rows) == 1) points(d$Time[my_rows], d$PolPop[my_rows],
    pch = 20, cex = 0.2, col = nga_col[i])

  axis(1, tck = 0.02, at = seq(min(d$Time[my_rows]),
    max(d$Time[my_rows]), by = 100), labels = FALSE)

  my_points <- which(d$NGA == NGA_keep[i] & d$MG == 1)
  points(d$Time[my_points], d$PolPop[my_points],
    bg = d$mg_col[my_points], col = d$mg_outline[my_points], pch = 21)
  my_points <- which(d$NGA == NGA_keep[i] & d$MG == 0)
  points(d$Time[my_points], d$PolPop[my_points],
    bg = d$mg_col[my_points], col = d$mg_outline[my_points], pch = 21)

}

dev.off()


# make a population vs social complexity figure with missigness shading

png("./temp/sc_pop.png", res = 300, units = "in", height = 5, width = 5)

d$logPop <- d$PolPop
d$Pop <- 10^d$logPop

d$na_col <- ifelse(d$MG_missing == 1, "gray", "black")

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
