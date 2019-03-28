
rm(list = ls())
source("../project_support.r")

dir_init("./temp")


# explore the missingness patterns

d <- read.csv("./input/RegrDat.csv", stringsAsFactors = FALSE)

d$MG_known <- 1 - d$MG_missing

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



png("./temp/Thirteen.png", res = 300, units = "in", height = 12, width = 8)

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
# points(symbol_line, symbol_anchor - 4, col = col.alpha("gray", 0.8), pch = 1)


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

  # my_points <- which(d$NGA == NGA_keep[i] & is.na(d$MG))
  # points(d$Time[my_points], d$PolPop[my_points],
  #   col = col.alpha("gray", 0.8), pch = 1)
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

png("./temp/Five.png", res = 300, units = "in", height = 4.8, width = 8)

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

  # my_points <- which(d$NGA == NGA_keep[i] & is.na(d$MG))
  # points(d$Time[my_points], d$PolPop[my_points],
  #   col = col.alpha("gray", 0.8), pch = 1)
  my_points <- which(d$NGA == NGA_keep[i] & d$MG == 1)
  points(d$Time[my_points], d$PolPop[my_points],
    bg = d$mg_col[my_points], col = d$mg_outline[my_points], pch = 21)
  my_points <- which(d$NGA == NGA_keep[i] & d$MG == 0)
  points(d$Time[my_points], d$PolPop[my_points],
    bg = d$mg_col[my_points], col = d$mg_outline[my_points], pch = 21)

}

dev.off()

# this should be its own step...mb in the script that creates regrdat



