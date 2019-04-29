
rm(list = ls())
source("../project_support.r")

dir_init("./temp")

######

print("Load regression dataset")

d <- read.csv("./input/RegrDat.csv", stringsAsFactors = FALSE)



print("Find when writing first appeared for each NGA")

writing <- d %>%
  filter(Writing == 1) %>%
  group_by(NGA) %>%
  filter(row_number()==1) %>%
  dplyr::select(NGA, Time, Writing) %>%
  dplyr::rename(`Writing First Recorded` = Time)

# Find when moralising gods first appeared for each NGA
MG <- d %>%
  filter(MoralisingGods == 1) %>%
  group_by(NGA) %>%
  filter(row_number()==1) %>%
  dplyr::select(NGA, Time, MoralisingGods) %>%
  dplyr::rename(`Moralizing Gods First Recorded` = Time)

# join writing and moralising gods data
MGwriting <- left_join(writing, MG)

# remove Valley of Oaxaca which is always NA for MGs
MGwriting <- MGwriting %>%
  filter(NGA != "Valley of Oaxaca") %>%
  # calculate time difference between writing and MGs first appearing (negative means MGs appeared first)
   mutate(TimeDifference = `Moralizing Gods First Recorded` -  `Writing First Recorded`) %>%
  # add whether NGA was used in analysis or not
   mutate(`Social complexity data available \nboth before and after the appearance \nof MG` = if_else(NGA == "Deccan" | NGA == "Kachi Plain" | NGA == "Kansai" | NGA == "Konya Plain" | NGA == "Latium" | NGA == "Middle Yellow River Valley" | NGA == "Niger Inland Delta" | NGA == "Orkhon Valley" | NGA == "Paris Basin" | NGA == "Sogdiana" | NGA == "Susiana" | NGA == "Upper Egypt", "Yes", "No"))

# extract just NGAs used in analysis
MGwritinganl <- MGwriting[MGwriting$`Social complexity data available 
both before and after the appearance 
of MG` == "Yes",]

# calculate the mean difference between the appearance of writing and moralising gods of just NGAs included in analysis
mean(MGwritinganl$TimeDifference)

# calculate the median difference between the appearance of writing and moralising gods of just NGAs included in analysis
median(MGwritinganl$TimeDifference)

# calculate the mean difference between the appearance of writing and moralising gods across all NGAs  
mean(MGwriting$TimeDifference)

# calculate the median difference between the appearance of writing and moralising gods across all NGAs  
median(MGwriting$TimeDifference)

# set the time difference between writing and moralizing gods to 0 for Kachi Plain
MGwritingfilt <- MGwriting %>%
  mutate(TimeDifference = if_else(NGA == "Kachi Plain", 0, as.numeric(TimeDifference)))
# calculate the mean difference between the appearance of writing and moralising gods across all NGAs with the difference in Kachi Plain set to 0
mean(MGwritingfilt$TimeDifference)

# set the time difference between writing and moralizing gods to 0 for Kachi Plain for subset of NGAs used in analysis
MGwritinganlfilt <- MGwritinganl %>%
  mutate(TimeDifference = if_else(NGA == "Kachi Plain", 0, as.numeric(TimeDifference)))
# calculate the mean difference between the appearance of writing and moralising gods NGAs used in analysis with the difference in Kachi Plain set to 0
mean(MGwritinganlfilt$TimeDifference)

# save output
ggplot(MGwriting, aes(x = `Writing First Recorded`, y = `Moralizing Gods First Recorded`, colour = `Social complexity data available \nboth before and after the appearance \nof MG`)) +
  geom_abline(slope = 1, intercept = 0, color = 'gray50', size = 1) +
  geom_abline(slope = 1, intercept = 100, color = 'gray60', linetype = "dashed") +
  geom_abline(slope = 1, intercept = -100, color = 'gray60', linetype = "dashed") +
  geom_point(size = 3) +
  geom_text(aes(x = `Writing First Recorded`, y = `Moralizing Gods First Recorded`, label=ifelse(NGA == "Susiana" | NGA == "Kachi Plain", NGA,'')), vjust = -1.25, inherit.aes = FALSE, size = 4, color = "#0072B2") +
  labs(color = "Social complexity data available \nboth before and after the appearance \nof MG") +
  xlab("Writing First Recorded") +
  ylab("Moralizing Gods First Recorded") +
  scale_color_manual(values=c("#E69F00", "#0072B2")) + 
  scale_y_continuous(limits = c(-4000,2000)) +
  scale_x_continuous(limits = c(-4000,2000)) +
  theme_bw() +
  theme(
    legend.position = c(0.70, 0.2),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 13),
    axis.text = element_text(colour = "black", size=11),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    axis.line.x = element_line(color="black", size = 0.5),
    axis.line.y = element_line(color="black", size = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank())

ggsave("./temp/WritingMG.png",width = 6,height = 5.5, dpi = 300)


print("show transition to MG = 1 in 12 NGAs")

d$MG_original <- as.character(d$MG)
d$MG_original[d$MG_missing == 1] <- "NA"
d$text_color <- ifelse(d$MG_missing == 1, gray(0.3), "black")
d$cell_color <- case_when(
  d$Writing == 1 & d$MG_original %in% c("1", "NA") ~ "#d6e6a5",
  d$Writing == 1 & d$MG_original == "0" ~ "#FFA8A0",
  d$Writing == 0 ~ "white"
)

# calculate "time of first appearance" and subtract off

NGAs_short <- c("Deccan", "Kachi Plain", "Kansai", "Konya Plain",
  "Latium", "Middle Yellow River Valley", "Niger Inland Delta", 
  "Orkhon Valley", "Paris Basin", "Sogdiana", "Susiana", "Upper Egypt")

d$time_to_first_obs <- NA

for (i in 1:length(NGAs_short)) {
  nga_rows <- which(d$NGA == NGAs_short[i])
  nga_first_mg_row <- nga_rows[min(which(d$MG[nga_rows] == 1))]
  d$time_to_first_obs[nga_rows] <- d$Time[nga_rows] - d$Time[nga_first_mg_row]
}

png("./temp/missingness_table.png", res = 300, units = "in", height = 5.5, width = 11)

par(mar = c(5.1, 12, 4.1, 2.1))

plot(1, 1, type = "n", frame.plot = FALSE,
  xlim = c(-1010, 100), ylim = c(1, 12.5), axes = FALSE, ann = FALSE)
axis(2, at = 12:1, labels = NGAs_short, las = 1)

years <- seq(-1000, 100, by = 100)

for (i in 1:length(NGAs_short)) {

  my_rows <- which(d$NGA == NGAs_short[i] &
    d$time_to_first_obs <= 100 & d$time_to_first_obs >= -1000)

  for (j in 1:length(my_rows)) {
    rect(
      d$time_to_first_obs[my_rows[j]] - 50,
      (13 - i) - 0.5,
      d$time_to_first_obs[my_rows[j]] + 50,
      (13 - i) + 0.5,
      col = d$cell_color[my_rows[j]], border = NA
    )
  }

  text(d$time_to_first_obs[my_rows], 13 - i, labels = d$MG_original[my_rows],
    col = d$text_color[my_rows])

}

abline(h = 1:13 - 0.5, col = "gray")

axis(3, years)
mtext(text = "years until first appearance of moralizing god", line = 2.3)

dev.off()



print("show MG missingness pattern by SC")

d$MG_known <- 1 - d$MG_missing

d$MG_code <- case_when(
  d$MG_known == 1 & d$MG == 0 ~ 1,
  d$MG_known == 1 & d$MG == 1 ~ 2,
  d$MG_known == 0 ~ 3
)

# a barbell plot

dm <- d[, c("NGA", "OriginalNGA", "Writing", "PolPop", "MG", "Mean", "Lag1", "Lag2", "MG_known")]

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



print("calculate descriptives")

# calc missingness patterns

m1 <- lm(PolPop ~ MG_known, data = dm)

pop1 <- 10^(coef(m1)[1] + coef(m1)[2] * 1) # 3.9 million people
pop0 <- 10^(coef(m1)[1] + coef(m1)[2] * 0) # 7756 ppl

expect_true(abs(pop1 - 3891008) < 10000)
expect_true(abs(pop0 - 7756) < 100)

m2 <- glm(MG_known ~ Writing, data = dm, family = "binomial")

table(d$MG_known, d$Writing)

pr_missing1 <- logistic(coef(m2)[1] + coef(m2)[2]) # 0.77
pr_missing0 <- logistic(coef(m2)[1]) # 0.03

m3 <- glm(Writing ~ MG_known, data = dm, family = "binomial")

pr_read1 <- logistic(coef(m3)[1] + coef(m3)[2]) # 0.95
pr_read0 <- logistic(coef(m3)[1]) # 0.17

exp(coef(m2)[2]) # OR: 99

m3 <- lm(PolPop ~ Mean, data = dm)

10^(coef(m1)[1] + coef(m1)[2] * 0.2) # 27k
10^(coef(m1)[1] + coef(m1)[2] * 0.4) # 93k
10^(coef(m1)[1] + coef(m1)[2] * 0.6) # 323k

length(unique(dm$NGA[which(dm$Mean < 0.4)])) # 22

sum(dm$MG_known == 1 & dm$Mean < 0.4) # 8 rows

sum(dm$MG[dm$MG_known == 1 & dm$Mean < 0.4]) # 3 rows


print("make a population vs social complexity figure with missigness shading")

expect_equal(nrow(dm), 801)

dm$logPop <- dm$PolPop

dm$na_col <- ifelse(dm$MG_known == 0, "gray", "black")

cor(dm$Mean, dm$logPop) # 0.94

png("./temp/sc_pop.png", res = 300, units = "in", height = 5, width = 5)

plot(dm$Mean, dm$logPop, col = col_alpha(dm$na_col, 0.7), pch = 16,
  xlim = c(0, 1), xlab = "social complexity",
  ylab = "log10 population", ylim = c(1.5, 8.5), frame.plot = FALSE)

abline(lm(logPop ~ Mean, data = dm), lty = 2)

rect(0.00, 6, 0.4, 7.5)
points(0.05, 7, col = "black", pch = 20)
text(0.05, 7, label = "has MG data", col = "black", pos = 4)
points(0.05, 6.5, col = "gray", pch = 20)
text(0.05, 6.5, label = "MG unknown", col = gray(0.4), pos = 4)

dev.off()

#########

dir_init("./output")

files <- list.files("./temp", full.names = TRUE)
file.copy(files, "./output")
