
rm(list = ls())
source("../project_support.r")

dir_init("./temp")

############

print("load the regression data frame")

d <- read.csv("./input/RegrDat.csv", stringsAsFactors = FALSE)

# add centering
d$Mean_c <- d$Mean - 0.5



print("compile the stan models")

m1_bin <- stan_model(file = "./stan/original.stan")
m2_bin <- stan_model(file = "./stan/revised.stan")



print("fit the original model, with NA's as 0's")

dm <- d[, c("MG", "Mean_c", "Lag1", "Lag2", "Phylogeny", "Space")]

drop <- which(is.na(dm$Lag1) | is.na(dm$Lag2))
dm <- dm[-drop, ]

# 801 in the original analysis
expect_equal(nrow(dm), 801)

dm <- as.list(dm)
dm$N <- length(dm$MG)

m1 <- sampling(m1_bin, data = dm)

save(m1, file = "./temp/m1.rdata")



print("re-impute NAs at observed rate and fit again")

dm <- d[, c("NGA", "MG", "Mean_c", "Lag1", "Lag2",
  "Phylogeny", "Space", "MG_missing")]

drop <- which(is.na(dm$Lag1) | is.na(dm$Lag2))
if (length(drop) > 0) dm <- dm[-drop, ]

mg_known_present <- sum(dm$MG == 1)
mg_known_absent <- sum(dm$MG == 0 & dm$MG_missing == 0)
mg_known <- sum(dm$MG_missing == 0)
mg_missing <- sum(dm$MG_missing == 1)
missing_rows <- which(dm$MG_missing == 1)

expect_equal(mg_known_present, 299)
expect_equal(mg_known, 311)

imputation_prob <- mg_known_present / mg_known # occurance of 1's in known data

# re-impute missing MG data

set.seed(1234)
dm$MG[missing_rows] <- rbinom(mg_missing, 1, prob = imputation_prob)

# recalculate lag terms accordingly

NGAs <- sort(unique(dm$NGA))

for (i in 1:length(NGAs)) {
  nga_rows <- which(dm$NGA == NGAs[i])
  if (length(nga_rows) > 1) {
    for (j in 2:length(nga_rows)) {
      dm$Lag1[nga_rows[j]] <- dm$MG[nga_rows[j - 1]]
    }
  }
  if (length(nga_rows) > 2) {
    for (j in 3:length(nga_rows)) {
      dm$Lag2[nga_rows[j]] <- dm$MG[nga_rows[j - 2]]
    }
  }
}

dm <- as.list(dm)
dm$N <- length(dm$MG)

m1_alt1 <- sampling(m1_bin, data = dm)

save(m1_alt1, file = "./temp/m1_alt1.rdata")



print("re-impute NAs using coin flip and fit again")

dm <- d[, c("NGA", "MG", "Mean_c", "Lag1", "Lag2",
  "Phylogeny", "Space", "MG_missing")]

drop <- which(is.na(dm$Lag1) | is.na(dm$Lag2))
if (length(drop) > 0) dm <- dm[-drop, ]

mg_missing <- sum(dm$MG_missing == 1)

imputation_prob <- 0.5 # principle of indifference

# re-impute missing MG data

set.seed(1234)
dm$MG[missing_rows] <- rbinom(mg_missing, 1, prob = imputation_prob)

# recalculate lag terms accordingly

NGAs <- sort(unique(dm$NGA))

for (i in 1:length(NGAs)) {
  nga_rows <- which(dm$NGA == NGAs[i])
  if (length(nga_rows) > 1) {
    for (j in 2:length(nga_rows)) {
      dm$Lag1[nga_rows[j]] <- dm$MG[nga_rows[j - 1]]
    }
  }
  if (length(nga_rows) > 2) {
    for (j in 3:length(nga_rows)) {
      dm$Lag2[nga_rows[j]] <- dm$MG[nga_rows[j - 2]]
    }
  }
}

dm <- as.list(dm)
dm$N <- length(dm$MG)

m1_alt2 <- sampling(m1_bin, data = dm)

save(m1_alt2, file = "./temp/m1_alt2.rdata")



print("analyze complete cases only with revised model")

dm <- d[, c("NGA", "MG", "Mean_c", "Phylogeny", "Space", "MG_missing")]

# drop all "unknown" outcomes
drop <- which(dm$MG_missing == 1)
dm <- dm[-drop, ]

expect_equal(nrow(dm), 336)
# 311 cases are left from the 801 originally
# but 25 more because we don't have to drop Lag1 = NA cases

# 26 NGAs in the reduced set
NGAs <- c(
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

dm$nga <- match(dm$NGA, NGAs)

dm <- as.list(dm)
dm$N <- length(dm$MG)
dm$N_nga <- length(unique(dm$nga))

m2 <- sampling(m2_bin, data = dm)

save(m2, file = "./temp/m2.rdata")

dm <- d[, c("MG", "Mean_c", "Phylogeny", "Space")]



print("fit revised model to original data")

dm <- d[, c("NGA", "MG", "Mean_c", "Lag1", "Lag2", "Phylogeny", "Space")]

drop <- which(is.na(dm$Lag1) | is.na(dm$Lag2))
if (length(drop) > 0) dm <- dm[-drop, ]

expect_equal(nrow(dm), 801)

NGAs <- sort(unique(dm$NGA))
dm$nga <- match(dm$NGA, NGAs)

dm <- as.list(dm)
dm$N <- length(dm$MG)
dm$N_nga <- length(unique(dm$nga))

expect_equal(length(unique(dm$nga)), 28)

m2_alt1 <- sampling(m2_bin, data = dm)

save(m2_alt1, file = "./temp/m2_alt1.rdata")



print("fit revised model to writing-only data per Savage, et al.'s reply")

dm <- d[, c("NGA", "MG", "Mean_c", "Phylogeny", "Space", "Writing", "MG_missing")]

# drop all "unknown" outcomes - but only for populations without writing
drop <- which(dm$MG_missing == 1 & dm$Writing == 0)
dm <- dm[-drop, ]

expect_equal(nrow(dm), 422)

# 27 NGAs in the reduced set - now including Valley of Oaxaca
NGAs <- c(
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
  "Upper Egypt",                "Yemeni Coastal Plain",
  "Valley of Oaxaca"
)

dm$nga <- match(dm$NGA, NGAs)

dm <- as.list(dm)
dm$N <- length(dm$MG)
dm$N_nga <- length(unique(dm$nga))

m2_alt2 <- sampling(m2_bin, data = dm)

save(m2_alt2, file = "./temp/m2_alt2.rdata")


############

dir_init("./output")

files <- list.files("./temp", full.names = TRUE)
file.copy(files, "./output")
