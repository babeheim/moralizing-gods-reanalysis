
rm(list = ls())
source("../project_support.r")

dir_init("./temp")

############

# load the data

d <- read.csv("./input/RegrDat.csv", stringsAsFactors = FALSE)

# add centering
d$Mean_c <- d$Mean - 0.5

# define models for analysis

original_model <- alist(
  MG <- dbinom(1, p),
  logit(p) <- a +
    b_sc * Mean_c +
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

# rerun the original analysis with NA's as 0's

dm <- d[, c("MG", "Mean_c", "Lag1", "Lag2", "Phylogeny", "Space")]

drop <- which(is.na(dm$Lag1) | is.na(dm$Lag2))
dm <- dm[-drop, ]

# 801 in the original analysis
expect_equal(nrow(dm), 801)

m1 <- map2stan(original_model, data = dm)

save(m1, file = "./temp/m1.rdata")

# these should be comparable to the original results

# original model but re-assign NAs according to another missingness rule:

dm <- d[, c("NGA", "MG", "Mean_c", "Lag1", "Lag2", "Phylogeny", "Space")]

imputation_prob <- 311 / 323 # occurance of 1's in known data

dm$MG[dm$MG_missing == 1] <- rbinom(sum(dm$MG_missing == 1), 1, prob = imputation_prob)

# recalculate lag terms accordingly

NGAs <- sort(unique(dm$NGA))

for (i in 1:length(NGAs)) {
  my_rows <- which(dm$NGA == NGAs[i])
  if (length(my_rows) > 1) {
    for (j in 2:length(my_rows)) {
      dm$Lag1[my_rows[j]] <- dm$MG[my_rows[j - 1]]
    }
  }
  if (length(my_rows) > 2) {
    for (j in 3:length(my_rows)) {
      dm$Lag2[my_rows[j]] <- dm$MG[my_rows[j - 2]]
    }
  }
}

drop <- which(is.na(dm$Lag1) | is.na(dm$Lag2))
if (length(drop) > 0) dm <- dm[-drop, ]

m1_alt1 <- map2stan(original_model, data = dm)

save(m1_alt1, file = "./temp/m1_alt1.rdata")

# original model but re-assign NAs according to another missingness rule:

dm <- d[, c("NGA", "MG", "Mean_c", "Lag1", "Lag2", "Phylogeny", "Space")]

imputation_prob <- 0.5 # principle of indifference

dm$MG[dm$MG_missing == 1] <- rbinom(sum(dm$MG_missing == 1), 1, prob = imputation_prob)

# recalculate lag terms accordingly

NGAs <- sort(unique(dm$NGA))

for (i in 1:length(NGAs)) {
  my_rows <- which(dm$NGA == NGAs[i])
  if (length(my_rows) > 1) {
    for (j in 2:length(my_rows)) {
      dm$Lag1[my_rows[j]] <- dm$MG[my_rows[j - 1]]
    }
  }
  if (length(my_rows) > 2) {
    for (j in 3:length(my_rows)) {
      dm$Lag2[my_rows[j]] <- dm$MG[my_rows[j - 2]]
    }
  }
}

drop <- which(is.na(dm$Lag1) | is.na(dm$Lag2))
if (length(drop) > 0) dm <- dm[-drop, ]

m1_alt2 <- map2stan(original_model, data = dm)

save(m1_alt2, file = "./temp/m1_alt2.rdata")

# now analyze complete cases only; drop the NA's!

dm <- d[, c("NGA", "MG", "Mean_c", "Phylogeny", "Space", "MG_missing")]

# drop all "unknown" outcomes
drop <- which(dm$MG_missing == 1)
dm <- dm[-drop, ]

expect_equal(nrow(dm), 336)
# 311 cases are left from the 801 originally
# but 25 more because we don't have to drop Lag1 = NA cases

# 26 NGAs in the reduced set
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

dm$nga <- match(dm$NGA, NGAs)

m2 <- map2stan(revised_model, data = dm)

save(m2, file = "./temp/m2.rdata")

# now the revised model on the original data with NA's as 0's

dm <- d[, c("NGA", "MG", "Mean_c", "Phylogeny", "Space")]

NGAs <- sort(unique(dm$NGA))
dm$nga <- match(dm$NGA, NGAs)

m2_alt1 <- map2stan(revised_model, data = dm)

save(m2_alt1, file = "./temp/m2_alt1.rdata")

############

dir_init("./output")

files <- list.files("./temp", full.names = TRUE)
file.copy(files, "./output")
