
rm(list = ls())
source("../project_support.r")

dir_init("./temp")

d <- read.csv("./input/RegrDat.csv", stringsAsFactors = FALSE)

dm <- d[, c("MG", "Mean", "Lag1", "Lag2", "Phylogeny", "Space")]

# original analysis dropped first and second observations for each NGA
# because lag terms had to be missing
drop <- which(is.na(dm$Lag1) | is.na(dm$Lag2))
dm <- dm[-drop, ]

nrow(dm) == 801

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

m1 <- map2stan(original_model, data = dm)

save(m1, file = "./temp/m1.rdata")

# these results should be comparable to published model



dm <- d[, c("NGA", "MG", "Mean", "Phylogeny", "Space", "MG_missing", "Lag1", "Lag2")]

# drop all "unknown" outcomes
drop <- which(dm$MG_missing == 1)
dm <- dm[-drop, ]

nrow(dm) == 336

# if we also dropped cases where Lag1 is NA we have
# 311 rows, because the 25 NGAs each have a first observation
# we get to keep those if no lag terms

# center social complexity at 0.5
dm$Mean_c <- dm$Mean - 0.5

# create an index for NGA varying effects

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

save(m2, file = "./temp/m2.rdata")

dir_init("./output")

file.copy(c("./temp/m1.rdata", "./temp/m2.rdata"), "./output")
