#####################################
# Financial Econometrics Term Paper #
#####################################

# Thore Petersen
# October 2019

library(data.table)
library(rugarch)
library(sn) # skew distributions
library(fitdistrplus) # extends MASS::fitdistr(); probably overkill, we'll see
library(ggplot2)
library(copula)


# Load data ---------------------------------------------------------------

load("./data/tmp/02_tmp.RData")


# DEU -- FRA pair ---------------------------------------------------------

# extract pairwise data for DEU and FRA, reshape/convert to numeric matrix
dt_DF <- rbindlist(list(dt["DEU"], dt["FRA"]))
c_dt_DF <- dcast(dt_DF, Date ~ iso3c, value.var = "eqtl")
c_dt_DF[, Date := NULL]
c_dt_DF <- as.matrix(c_dt_DF)

test <- fitCopula(copula = tCopula(), data = c_dt_DF)
