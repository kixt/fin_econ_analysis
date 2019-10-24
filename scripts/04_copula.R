#####################################
# Financial Econometrics Term Paper #
#####################################

# Thore Petersen
# October 2019

library(data.table)
library(ggplot2)
library(copula)


# Load data ---------------------------------------------------------------

load("./data/tmp/03_tmp.RData")


# Empirical quantiles -----------------------------------------------------

# ECDF of residuals
# copula::pobs() returns ecdf values scaled by n/(n+1), st. everything in unit cube
dt[, eqntl := pobs(res), by = iso3c]

# create data.table of empirical quantiles to check dependence
eqntls <- dcast(dt, Date ~ iso3c, value.var = "eqntl")


# DEU -- FRA pair ---------------------------------------------------------

# characteristic pattern for DEU--FRA, almost linear, with clustering in corners,
# i.e. strong tail dependence
ggplot(eqntls, aes(x = FRA, y = DEU)) + 
  geom_point()

# extract pairwise data for DEU and FRA, reshape/convert to numeric matrix
dt_DF <- rbindlist(list(dt["DEU"], dt["FRA"]))
c_dt_DF <- dcast(dt_DF, Date ~ iso3c, value.var = "eqntl")
rm(dt_DF)
c_dt_DF[, Date := NULL]
c_dt_DF <- as.matrix(c_dt_DF)

test <- fitCopula(copula = tCopula(), data = c_dt_DF)


