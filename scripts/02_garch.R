#####################################
# Financial Econometrics Term Paper #
#####################################

# Thore Petersen
# October 2019

library(data.table)
library(rugarch)
library(ggplot2)


# Load data ---------------------------------------------------------------

load("./data/clean/clean.RData")
setkey(dt, iso3c, Date)

markets <- unique(dt$iso3c)
n <- length(markets)


# Test various specifications ---------------------------------------------

specs <- vector("list", n)
names(specs) <- markets

fit <- vector("list", n)
names(fit) <- markets

## DEU
specs[["DEU"]] <- vector("list", 10)
specs[["DEU"]][[3]] <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(1, 1), include.mean = TRUE),
  distribution.model = "sstd"
)

fit[["DEU"]] <- vector("list", 10)

fit[["DEU"]][[3]] <- ugarchfit(specs[["DEU"]][[3]], dt["DEU", ret], solver = "hybrid")

plot(fit[["DEU"]][[3]])
