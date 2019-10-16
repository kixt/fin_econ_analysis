#####################################
# Financial Econometrics Term Paper #
#####################################

# Thore Petersen
# October 2019

library(data.table)
library(rugarch)
library(ggplot2)


# Load data ---------------------------------------------------------------

load("./data/clean.RData")
setkey(dt, iso3c, Date)

# Set up GARCH specification ----------------------------------------------

markets <- unique(dt$iso3c)

# set up recursive list for specifications and results of GARCH estimation
garch <- vector("list", 2)
names(garch) <- c("specs", "results")
garch[[1]] <- garch[[2]] <- vector("list", length(markets))
names(garch[[1]]) <- names(garch[[2]]) <- markets

test_spec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(1, 1), include.mean = TRUE),
  distribution.model = "std" 
)
# use e.g. fixed.pars = list(shape = 5) to set dof of t distribution.model, 
# otherwise estimted via ML

for(m in markets) {
  x <- dt[m, ret]
  garch[[2]][[m]] <- ugarchfit(test_spec, x, solver = "hybrid")
}; rm(m, x)



# Extract GARCH residuals -------------------------------------------------

dt[, res := unlist(lapply(garch[[2]], residuals, standardize = TRUE))]

