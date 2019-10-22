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

dfams <- c("norm", "t", "st")
specs <- vector("list", n)
names(specs) <- dfams

fit <- vector("list", length(dfams))
names(fit) <- dfams

# specify GJR-GARCH with different error distributions
for(d in dfams) {
  specs[[d]] <- ugarchspec(
    variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)),
    mean.model = list(armaOrder = c(1, 0), include.mean = TRUE),
    distribution.model = switch(d, norm = "norm", t = "std", st = "sstd")
    )
  
  fit[[d]] <- vector("list", n)
  names(fit[[d]]) <- markets
  
  for(m in markets) {
    fit[[d]][[m]] <- ugarchfit(specs[[d]], dt[m, ret], solver = "hybrid")
  }
}

# save(fit, file = "./data/tmp/gjr_garch_fit.RData")
# load("./data/tmp/gjr_garch_fit.RData")

# extract information criteria
ic <- expand.grid(markets, c("aic", "bic"), KEEP.OUT.ATTRS = FALSE)
setDT(ic)
colnames(ic) <- c("iso3c", "ic")
setkey(ic, iso3c, ic)

ic[, eval(dfams) := NA_real_]
for(m in markets) {
  for(d in dfams) {
    ic_vals <- infocriteria(fit[[d]][[m]])
    ic[m, eval(d) := ic_vals[1:2]]
  }
}; rm(ic_vals)

ic[, pref := colnames(.SD)[which.min(.SD)], by = c("iso3c", "ic")]
# skewed-t preferred in all countries according to both ICs


# Export GARCH results ----------------------------------------------------

# extract standardized residuals, i.e. epsi_t / sig_t|sig_{t-1}, and sigma_t
dt[, res := unlist(lapply(fit[["st"]], residuals, standardize = TRUE))]
dt[, sig := unlist(lapply(fit[["st"]], sigma))]

save(dt, markets, n, file = "./data/tmp/02_tmp.RData")
