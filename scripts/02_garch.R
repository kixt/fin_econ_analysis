#####################################
# Financial Econometrics Term Paper #
#####################################

# Thore Petersen
# October 2019

# fit AR(1)-GJR-GARCH(1, 1) to daily nominal stock returns
# esimate model with skew-t and skew-generalised error distribution residuals
# select model by country according to AIC, BIC

library(data.table)
library(rugarch)
library(ggplot2)


# Load data ---------------------------------------------------------------

load("./data/clean/clean.RData")
setkey(dt, iso3c, Date)

markets <- unique(dt$iso3c)
n <- length(markets)


# Test various specifications ---------------------------------------------

dfams <- c("norm", "t", "st", "sged")
specs <- vector("list", n)
names(specs) <- dfams

fit <- vector("list", length(dfams))
names(fit) <- dfams

# specify GJR-GARCH with different error distributions
for(d in dfams) {
  specs[[d]] <- ugarchspec(
    variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)),
    mean.model = list(armaOrder = c(1, 0), include.mean = TRUE),
    distribution.model = 
      switch(d, norm = "norm", t = "std", st = "sstd", sged = "sged")
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
# skew generalised error distr preferred in most countries, better accounts for 
# kurtosis/peakedness in the mode

pref_d <- ic[ic == "aic", .(iso3c, pref)]


# Export GARCH results ----------------------------------------------------

# extract standardized residuals, i.e. epsi_t / sig_t|sig_{t-1}, and sigma_t,
# extract it from the fit using the preferred distribution by country
dt[, res := as.numeric(residuals(fit[[pref_d[.GRP, pref]]][[.GRP]], standardize = TRUE)), 
   by = iso3c]
dt[, sig := as.numeric(sigma(fit[[pref_d[.GRP, pref]]][[.GRP]])), 
   by = iso3c]

save(dt, markets, n, pref_d, file = "./data/tmp/02_tmp.RData")
