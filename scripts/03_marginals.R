#####################################
# Financial Econometrics Term Paper #
#####################################

# Thore Petersen
# October 2019

library(data.table)
library(sn) # skew distributions
library(sgt) # skew generalised t and related
library(fitdistrplus) # extends MASS::fitdistr(); probably overkill, we'll see
library(ggplot2)


# Load data ---------------------------------------------------------------

load("./data/tmp/02_tmp.RData")
setkey(dt, iso3c, Date)


# Fitting distributions -- Preliminaries ----------------------------------

dfams <- c("st", "sged")

## wrappers around functions for skew-t and skew-gen.-t for fitting

pst <- function(q, xi=0, omega=1, alpha=0, nu=Inf, dp=NULL, method=0) {
  # wrapper for sn::pst() to use in fitdistrplus::fitdist(), which has specific 
  # requirements for the pst() function:
  #   - first argument must be named q
  #   - return NULL when q = NULL (sn::pst() throws an error)
  #   - names must match exactly (sn::pst() allows ... for functionality I don't 
  #   need, thus does not throw "unused arguments" error)
  
  if(length(q) == 0) return(NULL)
  sn::pst(x = q, xi, omega, alpha, nu, dp, method)
}

psged <- function(q, mu = 0, sigma = 1, lambda = 0, p = 2) {
  # wrapper around sgt::psgt() for fitdistrplus::fitdist, specific requirements 
  # for arguments, and fix q = Inf to use SGED, not SGT
  
  if(length(q) == 0) return(NULL)
  sgt::psgt(quant = q, mu = mu, sigma = sigma, lambda = lambda, p = p, q = Inf)
}

dsged <- function(x, mu = 0, sigma = 1, lambda = 0, p = 2) {
  # wrapper around sgt::dsgt()
  
  if(length(x) == 0) return(NULL)
  sgt::dsgt(x = x, mu = mu, sigma = sigma, lambda = lambda, p = p, q = Inf)
}

## starting values for the various distributions
svals <- vector("list", length(dfams))
names(svals) <- dfams

svals[["st"]] <- list(xi = 0, omega = 1, alpha = 0, nu = 5)
svals[["sged"]] <- list(mu = 0, sigma = 1, lambda = 0, p = 1)


# Fitting distributions ---------------------------------------------------

fit <- vector("list", length(dfams))
names(fit) <- dfams

for(d in dfams) {
  fit[[d]] <- vector("list", n)
  names(fit[[d]]) <- markets
  
  for(m in markets) {
    fit[[d]][[m]] <- fitdist(dt[m, res], d, start = svals[[d]])
  }
}


# Compare fits of distributions -------------------------------------------

# create data.table for the values of the information criteria of the fits
ic <- expand.grid(markets, c("aic", "bic"), KEEP.OUT.ATTRS = FALSE)
setDT(ic)
colnames(ic) <- c("iso3c", "ic")
setkey(ic, iso3c, ic)

# extract values for AIC, BIC from the fit results
ic[, eval(dfams) := NA_real_]
for(m in markets) {
  for(d in dfams) {
    
    # not too firm with do.call() yet, but this works
    ic_vals <- do.call(function(x) c(aic = x[[m]]$aic, bic = x[[m]]$bic), 
                       list(x = fit[[d]]))
    ic[m, eval(d) := ic_vals]
  }
}; rm(d, m, ic_vals)

# determine which family gives the minimum ic
ic[, pref := colnames(.SD)[which.min(.SD)], by = c("iso3c", "ic")]
# both information criteria agree for all countries

# agrees with preferred error distribution in rugarch (sensibly)
pref_d[, pref_fit := ic[ic == "aic", pref]]


# Significance tests of parameter estimates -------------------------------

dparams <- vector("list", length(dfams))
names(dparams) <- dfams

for(d in dfams) {
  dparams[[d]] <- expand.grid(markets, names(fit[[d]][[1]]$estimate), 
                              KEEP.OUT.ATTRS = FALSE)
  
  colnames(dparams[[d]]) <- c("iso3c", "param")
  setDT(dparams[[d]])
  setkey(dparams[[d]], iso3c, param)
  
  dparams[[d]][, c("est", "sd") := .(fit[[d]][[.GRP]]$estimate, fit[[d]][[.GRP]]$sd),
               by = iso3c]
  
  dparams[[d]][, tstat := est / sd]
  dparams[[d]][, pval := 1 - pnorm(abs(tstat))]
}


# Export fit results ------------------------------------------------------

fit <- fit_st

# calculate cdf values at the empirical observations, complicated looking call
# to have matching argument names
# split(unname(), names()) returns named list from named vector
dt[, 
   qntl := do.call(sn::pst, 
                   args = c(list(x = res), split(unname(fit[[.GRP]]$estimate), 
                                                 names(fit[[.GRP]]$estimate)))), 
   by = iso3c]

save(dt, n, markets, file = "./data/tmp/03_tmp.RData")


# Plot KDE and fitted parametric density ----------------------------------

gen_dens <- function(fit, grp_name, grp_id, n = 200) {
  # function to generate points to plot density estimate in fit, designed to 
  # generate data to overlay in ggplot() + facet_wrap(grp_name ~ .) instead of
  # stat_function()
  # 
  # Args:
  #   fit: object returned from fitdistrplus::fitdist()
  #   grp_name: column name of the group identifiers in original data
  #   grp_id: group identifier for which fit was estimated
  # Returns:
  #   a data.table with a grid of x, the corresponding density, and the group
  
  if(class(fit) != "fitdist") stop("fit must be class fitdist")
  
  dfn <- paste0("d", fit$distname)
  dargs <- split(unname(fit$estimate), names(fit$estimate))
  
  grid <- data.table(x = with(fit, seq(min(data), max(data), length = n)))
  grid[, density := do.call(dfn, c(list(x = x), dargs))]
  grid[, eval(grp_name) := grp_id]
  
  return(grid)
}

# create data.table for densities to overlay on plot
densities <- vector("list", n)
names(densities) <- markets
for(m in markets) {
  densities[[m]] <- gen_dens(fit_sged[[m]], "iso3c", m)
}
densities <- rbindlist(densities)

ggplot() +
  geom_density(aes(x = res), dt) +
  geom_line(aes(x = x, y = density), densities, colour = "red") +
  facet_wrap(iso3c ~ .)

# skew-t not very uniform, fit can't account for high probability around mode
ggplot(dt, aes(x = qntl, group = iso3c)) + 
  geom_histogram() +
  facet_wrap(iso3c ~.)
