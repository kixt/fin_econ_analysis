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

load("./data/clean/clean.RData")
setkey(dt, iso3c, Date)


# Characterise GARCH residuals --------------------------------------------

# ECDF of residuals
# copula::pobs() returns ecdf values scaled by n/(n+1), st. everything in unit cube
dt[, eqtl := pobs(res), by = iso3c]

# create data.table of empirical quantiles to check dependence
eqtls <- dcast(dt, Date ~ iso3c, value.var = "eqtl")

# characteristic pattern for DEU--FRA, almost linear, with clustering in corners,
# i.e. strong tail dependence
ggplot(eqtls, aes(x = FRA, y = DEU)) + 
  geom_point()


# Fitting distributions to residuals -- skew-t ----------------------------

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

# fit skewed-t-distribution to all residuals
fit_st <- vector("list", n)
names(fit_st) <- markets
for(m in markets) {
  fit_st[[m]] <- fitdist(dt[m, res], "st",
                         start = list(xi = 0, omega = 1, alpha = 0, nu = 5))
}


# Fitting distributions to residuals -- t ---------------------------------

# fit t-distribution to all residuals
fit_t <- vector("list", n)
names(fit_t) <- markets
for(m in markets) {
  fit_t[[m]] <- fitdist(dt[m, res], "t", start = list(df = 5))
}


# Compare fits of distributions -------------------------------------------

# create data.table for the values of the information criteria of the fits
dfams <- c("st", "t")
ic_dens <- expand.grid(markets, c("aic", "bic"), KEEP.OUT.ATTRS = FALSE)
setDT(ic_dens)
colnames(ic_dens) <- c("iso3c", "ic_dens")
setkey(ic_dens, iso3c, ic_dens)

# extract values for AIC, BIC from the fit results
ic_dens[, eval(dfams) := NA_real_]
for(m in markets) {
  for(d in dfams) {
    fit_list <- as.name(paste0("fit_", d))
    # not too firm with do.call() yet, but this works
    ic_vals <- do.call(function(x) c(aic = x[[m]]$aic, bic = x[[m]]$bic), 
                       list(x = fit_list))
    ic_dens[m, eval(d) := ic_vals]
  }
}

# determine which family gives the minimum ic_dens
ic_dens[, pref := colnames(.SD)[which.min(.SD)], by = c("iso3c", "ic_dens")]


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
  densities[[m]] <- gen_dens(fit_st[[m]], "iso3c", m)
}
densities <- rbindlist(densities)

ggplot() +
  geom_density(aes(x = res), dt) +
  geom_line(aes(x = x, y = density), densities, colour = "red") +
  facet_wrap(iso3c ~ .)





