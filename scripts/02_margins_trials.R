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

load("./data/clean.RData")
setkey(dt, iso3c, Date)


# Set up GARCH specification ----------------------------------------------

markets <- unique(dt$iso3c)
n <- length(markets)

# set up lists for specifications and results of GARCH estimation
garch_spec <- vector("list", n)
garch_fit <- vector("list", n)
names(garch_spec) <- names(garch_fit) <- markets

test_spec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(1, 1), include.mean = TRUE),
  distribution.model = "std" 
)
# use e.g. fixed.pars = list(shape = 5) to set dof of t distribution.model, 
# otherwise estimted via ML

# estimate GARCH models
for(m in markets) {
  garch_spec[[m]] <- test_spec
  
  # select rows in data.table by key, fast binary search
  garch_fit[[m]] <- ugarchfit(garch_spec[[m]], dt[m, ret], solver = "hybrid")
}


# Extract GARCH residuals -------------------------------------------------

# extract standardized residuals, i.e. epsi_t / sig_t|sig_{t-1}, and sigma_t
dt[, res := unlist(lapply(garch_fit, residuals, standardize = TRUE))]
dt[, sig := unlist(lapply(garch_fit, sigma))]

# # KDE for residuals
# res_kde <- vector("list", n)
# names(res_kde) <- markets
# 
# for(m in markets) {
#   res_kde[[m]] <- density(dt[m, res])
# }

ggplot(dt) +
  geom_density(aes(x = res)) +
  facet_wrap(iso3c ~ .)

# ECDF of residuals
# copula::pobs() returns ecdf values scaled by n/(n+1), st. border cases fall into unit cube
dt[, eqtl := pobs(res), by = iso3c]

# create data.table of empirical quantiles to check dependence
eqtls <- dcast(dt, Date ~ iso3c, value.var = "eqtl")

# characteristic pattern for DEU--FRA, almost linear, with clustering in corners,
# i.e. strong tail dependence
ggplot(eqtls, aes(x = FRA, y = GRC)) + 
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

gen_dens <- function(fit, grp_name, grp_id, n = 100) {
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

