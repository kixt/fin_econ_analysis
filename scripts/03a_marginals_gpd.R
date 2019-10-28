#####################################
# Financial Econometrics Term Paper #
#####################################

# Thore Petersen
# October 2019

# fit generalised Pareto distribution to lower and upper tail of GARCH-filtered 
# residuals, use ECDF between the tails
# fit not very good, shape parameter estimate with high SE in all countries
# tail behaviour potentially with less systematic overprediction of extreme tail 
# events compared to ST/SGED

library(data.table)
library(POT)
library(ggplot2)


# Load data ---------------------------------------------------------------

load("./data/tmp/02_tmp.RData")
setkey(dt, iso3c, Date)


# Define classes ----------------------------------------------------------

# define class for fitting results of single tail; POT package written for maxima, 
# all results for lower tail are based on -1 * (lower tail), including threshold
tail_fit <- setClass(
  "tail_fit", 
  slots = list(distr = "character", tail = "character", thresh = "numeric", 
               params = "list")
  )

# define class to combine objects for each tail
tails <- setClass(
  "tails", 
  slots = list(lower = "tail_fit", upper = "tail_fit", quantile = "numeric")
  )


# Identify tails and fit GPD per tail -------------------------------------

q <- 0.05
tail_threshs <- tapply(dt$res, dt$iso3c, quantile, probs = c(q, 1 - q))

dt[, tail := "none"]
for(m in markets) {
  dt[iso3c == m & res <= tail_threshs[[m]][1], tail := "lower"]
  dt[iso3c == m & res >= tail_threshs[[m]][2], tail := "upper"]
}
setkey(dt, iso3c, tail) # set tail as key for faster subsetting later

fit_gpd <- vector("list", n)
names(fit_gpd) <- markets

fitted_tails <- vector("list", n)
names(fitted_tails) <- markets

for(m in markets) {
  fit_gpd[[m]] <- vector("list", 2)
  names(fit_gpd[[m]]) <- c("lower", "upper")
  
  # fit GPD per country and tail
  for(t in c("lower", "upper")) {
    x <- dt[.(m, t), res]
    if(t == "lower") x <- -1 * x
    fit_gpd[[m]][[t]] <- fitgpd(x, min(x))
  }

  # extract relevant fitting results into objects of classes defined above
  fitted_tails[[m]] <- tails(
    quantile = q,
    lower = tail_fit(
      distr = "gpd",
      tail = "lower",
      thresh = fit_gpd[[m]][["lower"]]$threshold,
      params = as.list(fit_gpd[[m]][["lower"]]$fitted.values)
    ),
    upper = tail_fit(
      distr = "gpd",
      tail = "upper",
      thresh = fit_gpd[[m]][["upper"]]$threshold,
      params = as.list(fit_gpd[[m]][["upper"]]$fitted.values)
    )
  )
}


# Evaluate CDF of combination of truncated distributions ------------------

pgpd_ecdf_mix <- function(dat, tails, x = "res", t = "tail") {
  # function to calculate CDF combined from parametric GPD estimatations for tails,
  # and ecdf for middle observations
  # Args:
  #   dt: data.table with a column containing the random variable of interest and
  #     a factor column to identify the tail
  #   tails: an object of class tails containing the results of fitting a GPD
  #   x: column of the random variable
  #   t: column of the factor identifying the tail, must be {lower, none, upper}

  q <- tails@quantile
  
  # for deployment in dt[, j = , by = ], makes implicit use of `:=` in .SD
  dt <- copy(dat)
  
  xl <- dt[get(t) == "lower", get(x)]
  xm <- dt[get(t) == "none", get(x)]
  xu <- dt[get(t) == "upper", get(x)]
  
  # calculate cdf for lower tail, use parameter estimates from tails
  pl <- do.call(
    "pgpd", 
    args = c(list(q = (-1) * xl, loc = tails@lower@thresh), tails@lower@params)
    )
  pl <- pl * q # map (0, 1) -> (0, q)
  dt[get(t) == "lower", p_est := pl]
  
  # calculate ecdf in the middle, no real need for parametric fit here
  Fhat <- ecdf(xm)
  pm <- Fhat(xm) * (1 - 2 * q) + q # map (0, 1) -> (q, 1-q)
  dt[get(t) == "none", p_est := pm]
  
  pu <- do.call(
    "pgpd", 
    args = c(list(q = xu, loc = tails@upper@thresh), tails@upper@params)
    )
  pu <- 1 + q * (pu - 1) # map (0, 1) -> (1-q, 1)
  dt[get(t) == "upper", p_est := pu]
  
  return(dt[, p_est])
}


dt[, p_est := pgpd_ecdf_mix(.SD, fitted_tails[[.GRP]]), by = iso3c]

# tail fit not very good either, also high SE in parameter estimates of shape
ggplot(dt, aes(x = p_est)) + 
  geom_histogram(bins = 100) +
  facet_wrap(iso3c ~ .)


