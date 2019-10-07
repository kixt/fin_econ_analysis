#####################################
# Financial Econometrics Term Paper #
#####################################

# Thore Petersen
# October 2019

## Test univariate Markov switching GARCH


library(data.table)
library(tseries)


# Import data -------------------------------------------------------------

test_imp <- read.csv("./data/lse.csv")

dat <- data.table(
  date = as.Date(as.character(test_imp$Date)),
  x = test_imp$Price
)

# some prices too small by factor 100
dat[x / 100 < 1, x := x * 100]

# calculate returns
dat[, r := 100 * c(NA_real_, diff(log(x)))]
dat <- zoo::na.trim(dat)

rm(test_imp)


# Basic plots and descriptives --------------------------------------------

# library(ggplot2)
# ggplot(dat, aes(x = date, y = r)) +
#   geom_line()
# 
# ggplot(dat, aes(x = r)) +
#   geom_histogram()
# 
# acf(dat$r)
# acf(dat$r^2)


# ARCH(1) -----------------------------------------------------------------

# benchmark estimation of ARCH(1)
bench <- garch(dat$r, c(0, 1))


qmle_ll <- function(theta, x) {
  # function to calculate negative log-likelihood of ARCH(1) model for x
  #
  # Args:
  #   theta: vector of parameters, first element is a0, second element is a1
  #   x: series for which the ARCH(1) is to be estimated
  # Returns:
  #   the negative log-likelihood for minimisation
  
  n <- length(x)
  xbar <- mean(x)
  
  dat <- data.table(x = x)
  dat[, diff2 := (x - xbar)^2]
  
  a0 <- theta[1]
  a1 <- theta[2]
  
  uv <- mean(dat[, diff2])
  cv <- a0 + a1 * uv
  
  ll <- -1 / 2 * (log(cv) + dat[1, diff2] / cv)
  
  for(i in 2:n) {
    cv <- a0 + a1 * dat[i - 1, diff2]
    ll <- ll - 1 / 2 * (log(cv) + dat[i, diff2] / cv)
  }
  return(-ll)
}

own_ll <- optim(c(1.5, 0.4), qmle_ll, x = dat$r)


# Call Julia code ---------------------------------------------------------

library(XRJulia)
findJulia(test = TRUE)

library(readtext)

# .jl file is OK, problem with R connection
arch1_nll.jl <- readtext("./scripts/julia_src/arch1_nll.jl")[[2]]

arch1_nll <- juliaEval(arch1_nll.jl)
