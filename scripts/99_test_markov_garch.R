#####################################
# Financial Econometrics Term Paper #
#####################################

# Thore Petersen
# October 2019

## Test univariate Markov switching GARCH


library(data.table)
library(tseries)


# Import data -------------------------------------------------------------

test_imp <- read.csv("./data/raw/lse.csv")

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
bench <- garch(dat$r, c(0, 3))

plot(bench)

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
  dat[, diff2 := (x)^2]
  
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

own_ll_bounded <- optim(c(1, 0.1), qmle_ll, x = dat$r, lower = c(0, 0), method = "L-BFGS-B")


# Call Julia code ---------------------------------------------------------

library(XRJulia)
findJulia(test = TRUE)

library(readtext)

# read .jl source containing only the function
arch.jl <- readtext("./scripts/julia_src/arch.jl")[[2]]

# evaluate the code and load the Optim package
arch <- juliaEval(arch.jl)
jl_Optim <- juliaEval("using Optim")


# tell R that the evaluated function in the Julia workspace is a function usable 
# in R (argument order maintained from Julia)
arch_jl <- JuliaFunction(arch)

# evaluate function with inputs (keep track of types)
test <- arch_jl(dat$r, 3L, c(1, 0.5, 0.1, 0.01))

# get result
test2 <- juliaGet(test)








