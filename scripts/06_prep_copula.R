#####################################
# Financial Econometrics Term Paper #
#####################################

# Thore Petersen
# October, November 2019

library(data.table)
library(copula)


# Load data ---------------------------------------------------------------

load("./data/tmp/04_tmp.RData")


# Empirical quantiles -----------------------------------------------------

# create matrix of pseudo observations to use in copula functions
dis <- dcast(dt, Date ~ iso3c, value.var = "Fh")

dis[, Date := NULL]
dis <- as.matrix(dis)


# Estimate mixture copulae ------------------------------------------------

dis_ls <- vector("list", (n * (n - 1) / 2))
ls_i <- 1

for(i in 1:(n - 1)) {
  for(j in (i + 1):n) {
    dis_ls[[ls_i]] <- dis[, c(i, j)]
    names(dis_ls)[ls_i] <- paste(colnames(dis)[c(i, j)], collapse = "-")
    ls_i <- ls_i + 1
  }
}

# laod data from HPC calculation
load("./data/tmp/cop_mix_fit.RData")
load("./data/tmp/cop_mix_fit_vola.RData")

mc_vola_fit1 <- mc_vola_fit
mc_fit1 <- mc_fit

load("./data/tmp/cop_mix_fit_as.RData")
load("./data/tmp/cop_mix_fit_vola_as.RData")

# grad all parameter estimates, histogram of differences for diff. starting values
hist(
  unlist(
    lapply(
      mc_vola_fit1,
      function(x) lapply(x,
                         function(y) lapply(y[[2]],
                                            function(z) z@estimate))))
  -
    unlist(
      lapply(
        mc_vola_fit,
        function(x) lapply(x,
                           function(y) lapply(y[[2]],
                                              function(z) z@estimate)))),
  30
  )
# mostly stable to changes in starting values

# calculate AIC for all models, find minimum; I should rethink these nested lists
mc_vola_aic <- 
  lapply(
    mc_vola_fit, 
    function(x) lapply(x, 
                       function(y) lapply(y, 
                                          function(z) lapply(z, AIC))
                       )
    )
mc_vola_min_aic <- 
  lapply(
    mc_vola_aic, 
    function(x) lapply(x, function(y) which.min(unlist(y)))
    )

# extract best fitting copulas by AIC
best_cop <- vector("list", 2)
names(best_cop) <- c("low", "high")
for(i in seq_along(best_cop)) {
  best_cop[[i]] <- vector("list", 10)
  for(j in seq_along(best_cop[[i]])) {
    best_cop[[i]][[j]] <- unlist(mc_vola_fit[[i]][[j]])[mc_vola_min_aic[[i]][[j]]]
  }
}

save(best_cop, file = "./data/tmp/best_cop.RData")


# Asymmetric volatilities -------------------------------------------------

load("./data/tmp/cop_mix_fit_avola.RData")

mc_avola_aic <- 
  lapply(
    mc_avola_fit, function(x) lapply(x, function(y) lapply(y, AIC))
  )
mc_avola_min_aic <- 
  lapply(
    mc_avola_aic, 
    function(y) which.min(unlist(y))
  )

best_cop_avola <- vector("list", n * (n - 1))
for(i in seq_along(best_cop_avola)) {
  best_cop_avola[[i]] <- unlist(mc_avola_fit[[i]])[mc_avola_min_aic[[i]]]
}

save(best_cop_avola, file = "./data/tmp/best_cop_avola.RData")
