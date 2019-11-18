#####################################
# Financial Econometrics Term Paper #
#####################################

# Thore Petersen
# November 2019

## HPC Copula estimation

.libPaths(new = "./R_libs/libs")
library(data.table)
library(copula)
library(parallel)


# Functions ---------------------------------------------------------------

fit_mix_copula <- function(u, copulas = c("frank", "t", "gumbel", "clayton"),
                           optim.method = "L-BFGS-B", 
                           no_t_mix = c("yes", "no", "both")) {
  # function to fit a mixture copula on u, fits all copulas individually and
  # uses parameter estimates as starting values for estimation of the mixture,
  # mixture estimated by default with SANN algorithm, more robust to fucky
  # gradients than BFGS-like algorithms
  # Args:
  #   u: nxd matrix of pseudo-observations (in the spirit of the copula package)
  #   copulas: character vector of copulas to be included in the mixture
  #   optim.method: passed to optim()
  #   no_t_mix: should be excluded from the mixture copula estimation?
  # Returns:
  #   a list of the individual fits of all copulas and the fit of the mixture
  
  nc <- length(copulas)
  cops <- vector("list", nc)
  names(cops) <- copulas
  for(i in 1:nc) {
    # switch evaluates the expression to assign a copula-object to c_obj
    c_obj <- switch(copulas[i],
                    "frank" = frankCopula(),
                    "t" = tCopula(),
                    "gumbel" = gumbelCopula(),
                    "clayton" = claytonCopula())
    cops[[i]] <- do.call(fitCopula, args = list(c_obj, u, optim.method = "BFGS"))
  }
  
  # get the individually fit copula objects with parameter estimates
  if(match.arg(no_t_mix) == "yes") {
    mc <- mixCopula(lapply(cops[which(names(cops) != "t")], function(x) x@copula))
    mc_fit <- fitCopula(mc, u, method = "mpl", optim.method = optim.method)
  } else if(match.arg(no_t_mix) == "no") {
    mc <- mixCopula(lapply(cops, function(x) x@copula))
    mc_fit <- fitCopula(mc, u, method = "mpl", optim.method = optim.method)
  } else if(match.arg(no_t_mix) == "both") {
    mc_t <- mixCopula(lapply(cops, function(x) x@copula))
    mc_no_t <- mixCopula(lapply(cops[which(names(cops) != "t")], function(x) x@copula))
    
    mc_fit_t <- fitCopula(mc_t, u, method = "mpl", optim.method = optim.method)
    mc_fit_no_t <- fitCopula(mc_no_t, u, method = "mpl", optim.method = optim.method)
    
    mc_fit <- list(t = mc_fit_t, no_t = mc_fit_no_t)
  }
  res <- list(single_cops = cops, mix_cop = mc_fit)
  
  return(res)
}

fmc_error_wrap <- function(u, copulas = c("frank", "t", "gumbel", "clayton"),
                           optim.method = "L-BFGS-B", no_t_mix = "both") {
  # wrapper around fit_mix_copula() to include error handling capabilities, 
  # primarily for deployment in (par)lapply
  tryCatch(
    fit_mix_copula(u = u, copulas = copulas, optim.method = optim.method, 
                   no_t_mix = no_t_mix),
    error = function(e) e
    )
}


# Load data ---------------------------------------------------------------

load("./data/tmp/04_tmp.RData")
# load("./data/cop_mix_fit.RData")

# create matrix of pseudo observations to use in copula functions
dis <- dcast(dt, Date ~ iso3c, value.var = "Fh")

dis[, Date := NULL]
dis <- as.matrix(dis)


# Create cluster ----------------------------------------------------------

cl <- makeCluster(3)
clusterEvalQ(cl, .libPaths(new = "./R_libs/libs"))
clusterEvalQ(cl, library("copula")) # load copula package at nodes
clusterExport(cl, "fit_mix_copula")


# Estimate over both regimes ----------------------------------------------

dis_ls <- vector("list", (n * (n - 1) / 2))
ls_i <- 1

for(i in 1:(n - 1)) {
  for(j in (i + 1):n) {
    dis_ls[[ls_i]] <- dis[, c(i, j)]
    ls_i <- ls_i + 1
  }
}

mc_fit <- parLapply(cl, dis_ls, fmc_error_wrap)

save(mc_fit, file = "./data/results/cop_mix_fit.RData")


# Estimate mixture copulae by volatility regime ---------------------------

# very simple regime indicator
dt[, high_vola := ifelse(sig > quantile(sig, 0.9), TRUE, FALSE), by = iso3c]

dis_hv <- dcast(dt[high_vola == TRUE], Date ~ iso3c, value.var = "Fh")
dis_lv <- dcast(dt[high_vola == FALSE], Date ~ iso3c, value.var = "Fh")

# make list of matrices of pairwise residuals by regima and pair
dis_vola_ls <- vector("list", 2)
names(dis_vola_ls) <- c("low", "high")

for(i in seq_along(dis_vola_ls)) {
  dis_vola_ls[[i]] <- vector("list", (n * (n - 1) / 2))
  
  ls_i <- 1
  for(j in 1:(n - 1)) {
    for(k in (j + 1):n) {
      tmp <- as.matrix(switch(i, "1" = dis_lv, "2" = dis_hv)[, -1])
      dis_vola_ls[[i]][[ls_i]] <- na.exclude(tmp[, c(j, k)])
      ls_i <- ls_i + 1
    }
  }
}; rm(i, j, k, ls_i, tmp)

mc_vola_fit <- vector("list", 2)
names(mc_vola_fit) <- c("low", "high")

for(i in seq_along(mc_vola_fit)) {
  mc_vola_fit[[i]] <- parLapply(cl, dis_vola_ls[[i]], fmc_error_wrap)
}

stopCluster(cl)

save(mc_vola_fit, file = "./data/results/cop_mix_fit_vola.RData")
