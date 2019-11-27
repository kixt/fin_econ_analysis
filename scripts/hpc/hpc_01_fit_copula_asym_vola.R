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
                           no_t_mix = c("yes", "no", "both"), mc_t_w = NULL) {
  # function to fit a mixture copula on u, fits all copulas individually and
  # uses parameter estimates as starting values for estimation of the mixture,
  # mixture estimated by default with SANN algorithm, more robust to fucky
  # gradients than BFGS-like algorithms
  # Args:
  #   u: nxd matrix of pseudo-observations (in the spirit of the copula package)
  #   copulas: character vector of copulas to be included in the mixture
  #   optim.method: passed to optim()
  #   no_t_mix: should be excluded from the mixture copula estimation?
  #   mc_t_w: numeric vector, starting weights for mixture copula including t, 
  #   currently only implemented for no_t_mix = "both"
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
    mc_t <- mixCopula(lapply(cops, function(x) x@copula), w = mc_t_w)
    mc_no_t <- mixCopula(lapply(cops[which(names(cops) != "t")], function(x) x@copula))
    
    mc_fit_t <- tryCatch(
      fitCopula(mc_t, u, method = "mpl", optim.method = optim.method), 
      error = function(e) e
      )
    mc_fit_no_t <- tryCatch(
      fitCopula(mc_no_t, u, method = "mpl", optim.method = optim.method),
      error = function(e) e
      )
    
    mc_fit <- list(t = mc_fit_t, no_t = mc_fit_no_t)
  }
  res <- list(single_cops = cops, mix_cop = mc_fit)
  
  return(res)
}

fmc_error_wrap <- function(u, copulas = c("frank", "t", "gumbel", "clayton"),
                           optim.method = "L-BFGS-B", no_t_mix = "both", ...) {
  # wrapper around fit_mix_copula() to include error handling capabilities, 
  # primarily for deployment in (par)lapply
  tryCatch(
    fit_mix_copula(u = u, copulas = copulas, optim.method = optim.method, 
                   no_t_mix = no_t_mix, ...),
    error = function(e) e
    )
}


# Load data ---------------------------------------------------------------

load("./data/tmp/04_tmp.RData")


# Create cluster ----------------------------------------------------------

cl <- makeCluster(3)
clusterEvalQ(cl, .libPaths(new = "./R_libs/libs"))
clusterEvalQ(cl, library("copula")) # load copula package at nodes
clusterExport(cl, "fit_mix_copula")


# Estimate mixture copulae by volatility regime ---------------------------

# very simple regime indicator
dt[, high_vola := ifelse(sig > quantile(sig, 0.9), TRUE, FALSE), by = iso3c]

dis_ls <- vector("list", n * (n - 1))

ls_i <- 1
for(i in 1:n) {
  for(j in 1:n) {
    if(i != j) {
      ic <- markets[i]
      jc <- markets[j]
      tmp1 <- dt[ic][high_vola == TRUE, .(Date, Fh)]
      tmp2 <- dt[jc][high_vola == FALSE, .(Date, Fh)]
      tmpm <- merge(tmp1, tmp2, by = "Date", all = FALSE)
      colnames(tmpm)[2:3] <- c(ic, jc)
      
      dis_ls[[ls_i]] <- as.matrix(tmpm[, c(2, 3)])
      ls_i <- ls_i + 1
    }
  }
}

mc_avola_fit <- parLapply(cl, dis_ls, fmc_error_wrap)

# errors in 5 and 19, mix cops with t -> manually set starting weights

mc_avola_fit[[5]] <- fmc_error_wrap(
  dis_ls[[5]], optim.method = "L-BFGS-B", mc_t_w = c(0.1, 0.5, 0.2, 0.2)
  )

mc_avola_fit[[19]] <- fmc_error_wrap(
  dis_ls[[19]], optim.method = "L-BFGS-B", mc_t_w = c(0.2, 0.3, 0.2, 0.3)
)


stopCluster(cl)

save(mc_avola_fit, file = "./data/results/cop_mix_fit_avola.RData")
