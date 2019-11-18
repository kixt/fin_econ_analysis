#####################################
# Financial Econometrics Term Paper #
#####################################

# Thore Petersen
# November 2019

## HPC Copula goodness of fit tests

.libPaths(new = "./R_libs/libs")
library(data.table)
library(copula)
library(parallel)


# Load data ---------------------------------------------------------------

load("./data/04_tmp.RData")
load("./data/cop_mix_fit.RData")
load("./data/mix_cop_vola.RData")


# Create cluster ----------------------------------------------------------

cl <- makeCluster(10)
clusterEvalQ(cl, .libPaths(new = "./R_libs/libs"))
clusterEvalQ(cl, library("copula")) # load copula package at nodes


# Prepare data subsets ----------------------------------------------------

# create matrix of pseudo observations to use in copula functions
dens <- dcast(dt, Date ~ iso3c, value.var = "qntl")
raw_dat <- dcast(dt, Date ~ iso3c, value.var = "res")

dens[, Date := NULL]
dens <- as.matrix(dens)
raw_dat[, Date := NULL]
raw_dat <- as.matrix(raw_dat)

dens_ls <- vector("list", (n * (n - 1) / 2))
raw_dat_ls <- vector("list", (n * (n - 1) / 2))
ls_i <- 1

for(i in 1:(n - 1)) {
  for(j in (i + 1):n) {
    dens_ls[[ls_i]] <- dens[, c(i, j)]
    raw_dat_ls[[ls_i]] <- raw_dat[, c(i, j)]
    ls_i <- ls_i + 1
  }
}; rm(i, j, ls_i)

# very simple regime indicator
dt[, high_vola := ifelse(sig > quantile(sig, 0.9), TRUE, FALSE), by = iso3c]

qntls_hv <- dcast(dt[high_vola == TRUE], Date ~ iso3c, value.var = "qntl")
qntls_lv <- dcast(dt[high_vola == FALSE], Date ~ iso3c, value.var = "qntl")

raw_dat_hv <- dcast(dt[high_vola == TRUE], Date ~ iso3c, value.var = "res")
raw_dat_lv <- dcast(dt[high_vola == FALSE], Date ~ iso3c, value.var = "res")

# make list of matrices of pairwise residuals by regima and pair
dens_vola_ls <- vector("list", 2)
names(dens_vola_ls) <- c("low", "high")
raw_dat_vola_ls <- vector("list", 2)
names(raw_dat_vola_ls) <- c("low", "high")

for(i in seq_along(dens_vola_ls)) {
  dens_vola_ls[[i]] <- vector("list", (n * (n - 1) / 2))
  raw_dat_vola_ls[[i]] <- vector("list", (n * (n - 1) / 2))
  
  ls_i <- 1
  for(j in 1:(n - 1)) {
    for(k in (j + 1):n) {
      tmp1 <- as.matrix(switch(i, "1" = qntls_lv, "2" = qntls_hv)[, -1])
      tmp2 <- as.matrix(switch(i, "1" = raw_dat_lv, "2" = raw_dat_hv)[, -1])
      dens_vola_ls[[i]][[ls_i]] <- na.exclude(tmp1[, c(j, k)])
      raw_dat_vola_ls[[i]][[ls_i]] <- na.exclude(tmp2[, c(j, k)])
      ls_i <- ls_i + 1
    }
  }
}; rm(i, j, k, ls_i, tmp1, tmp2)


# GoF for general model ---------------------------------------------------

gof_wrap <- function(copula, x) {
  tryCatch(
    gofCopula(copula = copula, x = x, N = 1000, method = "SnB", 
              test.method = "single"),
    error = function(e) e
  )
}
clusterExport(cl, "gof_wrap")


cop_res <- setClass(
  "cop_res", slots = list(
    pairname = "character", regime = "character", x = "matrix", u = "matrix", 
    mix_t_cop = "fitCopula", mix_cop = "fitCopula", gof_mix_t = "list", 
    gof_mix = "list"
  )
)

general_fit <- vector("list", (n * (n - 1) / 2))
for(i in seq_along(general_fit)) {
  general_fit[[i]] <- cop_res(
    pairname = paste0(dimnames(dens_ls[[i]])[[2]], collapse = "-"),
    regime = "general",
    u = dens_ls[[i]],
    x = raw_dat_ls[[i]],
    mix_t_cop = mc_fit[[i]]$mix_cop$t,
    mix_cop = mc_fit[[i]]$mix_cop$no_t
  )
}


gof_mix_t <- parLapply(cl, general_fit, function(y) gof_wrap(y@mix_t_cop@copula, x = y@x))
gof_mix <- parLapply(cl, general_fit, function(y) gof_wrap(y@mix_cop@copula, x = y@x))

for(i in seq_along(general_fit)) {
  general_fit[[i]]@gof_mix_t <- list(gof_mix_t[[i]])
  general_fit[[i]]@gof_mix <- list(gof_mix[[i]])
}

stopCluster(cl)

save(general_fit, file = "./data/general_gof.RData")
