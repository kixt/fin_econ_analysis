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

load("./data/tmp/04_tmp.RData")
load("./data/tmp/best_cop.RData")


# Create cluster ----------------------------------------------------------

cl <- makeCluster(10)
clusterEvalQ(cl, .libPaths(new = "./R_libs/libs"))
clusterEvalQ(cl, library("copula")) # load copula package at nodes


# Prepare data subsets ----------------------------------------------------

# create matrix of pseudo observations to use in copula functions
dis <- dcast(dt, Date ~ iso3c, value.var = "Fh")
raw_dat <- dcast(dt, Date ~ iso3c, value.var = "res")

dis[, Date := NULL]
dis <- as.matrix(dis)
raw_dat[, Date := NULL]
raw_dat <- as.matrix(raw_dat)

dis_ls <- vector("list", (n * (n - 1) / 2))
raw_dat_ls <- vector("list", (n * (n - 1) / 2))
ls_i <- 1

for(i in 1:(n - 1)) {
  for(j in (i + 1):n) {
    dis_ls[[ls_i]] <- dis[, c(i, j)]
    raw_dat_ls[[ls_i]] <- raw_dat[, c(i, j)]
    ls_i <- ls_i + 1
  }
}; rm(i, j, ls_i)

# very simple regime indicator
dt[, high_vola := ifelse(sig > quantile(sig, 0.9), TRUE, FALSE), by = iso3c]

Fhs_hv <- dcast(dt[high_vola == TRUE], Date ~ iso3c, value.var = "Fh")
Fhs_lv <- dcast(dt[high_vola == FALSE], Date ~ iso3c, value.var = "Fh")

raw_dat_hv <- dcast(dt[high_vola == TRUE], Date ~ iso3c, value.var = "res")
raw_dat_lv <- dcast(dt[high_vola == FALSE], Date ~ iso3c, value.var = "res")

# make list of matrices of pairwise residuals by regima and pair
dis_vola_ls <- vector("list", 2)
names(dis_vola_ls) <- c("low", "high")
raw_dat_vola_ls <- vector("list", 2)
names(raw_dat_vola_ls) <- c("low", "high")

for(i in seq_along(dis_vola_ls)) {
  dis_vola_ls[[i]] <- vector("list", (n * (n - 1) / 2))
  raw_dat_vola_ls[[i]] <- vector("list", (n * (n - 1) / 2))
  
  ls_i <- 1
  for(j in 1:(n - 1)) {
    for(k in (j + 1):n) {
      tmp1 <- as.matrix(switch(i, "1" = Fhs_lv, "2" = Fhs_hv)[, -1])
      tmp2 <- as.matrix(switch(i, "1" = raw_dat_lv, "2" = raw_dat_hv)[, -1])
      dis_vola_ls[[i]][[ls_i]] <- na.exclude(tmp1[, c(j, k)])
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
    copula = "fitCopula"
  )
)


gof <- vector("list", 2)
names(gof) <- c("low", "high")

for(i in seq_along(gof)) {
  gof[[i]] <- vector("list", (n * (n - 1) / 2))
  
  for(j in seq_along(gof[[i]])) {
    gof[[i]][[j]] <- cop_res(
      pairname = paste0(dimnames(dis_ls[[i]])[[2]], collapse = "-"),
      regime = names(gof)[i],
      u = dis_vola_ls[[i]][[j]],
      x = raw_dat_vola_ls[[i]][[j]],
      copula = best_cop[[i]][[j]][[1]]
    )
  }
}

gof_res <- vector("list", 2)
names(gof_res) <- c("low", "high")

for(i in seq_along(gof_res)) {
  gof_res[[i]] <- vector("list", (n * (n - 1) / 2))
  
  gof_res[[i]] <- parLapply(cl, gof[[i]], function(y) gof_wrap(y@copula@copula, y@x))
}

stopCluster(cl)

save(gof, gof_res, file = "./data/results/gof_res.RData")
