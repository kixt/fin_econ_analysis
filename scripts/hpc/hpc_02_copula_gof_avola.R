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
load("./data/tmp/best_cop_avola.RData")


# Create cluster ----------------------------------------------------------

cl <- makeCluster(3)
clusterEvalQ(cl, .libPaths(new = "./R_libs/libs"))
clusterEvalQ(cl, library("copula")) # load copula package at nodes


# Prepare data subsets ----------------------------------------------------

# very simple regime indicator
dt[, high_vola := ifelse(sig > quantile(sig, 0.9), TRUE, FALSE), by = iso3c]

dis_ls <- vector("list", n * (n - 1))
raw_dat_ls <- vector("list", n * (n - 1))

ls_i <- 1
for(i in 1:n) {
  for(j in 1:n) {
    if(i != j) {
      ic <- markets[i]
      jc <- markets[j]
      # Fh
      tmp1 <- dt[ic][high_vola == TRUE, .(Date, Fh)]
      tmp2 <- dt[jc][high_vola == FALSE, .(Date, Fh)]
      tmpm <- merge(tmp1, tmp2, by = "Date", all = FALSE)
      colnames(tmpm)[2:3] <- c(ic, jc)
      dis_ls[[ls_i]] <- as.matrix(tmpm[, c(2, 3)])
      
      # raw data
      tmp1 <- dt[ic][high_vola == TRUE, .(Date, res)]
      tmp2 <- dt[jc][high_vola == FALSE, .(Date, res)]
      tmpm <- merge(tmp1, tmp2, by = "Date", all = FALSE)
      colnames(tmpm)[2:3] <- c(ic, jc)
      raw_dat_ls[[ls_i]] <- as.matrix(tmpm[, c(2, 3)])
      
      ls_i <- ls_i + 1
    }
  }
}; rm(tmp1, tmp2, tmpm, i, ic, jc, j, ls_i)


# GoF for general model ---------------------------------------------------

gof_wrap <- function(copula, x) {
  tryCatch(
    gofCopula(copula = copula, x = x, N = 1000, method = "SnB", 
              test.method = "single"),
    error = function(e) e
  )
}
clusterExport(cl, "gof_wrap")


cop_res_av <- setClass(
  "cop_res", slots = list(
    pairname = "character", regime = "character", x = "matrix", u = "matrix", 
    copula = "fitCopula"
  )
)

gof_av <- vector("list", n * (n - 1))

for(i in seq_along(gof_av)) {
  gof_av[[i]] <- cop_res_av(
    pairname = paste0(dimnames(dis_ls[[i]])[[2]], collapse = "-"),
    regime = "high-low",
    u = dis_ls[[i]],
    x = raw_dat_ls[[i]],
    copula = best_cop_avola[[i]][[1]]
  )
}

gof_av_res <- parLapply(cl, gof_av, function(y) gof_wrap(y@copula@copula, y@x))

stopCluster(cl)

save(gof_av, gof_av_res, file = "./data/results/gof_av_res.RData")
