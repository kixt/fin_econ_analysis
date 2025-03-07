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


# Concordance measures ----------------------------------------------------

conc_measures <- c("pearson", "kendall", "spearman")
conc <- vector("list", 3)
names(conc) <- c("all", "normal vola", "high vola")

for(i in seq_along(conc)) {
  conc[[i]] <- vector("list", length(conc_measures))
  names(conc[[i]]) <- conc_measures
}; rm(i)

Fhs <- dcast(dt, Date ~ iso3c, value.var = "Fh")

for(c in conc_measures) {
  conc[[1]][[c]] <- cor(Fhs[, -1], 
                        method = c,
                        use = "pairwise.complete.obs")
}
# Pearson and Spearman super close, Kendall a little smaller for all

## concordance measures by volatility regime
# very simple regime indicator
dt[, high_vola := ifelse(sig > quantile(sig, 0.9), TRUE, FALSE), by = iso3c]

Fhs_hv <- dcast(dt[high_vola == TRUE], Date ~ iso3c, value.var = "Fh")
Fhs_lv <- dcast(dt[high_vola == FALSE], Date ~ iso3c, value.var = "Fh")

for(c in conc_measures) {
  conc[[2]][[c]] <- cor(Fhs_lv[, -1], 
                        method = c,
                        use = "pairwise.complete.obs")
}

for(c in conc_measures) {
  conc[[3]][[c]] <- cor(Fhs_hv[, -1], 
                        method = c,
                        use = "pairwise.complete.obs")
}; rm(c)

conc[[1]]
conc[[2]]
conc[[3]]
# higher pairwise concordances in high volatility regimes


# Check symmetry of tail dependences --------------------------------------

# non-parametric estimators of pairwise dependence parameters
dep_l <- fitLambda(dis)
dep_u <- fitLambda(dis, lower.tail = FALSE)

(dep_l - dep_u) / dep_l
# only DEU-IRL, FRA-IRL, FRA-ITA with not too large relative differences in 
# lower/upper tail dependence, although both low

dep_hv <- array(dim = c(n, n, 2))
dep_lv <- array(dim = c(n, n, 2))

for(i in 1:(n - 1)) {
  for(j in (i + 1):n) {
    cols <- c(i, j) + 1
    # high volatility
    u <- as.matrix(Fhs_hv[, ..cols])
    u <- na.exclude(u)
    dep_hv[i, j, 1] <- dep_hv[j, i, 1] <- fitLambda(u)[1, 2]
    
    # low volatility
    u <- as.matrix(Fhs_lv[, ..cols])
    u <- na.exclude(u)
    dep_lv[i, j, 1] <- dep_lv[j, i, 1] <- fitLambda(u)[1, 2]
  }
}; rm(i, j, u, cols)

for(i in 1:(n - 1)) {
  for(j in (i + 1):n) {
    cols <- c(i, j) + 1
    # high volatility
    u <- as.matrix(Fhs_hv[, ..cols])
    u <- na.exclude(u)
    dep_hv[i, j, 2] <- dep_hv[j, i, 2] <- fitLambda(u, lower.tail = FALSE)[1, 2]
    
    # low volatility
    u <- as.matrix(Fhs_lv[, ..cols])
    u <- na.exclude(u)
    dep_lv[i, j, 2] <- dep_lv[j, i, 2] <- fitLambda(u, lower.tail = FALSE)[1, 2]
  }
}; rm(i, j, cols, u)

# lower non-parametric tail dependence measure in low-volatility regime compared 
# to high volatility -> good result, reproduces literature
# in high-vola regime typically higher upper tail dependence than lower tail
# all tail dependence values still much larger than implied by asymmetric copulae

