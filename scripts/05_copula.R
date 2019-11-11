#####################################
# Financial Econometrics Term Paper #
#####################################

# Thore Petersen
# October, November 2019

library(data.table)
library(ggplot2)
library(copula)


# Load data ---------------------------------------------------------------

load("./data/tmp/04_tmp.RData")


# Empirical quantiles -----------------------------------------------------

# create matrix of pseudo observations to use in copula functions
dens <- dcast(dt, Date ~ iso3c, value.var = "qntl")

dens[, Date := NULL]
dens <- as.matrix(dens)


# Concordance measures ----------------------------------------------------

conc_measures <- c("pearson", "kendall", "spearman")
conc <- vector("list", 3)
names(conc) <- c("all", "normal vola", "high vola")

for(i in seq_along(conc)) {
  conc[[i]] <- vector("list", length(conc_measures))
  names(conc[[i]]) <- conc_measures
}; rm(i)

qntls <- dcast(dt, Date ~ iso3c, value.var = "qntl")

for(c in conc_measures) {
  conc[[1]][[c]] <- cor(qntls[, -1], 
                        method = c,
                        use = "pairwise.complete.obs")
}
# Pearson and Spearman super close, Kendall a little smaller for all

## concordance measures by volatility regime
# very simple regime indicator
dt[, high_vola := ifelse(sig > quantile(sig, 0.9), TRUE, FALSE), by = iso3c]

qntls_hv <- dcast(dt[high_vola == TRUE], Date ~ iso3c, value.var = "qntl")
qntls_lv <- dcast(dt[high_vola == FALSE], Date ~ iso3c, value.var = "qntl")

for(c in conc_measures) {
  conc[[2]][[c]] <- cor(qntls_lv[, -1], 
                        method = c,
                        use = "pairwise.complete.obs")
}

for(c in conc_measures) {
  conc[[3]][[c]] <- cor(qntls_hv[, -1], 
                        method = c,
                        use = "pairwise.complete.obs")
}; rm(c)

conc[[1]]
conc[[2]]
conc[[3]]
# higher pairwise concordances in high volatility regimes


# Check symmetry of tail dependences --------------------------------------

# non-parametric estimators of pairwise dependence parameters
dep_l <- fitLambda(dens)
dep_u <- fitLambda(dens, lower.tail = FALSE)

(dep_l - dep_u) / dep_l
# only DEU-IRL, FRA-IRL, FRA-ITA with not too large relative differences in 
# lower/upper tail dependence, although both low

dep_hv <- array(dim = c(n, n, 2))
dep_lv <- array(dim = c(n, n, 2))

for(i in 1:(n - 1)) {
  for(j in (i + 1):n) {
    cols <- c(i, j) + 1
    # high volatility
    u <- as.matrix(qntls_hv[, ..cols])
    u <- na.exclude(u)
    dep_hv[i, j, 1] <- dep_hv[j, i, 1] <- fitLambda(u)[1, 2]
    
    # low volatility
    u <- as.matrix(qntls_lv[, ..cols])
    u <- na.exclude(u)
    dep_lv[i, j, 1] <- dep_lv[j, i, 1] <- fitLambda(u)[1, 2]
  }
}; rm(i, j, u, cols)

for(i in 1:(n - 1)) {
  for(j in (i + 1):n) {
    cols <- c(i, j) + 1
    # high volatility
    u <- as.matrix(qntls_hv[, ..cols])
    u <- na.exclude(u)
    dep_hv[i, j, 2] <- dep_hv[j, i, 2] <- fitLambda(u, lower.tail = FALSE)[1, 2]
    
    # low volatility
    u <- as.matrix(qntls_lv[, ..cols])
    u <- na.exclude(u)
    dep_lv[i, j, 2] <- dep_lv[j, i, 2] <- fitLambda(u, lower.tail = FALSE)[1, 2]
  }
}; rm(i, j, cols, u)

# lower non-parametric tail dependence measure in low-volatility regime compared 
# to high volatility -> good result, reproduces literature
# in high-vola regime typically higher upper tail dependence than lower tail
# all tail dependence values still much larger than implied by asymmetric copulae


# DEU -- FRA pair ---------------------------------------------------------

ecop1 <- empCopula(dens[, c("DEU", "FRA")], smoothing = "beta")
ecop1_plot <- wireframe2(ecop1, FUN = dCopula, screen = list(z = 10, x = -60), n.grid = 100)
print(ecop1_plot)

tcop <- fitCopula(tCopula(), dens[, c("DEU", "USA")], method = "mpl")

wireframe2(tcop@copula, FUN = dCopula, screen = list(z = 10, x = -65))

mixcop <- mixCopula(list(claytonCopula(1),
                         gumbelCopula(1.5),
                         frankCopula(1)))

test <- fitCopula(mixcop, dens[, c("DEU", "USA")], method = "mpl")

wireframe2(test@copula, FUN = dCopula, screen = list(z = 5, x = -70))


# Estimate mixture copulae ------------------------------------------------

# laod data from HPC calculation
load("./data/tmp/cop_mix_fit.RData")


lapply(mc_fit, function(x) x[[2]]@copula@w)



