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


# Concordance measures ----------------------------------------------------

conc_measures <- c("pearson", "kendall", "spearman")
conc <- vector("list", 3)
names(conc) <- c("all", "normal vola", "high vola")

for(i in seq_along(conc)) {
  conc[[i]] <- vector("list", length(conc_measures))
  names(conc[[i]]) <- conc_measures
}

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
}

conc[[1]]
conc[[2]]
conc[[3]]
# higher pairwise concordances in high volatility regimes


# Empirical quantiles -----------------------------------------------------

# create data.table of empirical quantiles to check dependence
eqntls <- dcast(dt, Date ~ iso3c, value.var = "eqntl")
dens <- dcast(dt, Date ~ iso3c, value.var = "qntl")

dens[, Date := NULL]
dens <- as.matrix(dens)


# DEU -- FRA pair ---------------------------------------------------------

ecop1 <- empCopula(dens[, c("DEU", "FRA")], smoothing = "beta")
ecop1_plot <- wireframe2(ecop1, FUN = dCopula, screen = list(z = 10, x = -60))
print(ecop1_plot)

tcop <- fitCopula(tCopula(), dens[, c("DEU", "FRA")], method = "mpl")

r <- rCopula(10, tcop@copula)

wireframe2(tcop@copula, FUN = dCopula, screen = list(z = 10, x = -65))
