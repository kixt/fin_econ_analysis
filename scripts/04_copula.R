#####################################
# Financial Econometrics Term Paper #
#####################################

# Thore Petersen
# October, November 2019

library(data.table)
library(ggplot2)
library(copula)


# Load data ---------------------------------------------------------------

load("./data/tmp/03_tmp.RData")


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


# DEU -- FRA pair ---------------------------------------------------------

# characteristic pattern for DEU--FRA, almost linear, with clustering in corners,
# i.e. strong tail dependence
ggplot(eqntls, aes(x = USA, y = DEU)) + 
  geom_point()

# extract pairwise data for DEU and FRA, reshape/convert to numeric matrix
dt_DF <- rbindlist(list(dt["DEU"], dt["FRA"]))
c_dt_DF <- dcast(dt_DF, Date ~ iso3c, value.var = "qntl")
rm(dt_DF)
c_dt_DF[, Date := NULL]
c_dt_DF <- as.matrix(c_dt_DF)

c_dt <- dcast(dt, Date ~ iso3c, value.var = "qntl")
c_dt[, Date := NULL]
c_dt <- as.matrix(c_dt)

test <- fitCopula(copula = tCopula(), data = c_dt_DF)

fitLambda(c_dt) - fitLambda(c_dt, lower.tail = F)
