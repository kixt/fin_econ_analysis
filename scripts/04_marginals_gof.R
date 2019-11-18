#####################################
# Financial Econometrics Term Paper #
#####################################

# Thore Petersen
# November 2019

## Compare goodness of fit for the different specification of marginals
# determine best fit and export data

library(data.table)
library(goftest)


# Load data ---------------------------------------------------------------

load("./data/tmp/03_tmp_gpd.RData")
load("./data/tmp/03_tmp_stsged.RData")


# Calculate MSE -----------------------------------------------------------

colnames(dt_gpd)[which(colnames(dt_gpd) == "Fh")] <- "Fh_gpd"
colnames(dt_stsged)[which(colnames(dt_stsged) == "Fh")] <- "Fh_stsged"

dt <- merge(dt_gpd, dt_stsged, 
            by = intersect(colnames(dt_gpd), colnames(dt_stsged)))
setkey(dt, iso3c, Date)
rm(dt_gpd, dt_stsged)

dt[, c("mse_gpd", "mse_stsged") := .(mean((eFh - Fh_gpd)^2), 
                                     mean((eFh - Fh_stsged)^2)), 
   by = .(iso3c, tail)]

# MSE in lower tail always lower for skew-t/skew-GED compared to GPD; in upper tail
# GPD typically better
dt[, mean(mse_gpd - mse_stsged), by = .(iso3c, tail)]


# Cramér-von Mises test ---------------------------------------------------

# apply CvM test in the tails (only there, otherwise mixture of ecdf part in 
# mixed GPD fit kind of cheating)

dt[tail == "lower", cvm.test(Fh_gpd / q), by = .(iso3c)]
dt[tail == "lower", cvm.test(Fh_stsged / q), by = .(iso3c)]

dt[tail == "upper", cvm.test((Fh_gpd - 1) / q + 1), by = .(iso3c)]
dt[tail == "upper", cvm.test((Fh_stsged - 1) / q + 1), by = .(iso3c)]


# Export data -------------------------------------------------------------

dt[, c("Fh_stsged", "mse_gpd", "mse_stsged") := NULL]
colnames(dt)[which(colnames(dt) == "Fh_gpd")] <- "Fh"

save(dt, markets, n, file = "./data/tmp/04_tmp.RData")

