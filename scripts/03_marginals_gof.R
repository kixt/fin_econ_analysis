#####################################
# Financial Econometrics Term Paper #
#####################################

# Thore Petersen
# November 2019

## Compare goodness of fit for the different specification of marginals

library(data.table)


# Load data ---------------------------------------------------------------

load("./data/tmp/03_tmp_gpd.RData")
load("./data/tmp/03_tmp_stsged.RData")


# Calculate MSE -----------------------------------------------------------

colnames(dt_gpd)[which(colnames(dt_gpd) == "qntl")] <- "d_gpd"
colnames(dt_stsged)[which(colnames(dt_stsged) == "qntl")] <- "d_stsged"

dt <- merge(dt_gpd, dt_stsged, 
            by = intersect(colnames(dt_gpd), colnames(dt_stsged)))
setkey(dt, iso3c, Date)
rm(dt_gpd, dt_stsged)

dt[, c("mse_gpd", "mse_stsged") := .(mean((eqntl - d_gpd)^2), 
                                     mean((eqntl - d_stsged)^2)), 
   by = .(iso3c, tail)]
