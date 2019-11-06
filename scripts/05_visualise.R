#####################################
# Financial Econometrics Term Paper #
#####################################

# Thore Petersen
# October, November 2019

library(data.table)
library(ggplot2)
library(MASS)
library(viridis)


# Load data ---------------------------------------------------------------

load("./data/tmp/03_tmp.RData")


# Heatmap of empirical quantiles ------------------------------------------

# calculate 2D KDEs for all country pairs, returns density estimates over grid 
# to plot in heatmap

kde_ls <- vector("list", n * (n - 1) / 2)
for(i in 1:(n - 1)) {
  for(j in (i + 1):n) {
    m1 <- markets[i]
    m2 <- markets[j]
    
    # reshape data and trim to pairwise complete observations
    dt_sub <- dcast(dt[.(c(m1, m2))], Date ~ iso3c, value.var = "qntl")
    dt_sub <- zoo::na.trim(dt_sub)
    
    # estimate bivariate KDE
    # need to include limits to make grid evenly spaced for geom_tiles()
    kde_tmp <- kde2d(dt_sub[, get(m1)], dt_sub[, get(m2)], n = 100, 
                     lims = c(0, 1, 0, 1))
    
    dt_tmp <- expand.grid(x = kde_tmp$x, y = kde_tmp$y)
    dt_tmp$z <- as.numeric(kde_tmp$z) # unwraps matrix
    
    dt_tmp$pair <- paste0(m1, "-", m2) # name of country pair
    
    # index for storage in list
    if(i == 1 && j == 2) {
      ls_i <- 1
    } else {
      ls_i <- ls_i + 1
    }
    
    kde_ls[[ls_i]] <- dt_tmp
  }
}; rm(m1, m2, i, j, ls_i, dt_tmp, dt_sub, kde_tmp)

eq_kde_dens <- rbindlist(kde_ls)

ggplot(eq_kde_dens, aes(x, y, fill = z)) +
  geom_tile() +
  scale_fill_viridis() +
  theme_minimal() +
  facet_wrap(pair ~ .)


# Heatmap of joint density estimate ---------------------------------------

support <- c(-8, 8)

kde_ls <- vector("list", n * (n - 1) / 2)
for(i in 1:(n - 1)) {
  for(j in (i + 1):n) {
    m1 <- markets[i]
    m2 <- markets[j]
    
    # reshape data and trim to pairwise complete observations
    dt_sub <- dcast(dt[.(c(m1, m2))], Date ~ iso3c, value.var = "res")
    dt_sub <- zoo::na.trim(dt_sub)
    
    # estimate bivariate KDE
    # need to include limits to make grid evenly spaced for geom_tiles()
    kde_tmp <- kde2d(dt_sub[, get(m1)], dt_sub[, get(m2)], n = 100, 
                     lims = c(support, support))
    
    dt_tmp <- expand.grid(x = kde_tmp$x, y = kde_tmp$y)
    dt_tmp$z <- as.numeric(kde_tmp$z) # unwraps matrix
    
    dt_tmp$pair <- paste0(m1, "-", m2) # name of country pair
    
    # index for storage in list
    if(i == 1 && j == 2) {
      ls_i <- 1
    } else {
      ls_i <- ls_i + 1
    }
    
    kde_ls[[ls_i]] <- dt_tmp
  }
}; rm(m1, m2, i, j, ls_i, dt_tmp, dt_sub, kde_tmp)

kde_dens <- rbindlist(kde_ls)

ggplot(kde_dens, aes(x, y, fill = z)) +
  geom_tile() +
  scale_fill_viridis() +
  theme_minimal() +
  facet_wrap(pair ~ .)

# not very instructive, cut off super heavy tails in limits for KDE, otherwise 
# plots hardly legible; other than that high variation per country, so dim 


