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


# Goodness of Fit ---------------------------------------------------------

load("./data/tmp/best_cop.RData")

cop_res <- setClass(
  "cop_res", slots = list(
    pairname = "character", regime = "character", x = "matrix", u = "matrix", 
    copula = "fitCopula"
  )
)

load("./data/tmp/gof_res.RData")


lapply(gof$low, function(x) x@copula@fitting.stats$convergence)

lapply(best_cop$low, function(x) x[[1]]@fitting.stats$convergence)



for(m1 in markets) {
  for(m2 in markets) {
    for(v in c(T, F)) {
      print(dt[, 
         length(intersect(
           .SD[which(iso3c == m1 & high_vola == v), Date], 
           .SD[which(iso3c == m2 & high_vola == v), Date]
           ))
         ])
    }
  }
}

