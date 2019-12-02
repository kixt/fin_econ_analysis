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
# load("./data/tmp/best_cop.RData")

cop_res <- setClass(
  "cop_res", slots = list(
    pairname = "character", regime = "character", x = "matrix", u = "matrix", 
    copula = "fitCopula"
  )
)

load("./data/tmp/gof_res.RData")
load("./data/tmp/gof_av_res.RData")


# Goodness of Fit ---------------------------------------------------------


lapply(gof$low, function(x) x@copula@fitting.stats$convergence)

# lapply(best_cop$high, function(x) x[[1]]@fitting.stats$convergence)



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

lapply(gof_res, function(x) lapply(x, function(y) y$p.value))
lapply(gof_av_res, function(y) y$p.value)


# Prepare table for LaTeX -------------------------------------------------

gof_res_row <- function(htest, pair, regime) {
  # function to generate row for LaTeX results table of GoF results
  # Args:
  #   htest: object of class htest containing the GoF results
  #   pair: character, name of country pair
  #   regime: character, volatility regime
  # Returns:
  #   data.table containing the results
  
  c <- stringr::str_extract(
    htest$method, 
    "t-copula|Gumbel|Frank|Clayton|Frnc, Gmbc, Clyc|Frnc, t-cp, Gmbc, Clyc"
    )
  if(!(c %in% c("Gumbel", "Frank", "Clayton"))) {
    c <- switch(
      c, 
      "t-copula" = "t",
      "Frnc, Gmbc, Clyc" = "m_nt",
      "Frnc, t-cp, Gmbc, Clyc" = "m_t"
    ) 
  }
  
  res <- data.table(
    pair = pair, regime = regime, pi_F = NA_real_, pi_G = NA_real_, 
    pi_C = NA_real_, pi_t = NA_real_, alpha = NA_real_, delta = NA_real_, 
    theta = NA_real_, rho = NA_real_, nu = NA_real_, p = NA_real_
  )
  
  res <- melt(res, id.vars = c("pair", "regime"))
  setkey(res, variable)
  
  if(c == "t") {
    res[
      c("pi_t", "rho", "nu", "p"), 
      value := c(1, htest$parameter, htest$p.value)
      ]
  } else if(c == "Gumbel") {
    res[
      c("pi_G", "delta", "p"), 
      value := c(1, htest$parameter, htest$p.value)
      ]
  } else if(c == "Clayton") {
    res[
      c("pi_C", "theta", "p"), 
      value := c(1, htest$parameter, htest$p.value)
      ]
  } else if(c == "Frank") {
    res[
      c("pi_G", "alpha", "p"), 
      value := c(1, htest$parameter, htest$p.value)
      ]
  } else if(c == "m_nt") {
    res[
      c("alpha", "delta", "theta", "pi_F", "pi_G", "pi_C", "p"), 
      value := c(htest$parameter, htest$p.value)
      ]
  } else if(c == "m_t") {
    res[
      c("alpha", "rho", "nu", "delta", "theta", "pi_F", "pi_t", "pi_G", "pi_C", "p"), 
      value := c(htest$parameter, htest$p.value)
      ]
  }
  
  # set all other weights of non-mixture copulas to zero
  res[grepl("pi+", variable) & is.na(value), value := 0]
  
  res <- dcast(res, pair + regime ~ variable)
  return(res)
}

res_av <- rbindlist(lapply(gof_av_res, function(x) gof_res_row(x, NA_character_, NA_character_)))

res_av[, pair := sapply(gof_av, function(x) x@pairname)]
res_av[, regime := sapply(gof_av, function(x) x@regime)]

res_av[, c1 := stringr::str_extract(pair, "^\\w{3}")]
res_av[, c2 := stringr::str_extract(pair, "\\w{3}$")]

res_av[, r1 := stringr::str_extract(regime, "high")]
res_av[, r2 := stringr::str_extract(regime, "low")]

done <- rep(NA_character_, length(unique(res_av$c1)))
for(c in unique(res_av$c1)) {
  res_av[!(c1 %in% done) & c2 == c, c("c1", "c2", "r1", "r2") := .(c2, c1, r2, r1)]
  done[which(is.na(done))[1]] <- c
}

setorder(res_av, c1, c2, r1)

res_av[, pair := paste(c1, c2, collapse = "-")]
res_av[, k := ifelse(r1 == "high", 2, 3)]

res_av[, c("c1", "c2", "r1", "r2", "regime") := NULL]
setcolorder(res_av, c("pair", "k"))
