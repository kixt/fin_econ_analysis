#####################################
# Financial Econometrics Term Paper #
#####################################

# Thore Petersen
# November 2019

## Tables for LaTeX

library(stargazer)
library(data.table)
library(moments)
library(countrycode)
library(stringr)
library(copula)


# Load data ---------------------------------------------------------------

load("./data/tmp/04_tmp.RData")


# Helper function ---------------------------------------------------------

caption_at_bottom <- function(x) {
  # take stargazer output and put the caption and label underneath the table
  # Args:
  #   x: character vector, output produced by stargazer
  # Returns:
  #   prints output like stargazer and invisibly returns the result
  
  cap <- grep("\\\\caption", x)
  lab <- grep("\\\\label", x)
  last <- grep("\\\\end\\{table", x)
  
  res <- paste(
    c(x[-last], x[cap], x[lab], x[last])[-c(cap, lab)], 
    collapse = "\n"
  )
  cat(res, "\n")
  invisible(res)
}


# Descriptives -- Returns -------------------------------------------------

returns_descript <- function(x) {
  # function to return basic descriptives for financial time series
  # Args:
  #   x: numeric vector of time data
  # Returns:
  #   a data.table with rows for mean, variance, skewness, kurtosis (all but 
  #   variance non-central), and the p-value of a Jarque-Bera test
  
  require(data.table)
  require(moments)
  
  res <- data.table(
    variable = c("Mean", "Variance", "Skewness", "Kurtosis", "Jarque-Bera p-value"),
    value = rep(NA_real_, 5)
  )
  
  res[, rownr := 1:.N]
  
  res[1:4, value := all.moments(x, 4)[-1]]
  res[variable == "Variance", value := moment(x, 2, central = TRUE)]
  res[variable == "Jarque-Bera p-value", value := jarque.test(x)$p.value]
  
  return(res)
}

descr_ret <- dt[, returns_descript(ret * 100), by = iso3c]
descr_ret[, country := countrycode(iso3c, "iso3c", "country.name")]
descr_ret2 <- dcast(descr_ret, variable + rownr ~ country, value.var = "value")
setorder(descr_ret2, rownr)
descr_ret2[, rownr := NULL]

descr_ret_sg <- stargazer(
  descr_ret2, 
  summary = FALSE, 
  label = "tab:data:ret_descript",
  rownames = FALSE,
  digits = 2
  )

descr_ret_sg <- caption_at_bottom(descr_ret_sg)

descr_ret_sg <- str_replace(
  descr_ret_sg, 
  "caption[{][}]", 
  paste("caption{Descriptive statistics for returns series in percent. ", 
        "The last row gives the $p$-value of a Jarque-Bera test for normality.}", 
        sep = "")
  )

descr_ret_sg <- str_replace(descr_ret_sg, "p-value", "$p$-value")
descr_ret_sg <- str_remove(descr_ret_sg, "variable")
cat(descr_ret_sg, sep = "\n")


# Cramer-von Mises marginals GoF ------------------------------------------

load("./data/tmp/cvm_gof_margins.RData")

cvm_res[, country := countrycode(iso3c, "iso3c", "country.name")]
cvm_res[, c("method", "iso3c") := NULL]
setcolorder(cvm_res, c(5, 3, 1, 2, 4))
setorder(cvm_res, country, tail)
colnames(cvm_res)[1:4] <- c("Country", "Tail", "Test statistic", "p-value")

cvm_gpd <- cvm_res[distr == "gpd"]
cvm_gpd[, distr := NULL]

cvm_gpd_sg <- stargazer(
  cvm_gpd, 
  summary = FALSE,
  label = "tab:margins_cvm",
  rownames = FALSE,
  digits = 3
  )

cvm_gpd_sg <- caption_at_bottom(cvm_gpd_sg)

cvm_gpd_sg <- str_replace(
  cvm_gpd_sg, 
  "caption[{][}]", 
  paste("caption{Results of Cramér-von Mises \\\\acs*{gof} tests of the fitted ", 
        "\\\\acs*{gpd} in the tails of the distributions.}", 
        sep = "")
)

cvm_gpd_sg <- str_replace(cvm_gpd_sg, "p-value", "$p$-value")
cvm_gpd_sg <- str_remove(cvm_gpd_sg, "Country")
cat(cvm_gpd_sg, sep = "\n")


# GoF results -------------------------------------------------------------

cop_res <- setClass(
  "cop_res", slots = list(
    pairname = "character", regime = "character", x = "matrix", u = "matrix", 
    copula = "fitCopula"
  )
)

load("./data/tmp/gof_res.RData")
load("./data/tmp/gof_av_res.RData")

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
      c("pi_F", "alpha", "p"), 
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

## asymmetric volatility
gof_tab_av <- rbindlist(
  lapply(gof_av_res, function(x) gof_res_row(x, NA_character_, NA_character_))
  )

gof_tab_h <- rbindlist(
  lapply(gof_res$high, function(x) gof_res_row(x, NA_character_, NA_character_))
  )

gof_tab_l <- rbindlist(
  lapply(gof_res$low, function(x) gof_res_row(x, NA_character_, NA_character_))
)

gof_tab_h[, pair := sapply(gof$high, function(x) x@pairname)]
gof_tab_h[, regime := "high-high"]

gof_tab_l[, pair := sapply(gof$low, function(x) x@pairname)]
gof_tab_l[, regime := "low-low"]

gof_tab_av[, pair := sapply(gof_av, function(x) x@pairname)]
gof_tab_av[, regime := sapply(gof_av, function(x) x@regime)]
# helper variables for sorting
gof_tab_av[, c1 := stringr::str_extract(pair, "^\\w{3}")]
gof_tab_av[, c2 := stringr::str_extract(pair, "\\w{3}$")]
gof_tab_av[, r1 := stringr::str_extract(regime, "high")]
gof_tab_av[, r2 := stringr::str_extract(regime, "low")]

# sort to have same pair, different regime in consecutive rows
done <- rep(NA_character_, length(unique(gof_tab_av$c1)))
for(c in unique(gof_tab_av$c1)) {
  gof_tab_av[
    !(c1 %in% done) & c2 == c, 
    c("c1", "c2", "r1", "r2") := .(c2, c1, r2, r1)
    ]
  done[which(is.na(done))[1]] <- c
}
setorder(gof_tab_av, c1, c2, r1)

gof_tab_av[, pair := paste0(c1, "-", c2)]
gof_tab_av[, k := ifelse(r1 == "high", 2, 3)]
gof_tab_av[, c("c1", "c2", "r1", "r2", "regime") := NULL]
setcolorder(gof_tab_av, c("pair", "k"))

gof_tab_h[, k := 4]
gof_tab_h[, regime := NULL]

gof_tab_l[, k := 1]
gof_tab_l[, regime := NULL]

gof_tab <- rbindlist(list(gof_tab_av, gof_tab_h), use.names = TRUE)
setorder(gof_tab, pair, k)

# all results in long for only DEU-FRA for presentation
gof_tab <- rbindlist(list(gof_tab_av, gof_tab_h, gof_tab_l), use.names = TRUE)
gof_tab <- gof_tab[pair == "DEU-FRA"]
gof_tab[, pair := NULL]
setorder(gof_tab, k)
gof_tab <- melt(gof_tab, id.vars = "k")
gof_tab <- dcast(gof_tab, variable ~ k)

# stargazer
gof_sg <- stargazer(
  gof_tab, 
  summary = FALSE,
  label = "tab:gof_res_av",
  rownames = FALSE,
  digits = 2
)

gof_sg <- caption_at_bottom(gof_sg)

gof_sg <- str_replace_all(gof_sg, "pi\\\\_(\\w{1})","$\\\\pi_\\1$")
gof_sg <- str_replace_all(gof_sg, "\\\\_","_")
gof_sg <- str_replace(gof_sg, "[&] k [&]","& $k$ &")
gof_sg <- str_replace(gof_sg, "pair", "Pair")
gof_sg <- str_replace(gof_sg, " p ", " $p$-value ")
gof_sg <- str_replace(gof_sg, "\\\\centering", "\\\\centering \\\\footnotesize")

for(s in c("alpha", "delta", "theta", "rho", "nu")) {
  gof_sg <- str_replace(gof_sg, s, paste0("$\\\\", s, "$"))
}

gof_sg <- str_replace(
  gof_sg, 
  "caption[{][}]", 
  paste("caption{Estimation results of best fitting copulas by market pair and", 
        " volatility regime. The table gives the parameter estimates of the ", 
        "single copulas (respective $\\\\pi$ equals unity), and for mixture ", 
        "copulas. \\\\acs*{gof} test $p$-values are from 1000 bootstrap ", 
        "replications of the $S_n^{(B)}$ test statistic of ", 
        "\\\\citet{Genest:2009}.}", 
        sep = "")
)

cat(gof_sg)


# Tail depedence ----------------------------------------------------------

# # no low vola regime
# td_tab <- data.table(
#   rbind(
#     data.table(
#       t(sapply(gof_av, function(x) c(pair = x@pairname, regime = x@regime))),
#       t(sapply(gof_av, function(x) lambda(x@copula@copula))),
#       sapply(gof_av, function(x) rho(x@copula@copula))
#     ),
#     data.table(
#       t(sapply(gof$high, function(x) c(pair = x@pairname, regime = x@regime))),
#       t(sapply(gof$high, function(x) lambda(x@copula@copula))),
#       sapply(gof$high, function(x) rho(x@copula@copula))
#     )
#   )
#   )
# 
# colnames(td_tab)[which(colnames(td_tab) == "V3")] <- "rho"
# 
# #td_tab[regime == "high", regime := "high-high"]
# # helper variables for sorting
# td_tab[, c1 := stringr::str_extract(pair, "^\\w{3}")]
# td_tab[, c2 := stringr::str_extract(pair, "\\w{3}$")]
# td_tab[, r1 := stringr::str_extract(regime, "high")]
# td_tab[, r2 := stringr::str_extract(regime, "low")]
# 
# # sort to have same pair, different regime in consecutive rows
# done <- rep(NA_character_, length(unique(td_tab$c1)))
# for(c in unique(td_tab$c1)) {
#   td_tab[
#     !(c1 %in% done) & c2 == c, 
#     c("c1", "c2", "r1", "r2") := .(c2, c1, r2, r1)
#     ]
#   done[which(is.na(done))[1]] <- c
# }
# setorder(td_tab, c1, c2, r1)
# 
# td_tab[, pair := paste0(c1, "-", c2)]
# td_tab[, k := ifelse(r2 == "low", 2, 3)]
# td_tab[, c("c1", "c2", "r1", "r2", "regime") := NULL]
# td_tab[is.na(k), k := 4]
# setcolorder(td_tab, c("pair", "k"))
# setorder(td_tab, pair, k)

# all regimes
td_tab <- data.table(
  rbind(
    data.table(
      t(sapply(gof_av, function(x) c(pair = x@pairname, regime = x@regime))),
      t(sapply(gof_av, function(x) lambda(x@copula@copula))),
      sapply(gof_av, function(x) rho(x@copula@copula))
    ),
    data.table(
      t(sapply(gof$high, function(x) c(pair = x@pairname, regime = x@regime))),
      t(sapply(gof$high, function(x) lambda(x@copula@copula))),
      sapply(gof$high, function(x) rho(x@copula@copula))
    ),
    data.table(
      t(sapply(gof$low, function(x) c(pair = x@pairname, regime = x@regime))),
      t(sapply(gof$low, function(x) lambda(x@copula@copula))),
      sapply(gof$low, function(x) rho(x@copula@copula))
    )
  )
)

colnames(td_tab)[which(colnames(td_tab) == "V3")] <- "rho"

# helper variables for sorting
td_tab[, c1 := stringr::str_extract(pair, "^\\w{3}")]
td_tab[, c2 := stringr::str_extract(pair, "\\w{3}$")]
td_tab[, r1 := stringr::str_extract(regime, "^\\w+.?(?=-)")] # match until "-"
td_tab[, r2 := stringr::str_extract(regime, "(?<=-)\\w+.$")] # match after "-"
td_tab[is.na(r1), c("r1", "r2") := stringr::str_extract(regime, "high|low")]

# sort to have same pair, different regime in consecutive rows
done <- rep(NA_character_, length(unique(td_tab$c1)))
for(c in unique(td_tab$c1)) {
  td_tab[
    !(c1 %in% done) & c2 == c, 
    c("c1", "c2", "r1", "r2") := .(c2, c1, r2, r1)
    ]
  done[which(is.na(done))[1]] <- c
}
setorder(td_tab, c1, c2, r1)

td_tab[, pair := paste0(c1, "-", c2)]
td_tab[r1 == "low" & r2 == "low", k := 1]
td_tab[r1 == "high" & r2 == "low", k := 2]
td_tab[r1 == "low" & r2 == "high", k := 3]
td_tab[r1 == "high" & r2 == "high", k := 4]
td_tab[, c("c1", "c2", "r1", "r2", "regime") := NULL]
setcolorder(td_tab, c("pair", "k"))
setorder(td_tab, pair, k)

# for presentation
td_tab <- td_tab[pair == "DEU-FRA"]
td_tab[, pair := NULL]
setorder(td_tab, k)
td_tab <- melt(td_tab, id.vars = "k")
td_tab <- dcast(td_tab, variable ~ k)

# stargazer
td_sg <- stargazer(
  td_tab, 
  summary = FALSE,
  label = "tab:cop_tail_dep",
  rownames = FALSE,
  digits = 3
)

td_sg <- caption_at_bottom(td_sg)

td_sg <- str_replace(td_sg, "pair", "Pair")
td_sg <- str_replace(td_sg, "\\\\centering", "\\\\centering")
td_sg <- str_replace(td_sg, " k ", " $k$ ")
td_sg <- str_replace(td_sg, " rho ", " Spearman's $\\\\rho$ ")
td_sg <- str_replace(td_sg, " lower ", " $\\\\lambda_l$ ")
td_sg <- str_replace(td_sg, " upper ", " $\\\\lambda_u$ ")

td_sg <- str_replace(
  td_sg, 
  "caption[{][}]", 
  paste("caption{Upper ($\\\\lambda_u$) and lower ($\\\\lambda_l$) tail ",
        "dependence, and Spearman's $\\\\rho$ implied by the fitted copulas.}",
        sep = "")
  )

cat(td_sg)










