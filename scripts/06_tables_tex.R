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

