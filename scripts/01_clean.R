#####################################
# Financial Econometrics Term Paper #
#####################################

# Thore Petersen
# September 2019

library(openxlsx)
library(data.table)
library(zoo)


# Import and prepare data -------------------------------------------------

mb_imp <- read.xlsx("./data/macrobond.xlsx", 
                    na.strings = "NaN", detectDates = TRUE)
setDT(mb_imp)
mb_tmp <- melt(mb_imp, id.vars = "Date")

# replace . in variable names with whitespace for easier regex matching
mb_tmp[, variable := gsub("[.]", " ", variable)]
# extract country information
mb_tmp[, iso3c := countrycode::countrycode(variable, "country.name", "iso3c")]
# extract currency information
mb_tmp[, currency := stringr::str_extract(variable, "EUR|GBP|USD|NOK")]

# order and calculate returns
setorder(mb_tmp, variable, Date)
mb_tmp[, ret := c(NA_real_, diff(log(value))), by = variable]

# trim leading missings
mb_tmp2 <- mb_tmp[, na.trim(.SD), by = variable]


# Subset to relevant series -----------------------------------------------

mb_fin <- mb_tmp2[currency != "GBP"]
mb_fin[, c("variable", "value") := NULL]

# Date of official EURO introduction and fixation of intra-European exch. rates
mb_fin <- mb_fin[Date >= as.Date("1999-01-01")]

# GBR data for EUR denominated index starts only mid 2005; GRC starts only early 2003

dt <- mb_fin
save(dt, file = "./data/clean.RData")
