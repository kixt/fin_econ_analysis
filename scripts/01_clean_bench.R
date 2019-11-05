#####################################
# Financial Econometrics Term Paper #
#####################################

# Thore Petersen
# November 2019

library(openxlsx)
library(data.table)
library(zoo)


# Import and prepare data -------------------------------------------------

mb_imp <- read.xlsx("./data/raw/macrobond.xlsx", sheet = "mb_bench", 
                    na.strings = "NaN", detectDates = TRUE)
setDT(mb_imp)
mb_tmp <- melt(mb_imp, id.vars = "Date")
mb_tmp[, Date := as.Date(Date, "%Y-%m-%d")]

# replace . in variable names with whitespace for easier regex matching
mb_tmp[, variable := gsub("[..]", " ", variable)]
# extract country information
mb_tmp[, iso3c := countrycode::countrycode(variable, "country.name", "iso3c")]
# extract currency information
mb_tmp[, currency := stringr::str_extract(variable, "EUR|USD")]
setkey(mb_tmp, iso3c, Date)

mb_tmp["DEU", index := "DAX30"]
mb_tmp["FRA", index := "CAC40"]
mb_tmp["IRL", index := "ISEQ"]
mb_tmp["ITA", index := "MIB"]
mb_tmp["USA", index := "DowJones"]

# calculate returns
mb_tmp[, ret := c(NA_real_, diff(log(value))), by = variable]

# trim leading missings
mb_tmp2 <- mb_tmp[, na.trim(.SD), by = variable]
setkey(mb_tmp2, iso3c, Date)


# Subset to relevant series -----------------------------------------------

mb_fin <- mb_tmp2[currency %in% c("EUR", "USD")]

mb_fin[, c("variable", "value") := NULL]
mb_fin[, c("currency") := NULL]

# Date of official EURO introduction and fixation of intra-European exch. rates
mb_fin <- mb_fin[Date >= as.Date("1999-01-01")]

dt <- mb_fin
save(dt, file = "./data/clean/clean.RData")