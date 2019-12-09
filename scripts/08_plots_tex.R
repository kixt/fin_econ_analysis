#####################################
# Financial Econometrics Term Paper #
#####################################

# Thore Petersen
# November 2019

## Plots for LaTeX

library(copula)
library(ggplot2)
library(data.table)


# Load data ---------------------------------------------------------------

load("./data/tmp/04_tmp.RData")

plot_scale <- 1.5
ts_line_size <- 0.1
full_wdth <- 10


# Labeller functions ------------------------------------------------------

ret_lab <- function(x) {
  
  #if(x == "sig") return("Volatility")
  #if(x == "ret") return("Return")
  
  iso3c <- vector("character", length(x))
  for(i in seq_along(x)) {
    iso3c[i] <- dt[index == x[i], iso3c][1]
  }
  
  if(any(x == "SP500")) x[which(x == "SP500")] <- "S&P500"
  paste0(countrycode::countrycode(iso3c, "iso3c", "country.name"), " (", x, ")")
}

lab_ret_vol <- function(x) {
  res <- x
  for(i in which(!(res %in% c("ret", "sig")))) {
    res[i] <- dt[index == x[i], iso3c][1]
    if(x[i] == "SP500") x[i] <- "S&P500"
    res[i] <- paste0(countrycode::countrycode(res[i], "iso3c", "country.name"), 
                     " (", x[i], ")")
  }
  for(i in which(res == "ret")) res[i] <- "Return"
  for(i in which(res == "sig")) res[i] <- "Volatility"
  return(res)
}


# Copula plots ------------------------------------------------------------

copulas <- c("Frank", "Clayton", "Gumbel", "t")
u <- vector("list", length(copulas))
n <- 1000

set.seed(314)
for(i in seq_along(u)) {
  cop <- switch(copulas[i],
                "Frank" = frankCopula(5),
                "Clayton" = claytonCopula(2),
                "Gumbel" = gumbelCopula(2),
                "t" = tCopula(0.8, 2))
  u[[i]] <- data.table(rCopula(n, cop))
  u[[i]][, fam := copulas[i]]
}

u <- rbindlist(u)
setcolorder(u, "fam")
colnames(u)[2:3] <- c("u1", "u2")

u[, c("x1", "x2") := .(qnorm(u1), qnorm(u2))]

cop_fams <- ggplot(u, aes(x = x1, y = x2, group = fam)) +
  geom_point() + 
  #geom_density_2d() +
  facet_wrap(~fam) +
  theme_minimal() +
  xlim(-4, 4) +
  ylim(-4, 4)

ggsave("../tex/figures/cop_fams.pdf",
       cop_fams,
       device = "pdf",
       width = full_wdth, height = 10, units = "cm",
       scale = plot_scale)


# Returns series ----------------------------------------------------------

rets <- ggplot(dt, aes(x = Date, y = ret, group = index)) +
  geom_line(size = ts_line_size) +
  facet_wrap(~index, ncol = 1, labeller = as_labeller(ret_lab)) +
  theme_minimal() +
  ylab("Return") +
  theme(axis.text.y = element_text(size = 6))

ggsave("../tex/figures/returns_ts.pdf",
       rets,
       device = "pdf",
       width = full_wdth, height = 7, units = "cm",
       scale = plot_scale)


# Volatility series -------------------------------------------------------

vols <- ggplot(dt, aes(x = Date, y = sig, group = index)) +
  geom_line(size = ts_line_size) +
  facet_wrap(~index, ncol = 1, labeller = as_labeller(ret_lab)) +
  theme_minimal() +
  ylab("Volatility") +
  theme(axis.text.y = element_text(size = 6))

ggsave("../tex/figures/volatilities_ts.pdf",
       vols,
       device = "pdf",
       width = full_wdth, height = 7, units = "cm",
       scale = plot_scale)


# Volatility and returns --------------------------------------------------

pdt_ret_vol <- melt(dt, 
                    id.vars = c("index", "Date"), 
                    measure.vars = c("ret", "sig"))

vols_rets <- ggplot(pdt_ret_vol, aes(x = Date, y = value, group = index)) +
  geom_line(size = ts_line_size) +
  facet_wrap(index~variable, 
             ncol = 1, 
             scales = "free_y", 
             labeller = as_labeller(lab_ret_vol, multi_line = FALSE)) +
  theme_minimal() +
  ylab(NULL) +
  theme(axis.text.y = element_text(size = 6)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3))

# for paper
ggsave("../tex/figures/vola_ret_ts.pdf",
       vols_rets,
       device = "pdf",
       width = full_wdth, height = 10, units = "cm",
       scale = plot_scale)

# for presentation
ggsave("../tex_pres/figures/vola_ret_ts.pdf",
       vols_rets,
       device = "pdf",
       width = full_wdth, height = 7.5, units = "cm",
       scale = 2)


# Non-uniform GARCH residuals ---------------------------------------------

load("./data/tmp/03_tmp_stsged.RData")

marg_hist <- ggplot(dt_stsged, aes(x = Fh, group = index)) +
  geom_histogram(bins = 50) +
  theme_minimal() +
  facet_wrap(~index, labeller = as_labeller(ret_lab)) +
  xlab("values of fitted CDF") +
  #xlab(expression(hat("F")^p~(epsilon~"*"))) +
  ylab(NULL)

ggsave("../tex/figures/st_sged_marg_hist.pdf",
       marg_hist,
       device = "pdf",
       width = full_wdth, height = 5, units = "cm",
       scale = plot_scale)
       

# GARCH residuals ---------------------------------------------------------

resids <- ggplot(dt, aes(x = Date, y = res, group = index)) +
  geom_line(size = ts_line_size) +
  facet_wrap(~index, ncol = 1, labeller = as_labeller(ret_lab)) +
  theme_minimal() +
  ylab("GARCH residuals") +
  theme(axis.text.y = element_text(size = 6))

ggsave("../tex/figures/garch_resid_ts.pdf",
       resids,
       device = "pdf",
       width = full_wdth, height = 6, units = "cm",
       scale = plot_scale)

