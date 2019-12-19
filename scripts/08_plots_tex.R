#####################################
# Financial Econometrics Term Paper #
#####################################

# Thore Petersen
# November 2019

## Plots for LaTeX

library(copula)
library(ggplot2)
library(data.table)
library(MASS)
library(viridis)


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

# paper
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

# presentation
cop_fams <- ggplot(u, aes(x = x1, y = x2, group = fam)) +
  geom_point() + 
  #geom_density_2d() +
  facet_wrap(~fam, nrow = 1) +
  theme_minimal() +
  xlim(-4, 4) +
  ylim(-4, 4)

ggsave("../tex_pres/figures/cop_fams.pdf",
       cop_fams,
       device = "pdf",
       width = full_wdth, height = 3.5, units = "cm",
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


# DEU-FRA motivating example ----------------------------------------------

dt[, high_vola := ifelse(sig > quantile(sig, 0.9), TRUE, FALSE), by = iso3c]

pdt <- dt[c("DEU", "FRA")]

dt_kde <- vector("list", 4)
for(i in 1:4) {
  pdt_tmp <- switch(
    i,
    "1" = rbindlist(list(pdt["DEU"][high_vola == FALSE],
                         pdt["FRA"][high_vola == FALSE])),
    "2" = rbindlist(list(pdt["DEU"][high_vola == TRUE],
                         pdt["FRA"][high_vola == FALSE])),
    "3" = rbindlist(list(pdt["DEU"][high_vola == FALSE],
                         pdt["FRA"][high_vola == TRUE])),
    "4" = rbindlist(list(pdt["DEU"][high_vola == TRUE],
                         pdt["FRA"][high_vola == TRUE]))
    )
  pdt_tmp2 <- dcast(pdt_tmp, Date ~ iso3c, value.var = "Fh")
  pdt_tmp2 <- na.exclude(pdt_tmp2)
  
  kde_tmp <- kde2d(pdt_tmp2[, DEU], pdt_tmp2[, FRA], n = 100, 
                   lims = c(0, 1, 0, 1))
  
  dt_tmp <- expand.grid(x = kde_tmp$x, y = kde_tmp$y)
  dt_tmp$z <- as.numeric(kde_tmp$z)
  dt_tmp$k <- i
  dt_kde[[i]] <- dt_tmp
}

dt_kde <- rbindlist(dt_kde)

empc_hm_sq <- ggplot(dt_kde, aes(x, y, fill = z)) +
  geom_tile() +
  scale_fill_viridis(name = "Density") +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  theme_minimal() +
  facet_wrap(~k, ncol = 2, labeller = as_labeller(function(x) paste0("k = ", x))) +
  xlab(expression(u["DEU ,t"])) +
  ylab(expression(u["FRA ,t"])) +
  theme(panel.spacing.x = unit(4, "mm"))

ggsave("../tex_pres/figures/deu_fra_hm_sq.pdf",
       empc_hm_sq,
       device = "pdf",
       width = full_wdth, height = 8, units = "cm",
       scale = plot_scale)

empc_hm_line <- ggplot(dt_kde, aes(x, y, fill = z)) +
  geom_tile() +
  scale_fill_viridis(name = "Density") +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  theme_minimal() +
  facet_wrap(~k, ncol = 4, labeller = as_labeller(function(x) paste0("k = ", x))) +
  xlab(expression(u["DEU ,t"])) +
  ylab(expression(u["FRA ,t"])) +
  theme(panel.spacing.x = unit(4, "mm"))

ggsave("../tex_pres/figures/deu_fra_hm_line.pdf",
       empc_hm_line,
       device = "pdf",
       width = full_wdth, height = 3, units = "cm",
       scale = plot_scale)
