Financial Econometrics Term Paper
Thore Petersen, Autumn/Winter 2019
NTNU Trondheim

01_clean_bench.R
- deprecated, use 01_clean.R instead -> use only total return indices
- cleans the indices Macrobond lists as benchmark equity indices per country, which is a mix between price and total return series

01_clean.R
- clean the total return indices of the country sample and export the data

02_garch.R
- imports data produced by 01_clean.R
- estimate the GARCH model to standardise the innovations
- exports the main dataset appended with the GARCH results

03_marginals.R
- soft-deprecated, use 03_marginals_gpd.R instead
- imports data produced by 02_garch.R
- fit parametric distribution (skew-t and skew-generalised error distribution) to the standardised residuals produced in 02_garch.R
- exports the main dataset appended with the fitted distribution function evaluated at the values of the GARCH residuals
- fit is pretty bad in tails, so use 03_marginals_gpd.R

03_marginals_gpd.R
- imports data produced by 02_garch.R
- fit generalised Pareto distribution to the tails of the GARCH residuals, use the ECDF in the middle
- fit is much better than that of skew-t and skew-GED in tails
- exports the main dataset appended with the concatenated distribution function evaluated at the values of the GARCH residuals

04_copula.R
- imports data produced by 03_marginals_gpd.R
- calculate and compare concordance measures on the values of the distribution function, and their empirical counterparts
- fit copulas...

05_visualise.R
- imports data produced by 03_marginals_gpd.R
- plot a heat map of pairwise distribution values on [0, 1]^2 (calculated as a 2D KDE)
- plot a heat map of a 2D KDE of the GARCH residuals
