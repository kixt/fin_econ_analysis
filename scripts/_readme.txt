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
- imports data produced by 02_garch.R
- fit parametric distribution (skew-t and skew-generalised error distribution) to the standardised residuals produced in 02_garch.R
- exports 03_tmp_stsged.RData

03_marginals_gpd.R
- imports data produced by 02_garch.R
- fit generalised Pareto distribution to the tails of the GARCH residuals, use the ECDF in the middle
- fit is much better than that of skew-t and skew-GED in tails
- exports 03_tmp_gpd.RData

04_marginals_gof.R
- imports 03_tmp_stsged.RData and 03_tmp_gpd.RData
- calculate MSE per tail, compare by tail
- exports dataset with distribution function evaluated at the standardised residuals

05_copula.R
- imports data produced by 04_marginals_gof.R
- calculate and compare concordance measures on the values of the distribution function, and their empirical counterparts
- fit copulas...

06_visualise.R
- imports data produced by 04_marginals_gof.R
- plot a heat map of pairwise distribution values on [0, 1]^2 (calculated as a 2D KDE)
- plot a heat map of a 2D KDE of the GARCH residuals
