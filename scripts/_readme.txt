Financial Econometrics Term Paper
Thore Petersen, Autumn/Winter 2019
NTNU Trondheim

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

05_dependence.R
- imports data produced by 04_marginals_gof.R
- calculate and interpret concordance measures on the values of the distribution function, and their empirical counterparts

06_prep_copula.R
- imports data produced by 04_marginals_gof.R, hpc/hpc_01_fit_copula.R, and hpc/hpc_01_fit_copula_asym_vola.R on High Performance Computing Cluster
- check convergence and determine best fitting copula by AIC
- exports best_cop.RData and best_cop_avola.RData

07_copula_gof.R
- imports data produced by 04_marginals_gof.R and hpc/hpc_02_copula_gof.R, and hpc/hpc_02_copula_gof_avola.R
- check GoF results
- produce table source code for LaTeX

08_plots_tex.R
- imports various datasets
- exports .pdf files of plots for paper and presentation

08_tables_tex.R
- imports various datasets
- export LaTeX source code of tables for paper

08_visualise.R
- imports data produced by 04_marginals_gof.R
- plot a heat map of pairwise distribution values on [0, 1]^2 (calculated as a 2D KDE)
- plot a heat map of a 2D KDE of the GARCH residuals

/hpc/
- contains scripts to estimate models on High Performance Computing Cluster

