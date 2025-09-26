# RLCs - a simple package to process rapid light curves

## Overview
This R package provides simple prcessing abilities for (sequential) rapid light curve datasets produced using a Walz WaterPAM fluorometer. It reads in a CSV file of RLC data, calculates key parameters (etr, npq, ynpq, yno), and fits the Eileers and Peters RLC model to the data to derive associated parameters (Fv/Fm, rETRmax, Ek, and alpha). If treatment identifier columns are provided, the package will take averages for each variable across treatment replicates and also fit a cumulative EnP model to the datasets, of predict_ul is set to TRUE, predictions of fit, upper and lower bounds are also provided for all models.

### Installation

This package can be directly installed from github using the devtools library in R, e.g.:

`library(devtools)`

`install_github("chrisjw18/rlcs")`

`library(rlcs)`

### Main Use (`fit_rlc()` function)

1. Make sure to curate your RLC dataset once downloaded as an excel sheet from the WaterPAM - i.e. separate out into columns, remove any rows that are not data rows, provide more friendly column names, deal with any missing values (leaving an empty cell is best), and save as a .CSV file.
3. Read in your CSV as a named object into R using standard CSV loaders e.g. `df1 <- read.csv('path/to/my/file.csv')`
5. Use the `rlc_fit()` function to automate processing, this requires the following arguments
- `dataframe` the name of the object containing the data
- `par` character string of the column name in your dataframe that contains the par levels - NB, this function can only process RLCs that have the same light levels
- `yield` character string of the column name in your dataframe that contains the yield values
- `f` character string of the column name in your dataframe that contains the f values
- `fm` character string of the column name in your dataframe that contains the fm values
- `sample` character string of the column name in your dataframe that contains the individual light curve sample names (should be unique per light curve)
- `treatment` character string of the column name in your dataframe that contains the treatment levels of your dataset if applicable (should be unique per set of light curves). If provided, averages are made per treatment level.
- `predict_ul` verbose (T/F), whether to provide rpedictions of fit, upper and lower bounds for e and p models

### Outputs

Once run, the function will provide the following outputs as a list named `rlc_results`:
1. `rlcs` containing the original data, with calculated variables (rETR, NPQ, YNPQ, YNO) and e and p model parameters (a,b,c) and derived parameters (rETRmax, alpha, Ek)
If `predict_ul == T`:
2. `fits` containing predicted fit, upper and lower bounds for each light curve e and p model calculated using predictNLS function
If `treatment` is provided:
3. `av_variables` containing averages of rETR, NPQ, YNPQ and YNO per treatment, including standard deviation and standard error
4. `av_params` containing averages of rETRmax, alpha, ek and fvfm per treatment, including standard deviation and standard error
5. `av_fits` containing e and p fit of etr ~ par with upper and lower bounds per treatment

### Further functions

1. `plot_rlcs()` allows to plot multiple rlcs on one plot, including errors - takes outputs of `rlc_fit()`, see `?plot_rlcs` in R
2. `plot_fits()` allows to plot multiple EnP fits with upper/lower bounds on one plot - takes outputs of `rlc_fit()`, see `?plot_fits` in R
3. `plot_params()` allows to plot average parmeters (Fv/Fm, rETRmax, alpha, Ek) ± SE - takes outputs of `rlc_fit()`, see `?plot_params` in R


