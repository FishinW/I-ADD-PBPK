# I-ADD PBPK model

This repository contains the R/mrgsolve implementation of the I-ADD PBPK model, release-profile fitting, model calibration/evaluation, local sensitivity analysis, and Monte Carlo simulations.

## Repository structure

- `R/models/`: mrgsolve model definitions.
- `R/analysis/`: release fitting, PBPK fitting/evaluation, sensitivity analysis, and Monte Carlo simulation scripts.
- `data/`: PK, device, and release-profile input data.

## Required R packages

```r
install.packages(c(
  "mrgsolve", "dplyr", "tidyr", "ggplot2", "purrr", "FME",
  "minpack.lm", "gridExtra", "patchwork", "data.table", "truncnorm"
))
```

## Suggested run order

0. `R/00_project_setup.R`
1. `R/analysis/01_fit_release_first_order.R`
2. `R/analysis/02_fit_release_weibull.R`
3. `R/analysis/03_fit_mouse_iv.R`
4. `R/analysis/04_fit_rat_iv.R`
5. `R/analysis/05_compare_mouse_rat_iv.R`
6. `R/analysis/06_simulate_rat_brain_iadd_mc.R`
7. `R/analysis/07_simulate_rat_kidney_iadd_mc.R`
8. `R/analysis/08_sensitivity_analysis.R`

`09_brain_iadd_discussion_scenarios.R` contains alternative brain I-ADD scenarios used for discussion-level comparisons.

## Reproducibility notes

- Run scripts from the repository root.
- Model equations and numerical parameter values were not intentionally changed during organization.
- Large analyses may take time because Monte Carlo simulations use repeated mrgsolve runs.
