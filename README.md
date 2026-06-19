# BayesRank

**Bayesian rank-based isotonic regression for the analysis of Caenorhabditis elegans toxicological assays**

`BayesRank` is an `R` package implementing the analysis from
[Presman, Anceschi, Huayta, Meyer, and Herring, (2026+) "*Order-Restricted Bayesian Ordinal Regression for the Modeling of Neuron Degeneration in Caenorhabditis elegans*"](https://arxiv.org/abs/PLACEHOLDER),
and allowing to reproduce the main results from the manuscript.

It fits a Bayesian rank-based model for ordinal regression under monotonicity constraints, extendsing existing rank-likelihood software to support random effects, parameter constraints, collapsed sampling, and parameter expansion.
While the package's code is specific to the target assay design, the implemented modeling framework readily generalizes to a much broader range of rank-based isotonic regression settings.

## Installation

The `BayesRank` package can be installed by running the following `R` commands

```r
# If the devtools R package is not already installed
# install.packages("devtools")

devtools::install_github("niccoloanceschi/BayesRank")
```

Alternatively, if `devtools::install_github()` produces warnings or fails, the `pak` package can be used:

```r
# If the pak R package is not already installed
# install.packages("pak")

pak::pak("niccoloanceschi/BayesRank")
```
## Documentation and Tutorial

The main functions are detailed in the package [manual](https://github.com/niccoloanceschi/BayesRank/blob/main/BayesRank_0.1.0.pdf).

A full [tutorial](https://niccoloanceschi.github.io/BayesRank/tutorial/C-elegans_Data_Analysis.html) with detailed workflow is provided in the repository.

## Example usage 

```{r}
library(BayesRank)

# Load raw data
data <- get_data()

# Format data to run Gibbs sampler
data_bayes <- get_data_BayesRank(data)

# Fit model
fit_bayes <- fit_BayesRank(data_bayes$y, data_bayes$X, data_bayes$W, data_bayes$Q,
                           iter_MC=10000, iter_Burn=5000)

# Process posterior samples
draws_mcmc <- postprocess_BayesRank(fit_bayes,data_bayes)

# Generate predictions
y_pred <- predict_BayesRank(data_bayes$X, data_bayes$W, data_bayes$Q, draws_mcmc)
```

## Structure

- `src/` contains the `Rcpp` code scripts
- `R/` contains the `R` source scripts
- `man/` contains the documentation files for the exported `R` functions (generated via `roxygen2`)
- `data/` contains the original assay data (in `.rda` format)
- `tutorial/` contains scripts to reproduce the main results and figures from the main paper
- `DESCRIPTION` is the package metadata file specifying authors, maintainers, version, and dependencies
- `NAMESPACE` defines both the imported and exported functions
- `LICENSE` contains the licensing terms of the package (MIT license)
- `.Rbuildignore` specifies files and directories that should be excluded when building the package
- `.gitignore `specifies files and directories ignored by Git version control
- `BayesRank.Rproj` is the RStudio project file for development and reproducibility

## Reference

[Presman, Anceschi, Huayta, Meyer, and Herring, (2026+) "*Order-Restricted Bayesian Ordinal Regression for the Modeling of Neuron Degeneration in Caenorhabditis elegans*"](https://arxiv.org/abs/PLACEHOLDER)
