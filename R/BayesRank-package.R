#' BayesRank: Bayesian Rank-Based Isotonic Regression for Caenorhabditis elegans Assays
#' 
#' Implements a Bayesian isotonic rank-likelihood model for ordinal regression
#' under monotonicity constraints, developed to analyze Caenorhabditis elegans neuron
#' damage scores from a developmental rotenone exposure assay from
#' [Bergemann et al. (2026) "\emph{Progeny effects of rotenone exposure 
#'  depend on parental toxicity}"](https://academic.oup.com/toxsci/article/209/3/kfag011/8475413).
#' The package reproduces the analysis from [Presman, Anceschi, Huayta,
#' Meyer, and Herring, (2026+) "\emph{Order-Restricted Bayesian Ordinal Regression for
#' the Modeling of Neuron Degeneration in Caenorhabditis
#' elegans}"](https://arxiv.org/abs/PLACEHOLDER).
#'
#' Relative to existing rank-likelihood implementations, it additionally supports
#' random effects, parameter constraints, collapsed sampling, and parameter expansion.
#' While the bundled data, design matrices, and tutorial are specific to
#' the C. elegans assay design described above, the underlying isotonic
#' rank-likelihood methodology readily generalizes to a broad range of rank-based
#' regression settings.
#'
#' @keywords internal
#' 
#' @useDynLib BayesRank, .registration = TRUE
#' @import dplyr
#' @import ggplot2
#' @import tidybayes
#' @import ggh4x
#' @import patchwork
#' @import viridis
#' @import Matrix
#' @importFrom tidyr pivot_longer separate
#' @importFrom scales breaks_width
#' @importFrom MASS mvrnorm
#' @importFrom RConics cubic
#' @importFrom ordinal clmm
#' @importFrom tibble tibble as_tibble
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
