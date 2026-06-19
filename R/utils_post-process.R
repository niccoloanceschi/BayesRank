
#' Post-process Gibbs sampler output after fitting on the C. elegans isotonic data
#'
#' Discards burn-in draws and derives quantities of interest from the raw MCMC
#' output of `fit_BayesRank`, intended for use by downstream plotting functions.
#'
#' @param fit_mcmc Output list from `fit_BayesRank`.
#' @param input_data List containing the original model inputs.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{alpha}{Post-burn-in incremental treatment coefficients. ("Low" and "High vs. Low").}
#'     \item{beta}{Post-burn-in cumulative treatment effects (with "High vs. Low" columns converted to "High")}
#'     \item{mu}{Post-burn-in replicate intercepts, including the constrained last replicate (sum-to-zero).}
#'     \item{eta}{Post-burn-in random effects.}
#'     \item{median_RE}{Posterior median of each random effect.}
#'     \item{icc}{Posterior draws of the intraclass correlation, rho2/(1+rho2).}
#'     \item{Zy}{Post-burn-in latent utilities, split by ordinal category.}
#'     \item{Zmin}{Per-draw minimum latent utility within each category.}
#'     \item{Zmax}{Per-draw maximum latent utility within each category.}
#'   }
#'
#' @note Assumes the specific column-name convention ("Low"/"High vs. Low,
#'   <Generation>, <Un-/Rechallenged>") produced by the paper's design matrix (i.e.
#'   not a general-purpose postprocessing routine).
#'
#' @export
#' 
postprocess_BayesRank <- function(fit_mcmc, input_data){
  
  n_burn <- fit_mcmc$hyper_params$iter_Burn
  
  alpha <- fit_mcmc$alpha[-c(1:n_burn),]

  # Compute cumulative effects
  beta <- alpha
  
  # Rename "High vs. Low" columns to "High"
  colnames(beta) <- gsub("High vs. Low", "High", colnames(beta))
  
  # Now overwrite the "High" columns with cumulative sums
  beta[, "High, F1, Un-Rechallenged"] <- alpha[, "Low, F1, Un-Rechallenged"] +
    alpha[, "High vs. Low, F1, Un-Rechallenged"]
  beta[, "High, P0, Un-Rechallenged"] <- alpha[, "Low, P0, Un-Rechallenged"] +
    alpha[, "High vs. Low, P0, Un-Rechallenged"]
  beta[, "High, F1, Rechallenged"]    <- alpha[, "Low, F1, Rechallenged"]    +
    alpha[, "High vs. Low, F1, Rechallenged"]
  beta[, "High, P0, Rechallenged"]    <- alpha[, "Low, P0, Rechallenged"]    +
    alpha[, "High vs. Low, P0, Rechallenged"]
  
  mu <- fit_mcmc$mu[-c(1:n_burn),]
  mu <- cbind(mu,-rowSums(mu))
  colnames(mu) <- paste('Replicate',c(1:ncol(mu)))
  
  icc <- fit_mcmc$rho2[-c(1:n_burn)]/(1+fit_mcmc$rho2[-c(1:n_burn)])
  
  eta <- fit_mcmc$eta[-c(1:n_burn),]
  
  median_RE <- apply(eta, 2, median) 
  
  ZZ <- fit_mcmc$Z[-c(1:n_burn),]
  
  Zy = lapply(1:max(input_data$y), function(ss) ZZ[,which(input_data$y==ss)])
  
  Zmin <- sapply(Zy, function(mat) apply(mat,1,min)) # category-wise minimum
  Zmax <- sapply(Zy, function(mat) apply(mat,1,max)) # category-wise maximum

  return(list(alpha=alpha,beta=beta, mu=mu, eta=eta, median_RE=median_RE,
              icc=icc, Zy=Zy, Zmin=Zmin, Zmax=Zmax))
}

#' Predict ordinal damage scores from posterior draws of the Bayesian rank-based isotonic model
#'
#' Computes posterior predictive draws of the ordinal response category by
#' combining posterior draws of treatment effects, replicate intercepts,
#' and random effects into a linear predictor. The resulting continuous score
#' is then mapped it to an ordinal category via the posterior estimates of the
#' cutpoints, estimated from the inferred separation between latent utility.
#'
#' @param X Design matrix for treatment coefficients.
#' @param W Design matrix for group-specific intercepts.
#' @param Q Incidence matrix assigning observations to random-effect groups.
#' @param vals_mcmc Post-processed MCMC output from `postprocess_bayes`
#'
#' @return A list with elements:
#'   \describe{
#'     \item{y_mcmc}{Matrix of predicted categories, one row per
#'       observation and one column per post-burn-in MCMC draw.}
#'     \item{y_median}{Vector of median predicted category per observation,
#'       taken across MCMC draws.}
#'   }
#'
#' @export
#' 
predict_BayesRank <- function(X, W, Q, vals_mcmc){
  
  delta_hat <- colMeans(vals_mcmc$Zmax[, -ncol(vals_mcmc$Zmax)])
  
  Xsp <- as(X, "sparseMatrix")
  Wsp <- as(W, "sparseMatrix")
  Qsp <- as(Q, "sparseMatrix")
  
  X_treat = (Xsp %*% t(vals_mcmc$alpha))
  W_mu    = (Wsp %*% t(vals_mcmc$mu))
  Q_eta   = (Qsp %*% t(vals_mcmc$eta))
  
  Z_mean <- X_treat + W_mu + Q_eta
  
  Y_pred_mcmc <- matrix(
    findInterval(Z_mean, delta_hat) + 1,
    nrow = nrow(Z_mean),
    ncol = ncol(Z_mean)
  )
  
  Y_pred <- apply(Y_pred_mcmc, 1, median)
  
  return(list(y_mcmc=Y_pred_mcmc,y_median=Y_pred))
}

#' Predict ordinal class from a fitted cumulative link mixed model
#'
#' Computes predicted ordinal categories from a fitted `clmm` model by
#' combining the fixed-effect linear predictor with subject-specific
#' random-effect, then mapping the result to a category via the
#' model's estimated cutpoints.
#'
#' @param fit_clmm A fitted `clmm` object (from package `ordinal`)
#' @param df_ord Data used to fit `clmm` frame (or one with the same format
#'
#' @return A factor of predicted categories, with levels matching
#'   `fit_clmm$y.levels`.
#'
#' @export
#' 
predict_clmm <- function(fit_clmm, df_ord){
  
  # 1. Extract model components
  beta <- coef(fit_clmm)              # fixed-effect coefficients (named vector)
  # Remove threshold parameters from coef (clmm stores them together)
  
  n_thresh <- length(fit_clmm$alpha)  # number of cutpoints
  beta_fixed <- beta[-(1:n_thresh)]   # the actual fixed-effect coefs
  
  cuts <- fit_clmm$alpha              # cutpoints
  
  re_blups <- ranef(fit_clmm)$UniqueWormID[, 1]
  names(re_blups) <- rownames(ranef(fit_clmm)$UniqueWormID)
  
  # 2. Build design matrix matching the model formula
  # Use the same columns that went into the model
  X_pred <- as.matrix(df_ord[, names(beta_fixed)])
  
  # 3. Compute linear predictor (with subject-specific RE)
  eta <- as.vector(X_pred %*% beta_fixed) + re_blups[as.character(df_ord$UniqueWormID)]
  
  # 4. Map to predicted class via cutpoints
  # clmm parameterization: P(Y <= ell) = link^{-1}(cuts[ell] - eta)
  # So predicted class = smallest ell such that eta < cuts[ell], or L if eta >= cuts[L-1]
  pred_class <- findInterval(eta, cuts) + 1
  pred_class <- factor(pred_class, levels = 1:length(fit_clmm$y.levels))
  
  return(y_pred=pred_class)
}


