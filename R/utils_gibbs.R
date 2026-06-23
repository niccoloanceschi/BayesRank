
#' Find the previous non-empty group of observations (by observed category)
#'
#' @param Jh List of index vectors, keyed by category.
#' @param start Index to start searching from (searching backward).
#' @param min_idx Index to stop at / return if nothing found.
#' 
#' @return Integer index of the previous non-empty element, or min_idx.
#' 
#' @noRd
#' 
prev_nonempty <- function(Jh, start, min_idx) {
  for(i in seq(start, min_idx)) if(length(Jh[[i]]) > 0) return(i)
  return(min_idx)
}

#' Find the next non-empty group of observations (by observed category)
#'
#' @param Jh List of index vectors, keyed by category.
#' @param start Index to start searching from.
#' @param max_idx Index to stop at / return if nothing found.
#' 
#' @return Integer index of the next non-empty element, or max_idx.
#' 
#' @noRd
#' 
next_nonempty <- function(Jh, start, max_idx) {
  for(i in seq(start, max_idx)) if(length(Jh[[i]]) > 0) return(i)
  return(max_idx)
}

#' Vectorized sampling from univariate truncated normal distributions
#' (inverse-CDF method)
#'
#' @param n Number of draws.
#' @param mean Mean of the underlying normal. Default 0.
#' @param sd SD of the underlying normal. Default 1.
#' @param lb Lower truncation bound. Default -Inf.
#' @param ub Upper truncation bound. Default Inf.
#' 
#' @return Numeric vector of length n.
#' 
#' @noRd
#' 
rtruncnorm_plain <- function(n, mean=0, sd=1, lb=-Inf, ub=Inf) {
  mean + sd * qnorm(runif(n, pnorm((lb - mean) / sd), pnorm((ub - mean) / sd)))
}

#' Vectorized sampling from multivariate truncated normal distributions
#'
#' Updates a block of correlated truncated-normal latent utilities,
#' by updating one coordinate at a time via inner Gibbs sampling sweeps.
#' Optionally dispatches to an Rcpp implementation when available for speed.
#'
#' @param Z_m Matrix of current latent utility values for the group.
#' @param Z_mean_m Matrix of conditional means.
#' @param muC_m Conditioning adjustment term from the marginalized random effect.
#' @param rho2_m Random-effect variance term for this group size.
#' @param sd_m Conditional SD for the truncated normal draws.
#' @param lb,ub Truncation bounds (ordinal category cutoffs).
#' @param n_iter Number of inner Gibbs sweeps. Default 10.
#' @param rcpp Use the Rcpp implementation? Default TRUE.
#' 
#' @return Updated matrix of latent utility values, same shape as Z_m.
#' 
#' @noRd
#' 
rMVTN_gibbs <- function(Z_m, Z_mean_m, muC_m, rho2_m, sd_m, lb, ub, n_iter=10, rcpp=T){
  
  if(rcpp){
    Z_m <- sample_MVTN_gibbs_cpp(Z_m=Z_m, Z_mean_m=Z_mean_m, muC_m=muC_m, 
                                 rho2_m=rho2_m, sd_m=sd_m, lb_h=lb, ub_h=ub, n_iter=n_iter)
  } else {
    m = nrow(Z_m)
    q = ncol(Z_m)
    
    Z_diff_m <- Z_m - Z_mean_m  # Compute once
    Z_sum_m <- colSums(Z_diff_m)  # Total sum
    rho2_muC_m = rho2_m * muC_m
    
    for(tt in 1:n_iter) {
      for(s in 1:m) {
        # Subtract current row from total to get sum of others
        muC_m <- Z_mean_m[s,] + rho2_muC_m + rho2_m * (Z_sum_m - Z_diff_m[s,])
        
        # Sample new value
        Z_m[s,] <- rtruncnorm_plain(q, mean=muC_m, sd=sd_m, lb=lb, ub=ub)
        
        # Update the running sum with the change
        Z_new_diff <- Z_m[s,] - Z_mean_m[s,]
        Z_sum_m <- Z_sum_m - Z_diff_m[s,] + Z_new_diff
        Z_diff_m[s,] <- Z_new_diff
      }
    }
  }
  return(Z_m)
}

#' Compute spike probability for a spike-and-slab prior update
#'
#' Computes the conditional posterior probability of the spike (point mass
#' at zero) component, given truncated-normal spike and slab parameters,
#' using a log-sum-exp for numerical stability.
#'
#' @param pi0 Prior spike probability.
#' @param mu0 Mean of the (truncated) spike component.
#' @param std0 SD of the spike component.
#' @param mu1 Mean of the (truncated) slab component.
#' @param std1 SD of the slab component.
#' 
#' @return Posterior spike probability.
#' 
#' @noRd
#' 
prob_spike <- function(pi0,mu0,std0,mu1,std1){
  
  log_pSpike <- log(pi0) + pnorm(mu0/std0, log.p = TRUE) -
    dnorm(0, mean = mu0, sd = std0, log = TRUE) 
  
  log_pSlab <- log(1 - pi0) + pnorm(mu1/std1, log.p = TRUE) -
    dnorm(0, mean = mu1, sd = std1, log = TRUE)
  
  log_pMax <- max(log_pSpike, log_pSlab)
  
  pSpike <- exp(log_pSpike - log_pMax)
  pSlab <- exp(log_pSlab - log_pMax)
  
  pi1 <- pSpike/(pSpike + pSlab)
  
  return(pi1)
}

#' Construct an asinh-warped grid for griddy Gibbs sampling of rho2
#'
#' Builds a grid of candidate values for the random-effect variance,
#' concentrated near the current mode (rho2_map) via an asinh transform.
#'
#' @param rho2_map Current mode/center of the grid.
#' @param rho2_max Upper bound of the grid.
#' @param n_grid Number of interior grid points. Default 100.
#' @param delta Spread parameter controlling concentration near rho2_map.
#'   Defaults to 10% of the distance to the nearer bound.
#'   
#' @return Numeric vector of grid values, with 0 as the first element.
#' 
#' @noRd
#' 
make_grid <- function(rho2_map, rho2_max, n_grid=100, delta = NULL) {
  stopifnot(rho2_map > 0, rho2_map < rho2_max, n_grid >= 4)
  
  half <- n_grid %/% 2
  if (is.null(delta)) delta <- 0.1 * min(rho2_map, rho2_max - rho2_map)
  
  phi_lo <- asinh(-rho2_map / delta)
  phi_hi <- asinh((rho2_max - rho2_map) / delta)
  
  phi_left  <- seq(phi_lo, 0, length.out = half + 1)[-1]
  phi_right <- seq(0, phi_hi, length.out = half + 1)[-1]
  
  phi  <- c(phi_left, phi_right)
  rho2 <- rho2_map + delta * sinh(phi)
  
  rho2 <- pmax(rho2, .Machine$double.eps)
  rho2 <- pmin(rho2, rho2_max)
  c(0,rho2)
}

#' Sample a value from a discrete grid via inverse-CDF interpolation
#'
#' @param r_vals Grid of candidate values.
#' @param p_vals Probability weights for each grid value.
#' 
#' @return A single interpolated draw from the grid.
#' 
#' @noRd
#' 
rCDFgrid <- function(r_vals, p_vals){
  
  cdf <- cumsum(p_vals)
  u   <- runif(1)
  k_star <- which(cdf[-1]>u)[1] 
  
  d_r <- (r_vals[k_star+1]-r_vals[k_star])
  d_cdf <- (u-cdf[k_star]) / (cdf[k_star+1]-cdf[k_star])
  
  r_vals[k_star] + d_r * d_cdf
}

#' Log unnormalized posterior of rho2 under random-effect marginalization
#'
#' @param rho2 Vector of candidate variance values.
#' @param sigma2 Residual variance.
#' @param a,b Inverse-gamma prior shape/rate for rho2.
#' @param v1,v2 Sufficient statistics by random-effect group size.
#' @param M Maximum random-effect group size.
#' 
#' @return Vector of log unnormalized posterior densities, one per rho2 value.
#' 
#' @noRd
#' 
logP_rho2_REmarg <- function(rho2, sigma2, a, b, v1, v2, M) {
  
  vec1M <- c(1:M)
  log_prior <- -(a + 1) * log(rho2) - b / rho2 # Prior: IG(a, b)
  r2_v2_s2 <- tcrossprod(vec1M, rho2) + sigma2
  log_det <- -0.5 * colSums(v2*log(r2_v2_s2)) # Determinant term
  log_exp <- -0.5 * colSums(v1 / (vec1M*r2_v2_s2)) # Quadratic form
  
  return(log_exp + log_det + log_prior) 
}

#' Draw rho2 from its marginal posterior via Griddy Gibbs
#'
#' @param r_vals Grid of candidate rho2 values.
#' @param s2 Residual variance.
#' @param a,b Inverse-gamma prior shape/rate for rho2.
#' @param v1,v2 Sufficient statistics by random-effect group size.
#' @param M Maximum random-effect group size.
#' 
#' @return A single sampled value of rho2.
#' 
#' @noRd
#' 
sample_rho2_REmarg <- function(r_vals, s2, a, b, v1, v2, M){
  
  logP_vals_REmarg <- c(-Inf, logP_rho2_REmarg(r_vals[-1], s2, a, b, v1, v2, M))
  
  p_vals_REmarg <- exp(logP_vals_REmarg-max(logP_vals_REmarg))
  p_vals_REmarg <- p_vals_REmarg/sum(p_vals_REmarg)
  
  rCDFgrid(r_vals,p_vals_REmarg) 
}

#' Find the saddle point for the Laplace approximation used in updating
#' the auxiliary variable in the parameter expansion for rho2
#' 
#' Solves a cubic equation for the positive real root used as the mode in
#' the saddle-point/Laplace approximation underlying the PX-DA update.
#'
#' @param r_shape,r_rate Inverse-gamma prior shape/rate for rho2.
#' @param sigma2 Residual variance.
#' @param R Number of random-effect groups.
#' @param mQ Mean group size.
#' @param cB Sufficient statistic (projected sum of squares).
#' 
#' @return The positive real root of the cubic (saddle point).
#' 
#' @noRd
#' 
get_saddle_point <- function(r_shape, r_rate, sigma2, R, mQ, cB){
  
  a3 <- mQ^2 * ((1+r_shape)+0.5*R)
  a2 <- mQ * (2*(1+r_rate)*sigma2 + 0.5*R*sigma2 - r_rate*mQ - 0.5*cB*sigma2)
  a1 <- sigma2 * ((1+r_rate) - 2*r_shape*mQ)
  a0 <- -r_rate * sigma2^2
  
  cubeS <- RConics::cubic(c(a3, a2, a1, a0))
  
  saddle <- Re(cubeS[abs(Im(cubeS)) < 1e-8 & Re(cubeS) > 0])
  
  return(saddle)
}

# Stan-style progress printer for Gibbs sampler iterations.
# Prints "Iteration: i / total [pct%] (Warmup|Sampling)" at ~10% intervals.
#
# @param: i Current iteration.
# @param: iter_MC Total iterations.
# @param: iter_Burn Warmup.
print_iter_progress <- function(i, iter_MC, iter_Burn) {
  pct <- floor(100 * i / iter_MC)
  phase <- if (i <= iter_Burn) "Warmup" else "Sampling"
  cat(sprintf("Iteration: %*d / %d [%3d%%]  (%s)\n",
              nchar(iter_MC), i, iter_MC, pct, phase))
}

#' Check and fill defaults for Gibbs sampler hyperparameters
#'
#' Validates a list of hyperparameters for `gibbs_sampler`, filling in
#' default values for any missing elements.
#'
#' @param hyper_params List of hyperparameters, with elements:
#'   \describe{
#'     \item{w0}{Prior SD for normla prior on the group-specific intercepts. Default 1.}
#'     \item{pi0}{Prior spike probability for treatment coefficients. Default 0.5.}
#'     \item{mu0}{Prior mean for the gaussian slab component of treatment coefficients. Default 0.}
#'     \item{lambda0}{Prior SD for the gaussian slab component of treatment coefficients. Default 1.}
#'     \item{r_shape}{Shape hyperparameter for the inverse-gamma prior on the random-effect variance (rho2). Default 2.}
#'     \item{r_rate}{Rate hyperparameter for the inverse-gamma prior on the random-effect variance (rho2). Default 2.}
#'     \item{iter_MH}{Inner sweeps for multivariate truncated normal Gibbs sampling. Default 10.}
#'     \item{rho2_GG}{Initial/center value for Griddy-Gibbs on rho2. Default 2.}
#'     \item{rho2_max}{Upper bound for Griddy-Gibbs on rho2. Default 5.}
#'     \item{n_grid_rho2}{Number of grid points for Griddy-Gibbs on rho2. Default 100.}
#'   }
#'   May be partially specified or NULL; missing elements are filled with defaults.
#'
#' @return A complete list of hyperparameters.
#'   
#' @export
#' 
set_hyperparameters <- function(hyper_params = NULL) {
  
  defaults <- list(
    w0 = 1, pi0 = 0.5, mu0 = 0, lambda0 = 1,
    r_shape = 2, r_rate = 2,
    iter_MH = 10, rho2_GG = 2, rho2_max = 5, n_grid_rho2 = 100
  )
  
  if (is.null(hyper_params)) hyper_params <- list()
  
  missing_names <- setdiff(names(defaults), names(hyper_params))
  hyper_params[missing_names] <- defaults[missing_names]
  
  hyper_params[names(defaults)]
}

#' Gibbs Sampler for Bayesian Isotonic Rank Likelihood Model
#'
#' Fits a Bayesian rank-based isotonic regression model via Gibbs sampling.
#' Optional analytic marginalization of random effects (after burn-in) and
#' parameter expansion for selected updates are available for improved mixing.
#'
#' @param y Ordinal response vector (values 1:K).
#' @param X Design matrix for treatment coefficients.
#' @param W Design matrix for group-specific intercepts.
#' @param Q Incidence matrix assigning observations to random-effect groups.
#' @param MC_seed Random seed. Default 0.
#' @param iter_MC Total MCMC iterations (must exceed iter_Burn). Default 10000.
#' @param iter_Burn Burn-in iterations. Default 5000.
#' @param hyper_params List of hyper-parameters (see `set_hyperparameters`)
#' @param marg_RE Marginalize random effects after burn-in? Default TRUE.
#' @param PX Use parameter expansion (PX-DA)? Default TRUE.
#' @param rcpp Use Rcpp for multivariate truncated normal sampling? (Only if `marg_RE` is TRUE) Default TRUE.
#'
#' @return List with posterior draws (including warmup) for 
#' \describe{
#'     \item{alpha}{Incremental treatment coefficients.}
#'     \item{beta}{Cumulative treatment effects.}
#'     \item{mu}{Replicate intercepts, excluding the constrained last replicate (sum-to-zero).}
#'     \item{eta}{Random effects.}
#'     \item{rho2}{Variance of random effects.}
#'     \item{Z}{Latent utilities of the rank likelihood.}
#'   }
#'  The list also includes hyper-parameters and MCMC settings.
#'  
#' @references
#' [Presman, Anceschi, Huayta, Meyer, and Herring, (2026+) "Order-Restricted
#' Bayesian Ordinal Regression for the Modeling of Neuron Degeneration in
#' Caenorhabditis elegans"](https://arxiv.org/pdf/2606.23358)
#' 
#' @export
#' 
fit_BayesRank <- function(y, X, W, Q, MC_seed=0, iter_MC=10000, iter_Burn=5000, 
                          hyper_params=NULL, marg_RE=T, PX=T, rcpp=T){
  
  if (iter_MC <= iter_Burn) {
    stop("Total iterations (iter_MC) must be greater than warmup (iter_Burn)")
  }
  
  # Set seed
  set.seed(MC_seed)
  
  # Unpack hyper-parameters 
  hyper_params <- set_hyperparameters(hyper_params)
  list2env(hyper_params, envir = environment())
  
  # Input processing -----------------------------------------------------------
  
  # Initialize dimensions 
  N <- length(y) # n. of observations
  K <- max(y) # n. of categories
  R <- ncol(Q) # n. of random effects
  C <- ncol(X) # n. of coeffs 
  G <- ncol(W) # n. of groups
  
  # modify group matrix W to enforce sum-0 intercepts
  W[, 1:(G-1)] <- W[, 1:(G-1)] - W[, G]
  W <- W[, 1:(G-1)]
  G <- G-1
  
  # Pre-compute worm & group counts
  nQ <- colSums(Q) # number of observations per RE group
  M = max(nQ) # max. number of observations per RE group
  HM <- 1 * t(outer(nQ, 1:M, "==")) # RE by group size
  nM <- rowSums(HM) # numerosity of RE group sizes
  mQ <- mean(nQ) # mean number of observations per RE group
  
  # Instantiate sparse matrices
  tQsp <- as(t(Q), "sparseMatrix")
  tXsp <- as(t(X), "sparseMatrix")
  Xsp <- as(X, "sparseMatrix")
  Wsp <- as(W, "sparseMatrix")
  Qsp <- as(Q, "sparseMatrix")
  Hsp <- as(HM, "sparseMatrix")
  
  # Initialize matrix products
  WtW <- t(W)%*%W
  QtW <- t(Q)%*%W
  QtX <- tQsp%*%X
  sumX2 <- colSums(X^2)
  
  # Index Pre-computation for Z-update -----------------------------------------
  
  # Z-update pre-computations: indices for bounds (shifted by 1)
  Jh = vector("list", K + 2)
  for(h in 2:(K+1)) Jh[[h]] = which(y == h-1)
  
  # Z-update pre-computations: indices of non-empty sets
  h_up  <- sapply(1:K, function(h) next_nonempty(Jh, h+2, K+2))
  h_low <- sapply(1:K, function(h) prev_nonempty(Jh, h,   1  ))
  
  if(marg_RE){
    
    # Initialize values for griddy gibbs
    rho2_vals = make_grid(rho2_GG,rho2_max,n_grid_rho2)
    
    # Composite index selection
    Fg = vector("list", R)
    for(g in 1:R) Fg[[g]] = which(Q[,g] == 1)
    
    qGroup <- matrix(0,K,M)
    
    # counting groups with m occurrences of y=h
    for(h in 1:K){
      for(g in 1:R){
        idx_len <- length(intersect(Fg[[g]], Jh[[h+1]]))
        if(idx_len>0){
          qGroup[h,idx_len] = qGroup[h,idx_len] + 1
        }
      }
    }
    
    # getting indices of groups with m occurrences of y=h
    Indx_hg <- Targ_hg <- Cond_hg <- vector("list", K)
    for(h in 1:K){
      Targ_hg[[h]] <- Cond_hg[[h]] <- Indx_hg[[h]] <- vector("list")
      for(m in 1:M){
        Indx_hg[[h]][[m]] = which(sapply(1:R, function(g) length(intersect(Fg[[g]], Jh[[h+1]]))==m))
        Targ_hg[[h]][[m]] = matrix(NA, nrow=m, ncol=qGroup[h,m])
        Cond_hg[[h]][[m]] = matrix(NA, nrow=M-m, ncol=qGroup[h,m])
      }
    }
    
    # filling matrices of target and conditioning indices with m occurrences of y=h by group
    for(h in 1:K){
      for(m in 1:M){
        if(qGroup[h,m]>0){
          idx_targ = sapply(Indx_hg[[h]][[m]], function(g) intersect(Fg[[g]], Jh[[h+1]]))
          Targ_hg[[h]][[m]] = idx_targ
          idx_cond = sapply(Indx_hg[[h]][[m]], function(g) {
            idx_hmg = setdiff(Fg[[g]], Jh[[h+1]]);
            c(idx_hmg,rep(N+1,M-m-length(idx_hmg)))})
          Cond_hg[[h]][[m]] = idx_cond
        }
      }
    }
    
  }
  
  # Initialize algorithm and outputs -------------------------------------------
  
  ## Z -- latent utilities
  cuts0 <- qnorm(cumsum(table(factor(y, levels = 1:K)))[-K]/N)
  lb0 <- c(-Inf, cuts0)
  ub0 <- c(cuts0,Inf)
  Z <- rep(0,N)
  for(h in 1:K){
    if(length(Jh[[h+1]]) > 0){
      Z[Jh[[h+1]]] <- rtruncnorm_plain(length(Jh[[h+1]]), mean=0, sd=1, lb=lb0[h], ub=ub0[h])
    }
  }
  
  sigma2 <- 1                                     # residual variance (fixed)
  s2inv <- 1 / sigma2                             # residual precision
  rho2 <- rho2_GG                                 # RE variance
  eta <- unname( t(Z%*%Q) / nQ )                  # random effects (RE) 
  beta <- unname( t(Z%*%W) / colSums(W>0) )       # group-specific intercepts
  treatment_coeffs <- matrix(0, nrow=C, ncol = 1) # coeffs
  
  # linear predictors
  X_treat = as.vector(Xsp %*% treatment_coeffs)
  W_beta = as.vector(Wsp %*% beta)
  Q_eta = as.vector(Qsp %*% eta)
  
  # Initialize outputs
  sigma2_vec <- numeric(iter_MC)
  rho2_vec <- numeric(iter_MC)
  beta_mat <- matrix(0, iter_MC, G)
  eta_mat <- matrix(0, iter_MC, R)
  Z_mat <- matrix(0, iter_MC, N)
  treatment_mat <- matrix(0, iter_MC, C)
  
  cat("SAMPLING FOR MODEL 'Bayesian Isotonic Rank Likelihood'\n\n")
  
  # Gibbs updates --------------------------------------------------------------
  
  for (i in 1:iter_MC) {
    
    # Print every 10% of total iterations, plus always the very first iteration
    if (i == 1 || i %% ceiling(iter_MC / 10) == 0) {
      print_iter_progress(i, iter_MC, iter_Burn)
    }
    
    # Regular updates : during burn-in or without RE marginalization
    regular_update = ( (i<iter_Burn) || (!marg_RE) )
    
    ## Z _ latent utilities ----
    if(regular_update){
      Z_mean <- X_treat + W_beta + Q_eta
      for(h in 1:K){
        if(length(Jh[[h+1]]) > 0){
          
          lb <- max(c(Z[Jh[[h_low[[h]]]]], -Inf), na.rm = TRUE)
          ub <- min(c(Z[Jh[[h_up[[h]]]]], +Inf), na.rm = TRUE)
          
          Z[Jh[[h+1]]] <- rtruncnorm_plain(length(Jh[[h+1]]),
            mean=Z_mean[Jh[[h+1]]], sd=sqrt(sigma2), lb=lb, ub=ub)
        }
      }
    } else {
      Z_mean <- c( X_treat + W_beta , 0 )
      Z0 <- c( Z , 0 )
      
      for(h in 1:K){
        if(length(Jh[[h+1]]) > 0){
          
          lb_h <- max(c(Z0[Jh[[h_low[[h]]]]], -Inf), na.rm = TRUE)
          ub_h <- min(c(Z0[Jh[[h_up[[h]]]]], +Inf), na.rm = TRUE)
          
          for(m in 1:M){  
            if(qGroup[h,m]>0){
              
              rho2_nQm <- rho2 / (sigma2 + rho2*(nQ[Indx_hg[[h]][[m]]]-1))
              sd_Qm <- sqrt(sigma2 + rho2_nQm * sigma2)
              
              if(m==1){
                muC_m <- Z_mean[Targ_hg[[h]][[m]]] + rho2_nQm *
                  colSums(matrix(Z0[Cond_hg[[h]][[m]]]-Z_mean[Cond_hg[[h]][[m]]],nrow=(M-m)))
                
                Z0[Targ_hg[[h]][[m]]] <- rtruncnorm_plain(qGroup[h,m], mean=muC_m, 
                                                          sd=sd_Qm, lb=lb_h, ub=ub_h) 
              } else {
                muC_m = rep(0,qGroup[h,m])
                if(m<M){ 
                  muC_m <- colSums(matrix(Z0[Cond_hg[[h]][[m]]]-Z_mean[Cond_hg[[h]][[m]]],nrow=(M-m))) 
                }
                
                Z_m <- matrix(Z0[Targ_hg[[h]][[m]]],nrow=m)
                Z_mean_m <- matrix(Z_mean[Targ_hg[[h]][[m]]],nrow=m)
                
                Z0[Targ_hg[[h]][[m]]] <- rMVTN_gibbs(Z_m, Z_mean_m, muC_m, 
                  rho2_nQm, sd_Qm, lb_h, ub_h, n_iter=iter_MH, rcpp=rcpp)
                
              }
            }
          }
        }
      }
      Z <- Z0[-(N+1)]
    }
    Z_mat[i,] <- Z
    
    # treatment coefficients ----
    if(regular_update){
      z_res0 <- Z - W_beta - Q_eta
    } else {
      z_res0 <- Z - W_beta
      rho2_nQ <- rho2 / (sigma2 + rho2*nQ)
      rQtX = rho2_nQ*QtX
      dXQrQtX = colSums(rho2_nQ*QtX^2)
    }
    
    for (k in 1:C) {
      treatment_coeffs[k] <- 0
      z_res <- z_res0 - Xsp%*%treatment_coeffs
      Sstat <- sumX2[k]
      Tstat <- sum(X[,k]*z_res)
      if(!regular_update){
        Sstat <- sumX2[k] - dXQrQtX[k]
        QZres <- as.vector(tQsp%*%z_res)
        Tstat <- sum(X[,k]*z_res) - sum(rQtX[,k]*QZres)
      }
      
      V1 <- 1 / (1/lambda0^2 + s2inv*Sstat)
      E1 <- V1 * (mu0/lambda0^2 + s2inv*Tstat)
      sqrtV1 = sqrt(V1)
      
      gamma1 <- rbinom(n = 1, size = 1, prob = prob_spike(pi0,mu0,lambda0,E1,sqrtV1))
      if (gamma1 != 1) {
        treatment_coeffs[k] <- rtruncnorm_plain(n=1, mean=E1, sd=sqrtV1, lb=0, ub=Inf)
        if(PX){
          Rstat <- s2inv*sum(z_res^2)
          if(!regular_update){ Rstat <- s2inv*sum(z_res^2) - s2inv*sum(rho2_nQ*QZres^2) }
          
          cC <- sqrtV1*s2inv*Tstat
          cD <- sqrtV1*mu0/lambda0^2
          cM <- exp(dnorm(cC+cD,log=T)-pnorm(cC+cD,log=T))
          cA <- Rstat - cC^2
          
          rPX_coef <- rgamma(n=1, shape=0.5*N, rate=0.5*(cA-cC*(cD+cM)))
          
          E1 <- V1 * (mu0/lambda0^2 + s2inv*Tstat*sqrt(rPX_coef))
          treatment_coeffs[k] <- rtruncnorm_plain(n=1, mean=E1, sd=sqrtV1, lb=0, ub=Inf)
        }
      }
    }
    treatment_mat[i, ] <- treatment_coeffs
    X_treat = as.vector(Xsp %*% treatment_coeffs)
    
    # beta _ group-specific intercepts ----
    if(regular_update){
      z_res <- Z - X_treat - Q_eta
      Wz_res <- as.vector( t(W) %*% z_res )
      Cov_beta <- solve(diag(G)/w0^2 + s2inv*WtW)
      loc_beta <- s2inv * Wz_res
    } else {
      z_res <- Z - X_treat
      Wz_res <- as.vector( t(W) %*% z_res )
      Qz_res <- as.vector( tQsp %*% z_res )
      rho2_nQ <- rho2 / (sigma2 + rho2*nQ)
      Cov_beta <- solve(diag(G)/w0^2 + s2inv*WtW - s2inv * crossprod(QtW,rho2_nQ*QtW))
      loc_beta <- s2inv * Wz_res - s2inv * colSums( (rho2_nQ*Qz_res) * QtW )
    }
    mu_beta <- Cov_beta%*%loc_beta
    beta <- as.matrix(mvrnorm(n=1, mu=mu_beta, Sigma=Cov_beta), ncol=1)
    beta_mat[i,] <- beta
    
    W_beta = as.vector(Wsp %*% beta)
    z_res <- Z - X_treat - W_beta
    Qz_res <- as.vector(tQsp %*% z_res)
    
    # rho2 _ RE variance (collapsed) ---- 
    if(!regular_update){
      rPX_rho2 <- 1
      if(PX){
        ss_proj <- s2inv*sum((1/nQ)*Qz_res^2)    
        ss_resid <- s2inv*sum(z_res^2) - ss_proj
        rho2_mode <- get_saddle_point(r_shape, r_rate, sigma2, R, mQ, ss_proj)
        ss_proj_corr <- ss_proj * sigma2 / (sigma2 + mQ*rho2_mode)
        rPX_rho2 <- rgamma(n=1, shape=0.5*N, rate=0.5*(ss_resid+ss_proj_corr))
      }
      HQz2 <- Hsp%*%(Qz_res^2)
      rho2 <- sample_rho2_REmarg(rho2_vals, sigma2, r_shape, r_rate, HQz2*rPX_rho2, nM, M)
      rho2_vec[i] <- rho2
    }
    
    # eta _ random effects (RE) ----
    sigma2_eta <- 1 / (1/rho2 + s2inv * nQ) 
    rPX_eta <-1
    if(PX) rPX_eta <- rgamma(n = 1, shape = 0.5*N,
                             rate = 0.5*(s2inv*sum(z_res^2)-s2inv^2*sum(sigma2_eta*Qz_res^2)))
    mu_eta <- sqrt(rPX_eta) * s2inv * sigma2_eta * Qz_res
    eta <- rnorm(R, mean = mu_eta, sd = sqrt(sigma2_eta))
    eta_mat[i,] <- eta
    Q_eta = as.vector(Qsp %*% eta)
  
    # rho2 _ RE variance (regular) ---- 
    if(regular_update){
      rho2 <- 1/rgamma(n = 1, shape = r_shape + R/2, rate = r_rate + (1/2) * sum(eta^2))
      rho2_vec[i] <- rho2
    } 
    
  }
  
  cat("\n\n")
  
  # Return output --------------------------------------------------------------
  
  # Copy column names
  colnames(treatment_mat) <- colnames(X)
  colnames(beta_mat) <- colnames(W)[-(G+1)]
  
  # Save mcmc settings
  mcmc_settings <- list(MC_seed = MC_seed, iter_MC = iter_MC, iter_Burn = iter_Burn,
                       marg_RE = marg_RE, PX = PX, rcpp = rcpp)
  hyper_params <- c(hyper_params, mcmc_settings)
  
  # Changing names to match paper syntax
  output <- list(
    alpha = treatment_mat,
    mu = beta_mat,
    eta = eta_mat,
    rho2 = rho2_vec,
    Z = Z_mat,
    hyper_params = hyper_params)
  
  return(output)
}

