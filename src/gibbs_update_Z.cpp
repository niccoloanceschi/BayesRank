#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' Vectorized sampling from univariate truncated normal distributions (inverse-CDF method)
//'
//' @param n Number of draws.
//' @param mean Vector of means of the underlying normal.
//' @param sd Vector of SDs of the underlying normal.
//' @param lb Lower truncation bound.
//' @param ub Upper truncation bound.
//'
//' @return Numeric vector of length n.
//'
//' @noRd
//' 
// [[Rcpp::export]]
arma::vec rtruncnorm_plain_cpp(int n, arma::vec mean, arma::vec sd, double lb, double ub) {
  vec u = runif(n);
  vec z_lb = (lb - mean) / sd;
  vec z_ub = (ub - mean) / sd;
  
  NumericVector p_lb = Rcpp::pnorm(as<NumericVector>(wrap(z_lb)));
  NumericVector p_ub = Rcpp::pnorm(as<NumericVector>(wrap(z_ub)));
  NumericVector p_unif = p_lb + as<NumericVector>(wrap(u)) * (p_ub - p_lb);
  
  NumericVector z_new = Rcpp::qnorm(p_unif);
  
  return mean + sd % as<vec>(z_new);
}

//' Vectorized sampling from multivariate truncated normal distributions
//'
//' Updates a block of correlated truncated-normal latent utilities,
//' by updating one coordinate at a time via inner Gibbs sampling sweeps.
//'
//' @param Z_m Matrix of current latent utility values for the group.
//' @param Z_mean_m Matrix of conditional means.
//' @param muC_m Conditioning adjustment term from the marginalized random effect.
//' @param rho2_m Random-effect variance term for this group size.
//' @param sd_m Conditional SD for the truncated normal draws.
//' @param lb_h,ub_h Truncation bounds (ordinal category cutoffs).
//' @param n_iter Number of inner Gibbs sweeps.
//'
//' @return Updated matrix of latent utility values, same shape as Z_m.
//'
//' @noRd
//' 
// [[Rcpp::export]]
arma::mat sample_MVTN_gibbs_cpp(arma::mat Z_m, arma::mat Z_mean_m, arma::vec muC_m, arma::vec rho2_m,
                                arma::vec sd_m, double lb_h, double ub_h, int n_iter) {
  
  int m = Z_m.n_rows;
  int q = Z_m.n_cols;
  
  mat Z_diff_m = Z_m - Z_mean_m;
  rowvec Z_sum_m = sum(Z_diff_m, 0);
  vec rho2_muC_m = rho2_m % muC_m;  // element-wise
  
  for(int tt = 0; tt < n_iter; tt++) {
    for(int s = 0; s < m; s++) {
      
      vec mu_m = trans(Z_mean_m.row(s)) + rho2_muC_m + 
        rho2_m % (trans(Z_sum_m) - trans(Z_diff_m.row(s)));  // element-wise 
      
      vec Z_new = rtruncnorm_plain_cpp(q, mu_m, sd_m, lb_h, ub_h);
      
      rowvec Z_new_diff = trans(Z_new) - Z_mean_m.row(s);
      Z_sum_m = Z_sum_m - Z_diff_m.row(s) + Z_new_diff;
      Z_diff_m.row(s) = Z_new_diff;
      Z_m.row(s) = trans(Z_new);
    }
  }
  
  return Z_m;
}