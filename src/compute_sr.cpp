#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Compute Sharpe Ratio
//'
//' Given a vector of portfolio weights, a vector of mean returns, and a covariance matrix,
//' this function computes the Sharpe ratio defined as
//' \deqn{\frac{w^T \mu}{\sqrt{w^T \Sigma w}}.}
//'
//' @param weights A numeric vector of portfolio weights.
//' @param mu A numeric vector of expected returns.
//' @param sigma A numeric covariance matrix.
//' @param do_checks Logical flag indicating whether to perform input checks (default = false).
//' @return A double representing the Sharpe ratio.
//' @export
// [[Rcpp::export]]
double compute_sr(const arma::vec& weights,
                 const arma::vec& mu,
                 const arma::mat& sigma,
                 const bool do_checks = false) {
 if (do_checks) {
   if (weights.n_elem == 0) {
     Rcpp::stop("weights must be non-empty");
   }
   if (mu.n_elem == 0) {
     Rcpp::stop("mu must be non-empty");
   }
   if (sigma.n_elem == 0) {
     Rcpp::stop("sigma must be provided");
   }
   if (weights.n_elem != mu.n_elem) {
     Rcpp::stop("weights and mu must be of the same length");
   }
   if (sigma.n_rows != sigma.n_cols) {
     Rcpp::stop("sigma must be a square matrix");
   }
   if (sigma.n_rows != weights.n_elem) {
     Rcpp::stop("The dimensions of sigma must match the length of weights");
   }
 }

 return arma::dot(weights, mu) / std::sqrt(arma::as_scalar(weights.t() * sigma * weights));
}

