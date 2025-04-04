#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Compute Optimal Portfolio Sharpe Ratio
//'
//' Given a mean vector \eqn{\mu}, a covariance matrix \eqn{\Sigma},
//' and an asset selection vector, this function computes the optimal Sharpe ratio defined as
//' \deqn{\sqrt{\mu^T \Sigma^{-1}\mu}}
//' over the selected assets.
//' If the provided asset selection vector has length less than N,
//' the function uses the subset of assets and computes the ratio using those assets.
//'
//' @param mu Mean vector.
//' @param sigma Covariance matrix.
//' @param selection Unsigned integer vector with asset indices.
//' @param do_checks Logical flag indicating whether to perform input checks (default = FALSE).
//' @return A scalar value corresponding to \eqn{\sqrt{\mu^T \Sigma^{-1}\mu}}.
//' @export
// [[Rcpp::export]]
double compute_mve_sr(const arma::vec& mu,
                const arma::mat& sigma,
                const arma::uvec& selection,
                const bool do_checks = false) {

 // Optional input checks.
 if (do_checks) {
   if (mu.n_elem == 0 || sigma.n_elem == 0) {
     Rcpp::stop("Mean vector mu and covariance matrix sigma must be supplied");
   }
   if (sigma.n_rows != sigma.n_cols) {
     Rcpp::stop("Covariance matrix sigma must be square");
   }
   if (mu.n_elem != sigma.n_rows) {
     Rcpp::stop("Mean vector mu must have a length equal to the dimensions of sigma");
   }
   if (selection.n_elem > 0 && arma::max(selection) >= mu.n_elem) {
     Rcpp::stop("Asset selection index out of bounds");
   }
 }

 // If selection is not supplied or is full, use full mu and sigma.
 if (selection.n_elem == 0 || selection.n_elem == mu.n_elem) {
   return std::sqrt( arma::dot(mu, arma::solve(sigma, mu)) );
 }

 // Otherwise, subset the inputs according to the asset selection.
 const arma::vec mu_sel = mu.elem(selection);
 const arma::mat sigma_sel = sigma.submat(selection, selection);

 return std::sqrt( arma::dot(mu_sel, arma::solve(sigma_sel, mu_sel)) );
}

