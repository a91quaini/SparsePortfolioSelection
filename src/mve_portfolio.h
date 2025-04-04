#ifndef MVE_PORTFOLIO_H
#define MVE_PORTFOLIO_H

#include <RcppArmadillo.h>

//' Compute Mean-Variance Efficient Portfolio Sharpe Ratio
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
                      const bool do_checks = false);

//' Compute Mean-Variance Efficient (MVE) Portfolio Weights
//'
//' Given a risk-aversion parameter \eqn{\gamma}, a first moment vector, a second moment matrix,
//' and a selection index vector, this function computes the mean variance efficient portfolio weights
//' \deqn{w = \frac{1}{\gamma}\, \text{second\_moment}^{-1}\, \text{first\_moment}},
//' over the selected assets.
//' If the provided asset selection vector has length less than N,
//' the function returns an N-length weight vector with zeros for the unselected assets.
//'
//' @param first_moment First moment vector.
//' @param second_moment Second moment matrix.
//' @param selection Unsigned integer vector with asset indices.
//' @param gamma Risk aversion parameter. Default is 1.
//' @param do_checks Logical flag indicating whether to perform input checks (default = FALSE).
//' @return An N-length vector of mean variance efficient weights.
//' @export
// [[Rcpp::export]]
arma::vec compute_mve_weights(const arma::vec& first_moment,
                              const arma::mat& second_moment,
                              const arma::uvec& selection,
                              const double gamma = 1.0,
                              const bool do_checks = false);

//' Compute Mean-Variance Efficient Sharpe Ratio with Cardinality Constraint
//'
//' This function takes as inputs the mean vector \code{mu}, the covariance matrix \code{sigma},
//' the maximum active portfolio cardinality \code{max_card},
//' and the fraction of combinations to evaluate \code{greedy_perc}.
//' With these inputs, it searches over all combinations of assets with cardinality from 1 up to \code{max_card}
//' and computes the square Sharpe ratio defined as
//' \eqn{\mu^T \, \sigma^{-1}\, \mu}.
//' It returns the highest square Sharpe ratio found along with the associated asset selection.
//' If \code{greedy_perc} is less than 1, then for each cardinality the search is performed over a random
//' subset (a fraction equal to \code{greedy_perc}) of all possible combinations.
//'
//' @param mu Mean vector.
//' @param sigma Coveriance matrix.
//' @param max_card Maximum cardinality to consider (from 1 up to the number of assets).
//' @param greedy_perc If less than 1, the fraction of combinations to evaluate for each cardinality;
//'                    if 1 or greater, all combinations are evaluated;
//'                    if less than 0, no combinations are evaluated.
//' @param do_checks Logical flag indicating whether to perform input checks (default = FALSE).
//' @return A list with \code{sqsr} (the optimal square Sharpe ratio) and \code{selection} (the asset indices of the optimal selection).
//' @export
// [[Rcpp::export]]
Rcpp::List compute_sparse_mve_sr(const arma::vec& mu,
                                 const arma::mat& sigma,
                                 unsigned int max_card = 1,
                                 const double greedy_perc = 1.0,
                                 const bool do_checks = false);

#endif // MVE_PORTFOLIO_H
