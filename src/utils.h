#ifndef UTILS_H
#define UTILS_H

#include <RcppArmadillo.h>
#include <vector>

//' Compute Portfolio Sharpe Ratio
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
                  const bool do_checks = false);

// Helper: Compute nCk (number of combinations of n elements in k places)
unsigned long long nCk(const unsigned int n, const unsigned int k);

// Helper: Recursively generate all combinations of k indices from {0,1,...,n-1}
void generate_combinations(const unsigned int n, const unsigned int k, const unsigned int offset,
                           arma::uvec &current, std::vector<arma::uvec> &all);

// Helper: Generate a random combination of k distinct indices from 0 to n-1
arma::uvec random_combination(const unsigned int n, const unsigned int k);

#endif // UTILS_H
