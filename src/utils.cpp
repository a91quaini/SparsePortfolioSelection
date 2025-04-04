#include "utils.h"
#include <algorithm>
#include <random>
#include <cmath>
#include <limits>

// Compute Portfolio Sharpe Ratio: (w^T mu) / sqrt(w^T sigma w)
double compute_sr(const arma::vec& weights,
                  const arma::vec& mu,
                  const arma::mat& sigma,
                  const bool do_checks) {

  // Optional input checks.
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

  // Compute the Sharpe ratio: (w^T mu) / sqrt(w^T sigma w)
  return arma::dot(weights, mu) / std::sqrt(arma::as_scalar(weights.t() * sigma * weights));
}


// Compute nCk (combinations count)
unsigned long long nCk(const unsigned int n, const unsigned int k) {
  if (k > n) return 0;
  if (k == 0 || k == n) return 1;
  unsigned long long res = 1;
  for (unsigned int i = 1; i <= k; i++) {
    res = res * (n - i + 1) / i;
  }
  return res;
}

// Recursively generate all combinations of k indices from {0,1,...,n-1}
void generate_combinations(const unsigned int n, const unsigned int k, const unsigned int offset,
                           arma::uvec &current, std::vector<arma::uvec> &all) {
  if (current.n_elem == k) {
    all.push_back(current);
    return;
  }
  for (unsigned int i = offset; i < n; i++) {
    current.insert_rows(current.n_elem, 1);
    current(current.n_elem - 1) = i;
    generate_combinations(n, k, i + 1, current, all);
    current.shed_row(current.n_elem - 1);
  }
}

// Generate a random combination of k distinct indices from 0 to n-1
arma::uvec random_combination(const unsigned int n, const unsigned int k) {
  std::vector<unsigned int> indices(n);
  for (unsigned int i = 0; i < n; i++) {
    indices[i] = i;
  }
  std::random_device rd;
  std::mt19937 g(rd());
  std::shuffle(indices.begin(), indices.end(), g);
  std::sort(indices.begin(), indices.begin() + k);
  arma::uvec comb(k);
  for (unsigned int i = 0; i < k; i++) {
    comb(i) = indices[i];
  }
  return comb;
}
