// Author: Alberto Quaini
#ifndef UTILS_H
#define UTILS_H

#include <RcppArmadillo.h>
#include <functional>  // <-- needed for std::function

// Set template T, representing both arma::mat or arma::vec in subsequent functions.
template<typename T>
T solve_sympd(const arma::mat& A, const T& b);

// Numerical knobs / utilities ---------------------------------------

// Symmetrize a covariance matrix and optionally add a ridge term proportional
// to its average diagonal. Optionally convert to second moments (Sigma + mu mu').
arma::mat stabilize_sigma_cpp(const arma::mat& sigma,
                              double ridge_epsilon = 1e-8);

// Normalize a weight vector by its L1 norm (gross exposure):
// identify the support sel = {i : |w_i| > tol}, keep w_out[-sel]=0,
// and scale w_out[sel] so that sum_i |w_out,i| = 1.
arma::vec normalize_weights_cpp(const arma::vec& w,
                                double epsilon = 1e-8,
                                int type = 1);  // type: 0 -> sum-to-1, 1 -> L1=1

// Robust Cholesky with escalating diagonal bump; returns upper-triangular U.
arma::mat safe_chol_cpp(const arma::mat& Q,
                        double ridge_epsilon = 1e-8);

// Design builder from moments (used by LASSO search).
Rcpp::List design_from_moments_cpp(const arma::vec& mu,
                                   const arma::mat& sigma,
                                   double n_obs,
                                   double ridge_epsilon = 1e-8);

// Combinations API ------------------------------------------------

// Streaming callback type: receives a size-k index vector (0..n-1)
using CombCallback = std::function<void(const arma::uvec&)>;

// Materialize all k-combinations of {0,...,n-1} into 'out' (cleared first).
// For large nCk this can be very memory-heavy; prefer for_each_combination.
void generate_combinations(unsigned int n, unsigned int k,
                           std::vector<arma::uvec>& out);

// Streaming version: invokes cb(cur) for each k-combination without storing all.
void for_each_combination(unsigned int n, unsigned int k,
                          const CombCallback& cb);

#endif // UTILS_H
