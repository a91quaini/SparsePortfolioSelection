// Author: Alberto Quaini
#ifndef UTILS_H
#define UTILS_H

#include <RcppArmadillo.h>
#include <string>
#include <vector>
#include <functional>  // <-- needed for std::function

// Set template T, representing both arma::mat or arma::vec in subsequent functions.
template<typename T>
T solve_sympd(const arma::mat& A, const T& b);

// ---- 0) Numerical knobs / utilities ---------------------------------------

// Small ridge used to stabilize covariance matrices.
constexpr double EPS_RIDGE = 1e-6;

// Symmetrize a covariance matrix and optionally add a ridge term proportional
// to its average diagonal.
arma::mat prep_covariance(const arma::mat& sigma,
                          double epsilon,
                          bool stabilize);

// Normalize a weight vector using either "absolute" (scale by |sum|) or
// "relative" (scale by max(|sum|, tol * L1)) strategies.
arma::vec normalize_weights(const arma::vec& w,
                            const std::string& mode = "relative",
                            double tol = 1e-6,
                            bool do_checks = false);

// Rcpp exports
arma::mat prep_covariance_cpp(const arma::mat& sigma,
                              double epsilon,
                              bool stabilize);

arma::vec normalize_weights_cpp(const arma::vec& w,
                                const std::string& mode,
                                double tol,
                                bool do_checks);

double eps_ridge_cpp();

// Robust Cholesky with escalating diagonal bump; returns upper-triangular U.
arma::mat safe_chol_cpp(const arma::mat& Q,
                        double base_bump = 1e-10,
                        unsigned int max_tries = 8);

// Design builder from moments (used by LASSO path).
Rcpp::List design_from_moments_cpp(const arma::vec& mu,
                                   const arma::mat& sigma,
                                   double n_obs,
                                   double epsilon,
                                   bool stabilize_sigma);

// ---- 2) Combinations API ------------------------------------------------

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
