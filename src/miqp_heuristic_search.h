// miqp_heuristic_search.h
#ifndef MIQP_HEURISTIC_SEARCH_H
#define MIQP_HEURISTIC_SEARCH_H

#include <RcppArmadillo.h>

Rcpp::List mve_miqp_search_cpp(const arma::vec& mu,
                               const arma::mat& sigma,
                               unsigned int k,
                               const arma::vec& fmin,
                               const arma::vec& fmax,
                               Rcpp::Nullable<arma::vec> x_start,
                               Rcpp::Nullable<arma::uvec> v_start,
                               int m = 1,
                               double gamma = 1.0,
                               bool exactly_k = true,
                               unsigned int expand_rounds = 5,
                               double expand_factor = 3.0,
                               double expand_tol = 1e-2,
                               double mipgap = 1e-4,
                               double time_limit = 100,
                               int threads = 0,
                               double ridge_epsilon = 0.0,
                               bool normalize_weights = false,
                               bool use_refit = false,
                               bool verbose = false,
                               bool do_checks = false);

#endif // MIQP_HEURISTIC_SEARCH_H
