// miqp_heuristic_search.h
// Mixed-integer quadratic heuristic for sparse MVE selection (Gurobi backend).
#ifndef MIQP_HEURISTIC_SEARCH_H
#define MIQP_HEURISTIC_SEARCH_H

#include <RcppArmadillo.h>

Rcpp::List mve_miqp_search_cpp(const arma::vec& mu,
                               const arma::mat& sigma,
                               unsigned int k,
                               bool exactly_k,
                               int m,
                               double gamma,
                               const arma::vec& fmin,
                               const arma::vec& fmax,
                               unsigned int expand_rounds,
                               double expand_factor,
                               double expand_tol,
                               double mipgap,
                               double time_limit,
                               int threads,
                               Rcpp::Nullable<arma::vec> x_start,
                               Rcpp::Nullable<arma::uvec> v_start,
                               bool compute_weights,
                               bool normalize_weights,
                               bool use_refit,
                               bool verbose,
                               double epsilon,
                               bool stabilize_sigma,
                               bool do_checks);

#endif // MIQP_HEURISTIC_SEARCH_H
