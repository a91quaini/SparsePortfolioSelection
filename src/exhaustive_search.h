// exhaustive_search.h
// Exhaustive (and sampled) search for best k-subset MVE Sharpe ratio.
#ifndef EXHAUSTIVE_SEARCH_H
#define EXHAUSTIVE_SEARCH_H

#include <RcppArmadillo.h>

// Search for the k-asset subset maximizing MVE Sharpe ratio (variance-based).
// Returns an R list: selection (uvec), weights (vec), sr (double), status ("EXHAUSTIVE"|"SAMPLED").
//
// Notes:
// - Indices are 0-based (consistent with utils/sharpe_ratio).
// - ridge_epsilon is applied once to sigma via stabilize_sigma_cpp; default is 0 (no regularization).
Rcpp::List mve_exhaustive_search_cpp(const arma::vec& mu,
                                     const arma::mat& sigma,
                                     unsigned int k,
                                     double ridge_epsilon = 0.0,
                                     bool enumerate_all = true,
                                     unsigned int max_samples = 0,
                                     bool dedup_samples = false,
                                     bool compute_weights = false,
                                     bool normalize_weights = false,
                                     int normalization_type = 1,
                                     bool do_checks = false);

#endif // EXHAUSTIVE_SEARCH_H
