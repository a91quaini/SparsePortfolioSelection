// exhaustive_search.h
// Exhaustive (and sampled) search for best k-subset MVE Sharpe ratio.
#ifndef EXHAUSTIVE_SEARCH_H
#define EXHAUSTIVE_SEARCH_H

#include <RcppArmadillo.h>

// Search for the k-asset subset maximizing MVE Sharpe ratio.
// Returns an R list: selection (uvec), weights (vec), sr (double), status ("EXHAUSTIVE"|"SAMPLED").
Rcpp::List mve_exhaustive_search_cpp(const arma::vec& mu,
                                     const arma::mat& sigma,
                                     unsigned int k,
                                     double epsilon,
                                     bool stabilize_sigma,
                                     bool do_checks,
                                     bool enumerate_all,
                                     unsigned int max_samples,
                                     bool dedup_samples,
                                     bool compute_weights);

#endif // EXHAUSTIVE_SEARCH_H
