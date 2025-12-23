// sharpe_ratio.h
// Utilities for Sharpe ratio and MVE weights mirroring Julia's SharpeRatio module.
#ifndef SHARPE_RATIO_H
#define SHARPE_RATIO_H

#include <RcppArmadillo.h>
#include "utils.h"

// Compute Sharpe ratio for given weights, allowing optional asset selection and
// covariance stabilization.
double compute_sr_cpp(const arma::vec& weights,
                      const arma::vec& mu,
                      const arma::mat& sigma,
                      const arma::uvec& selection,
                      double ridge_epsilon = 1e-8,
                      bool do_checks = false);

// Compute SR of the mean–variance efficient portfolio (sqrt(mu' Sigma^{-1} mu)).
double compute_mve_sr_cpp(const arma::vec& mu,
                          const arma::mat& sigma,
                          const arma::uvec& selection,
                          double ridge_epsilon = 1e-8,
                          bool do_checks = false);

// Compute MVE portfolio weights w = Sigma^{-1} mu with optional normalization.
arma::vec compute_mve_weights_cpp(const arma::vec& mu,
                                  const arma::mat& sigma,
                                  const arma::uvec& selection,
                                  double ridge_epsilon = 1e-8,
                                  bool normalize_weights = false,
                                  int normalization_type = 1,
                                  bool do_checks = false);

#endif // SHARPE_RATIO_H
