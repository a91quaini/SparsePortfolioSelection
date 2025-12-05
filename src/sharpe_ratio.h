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
                      double epsilon,
                      bool stabilize_sigma,
                      bool do_checks);

// Compute SR of the mean–variance efficient portfolio (sqrt(mu' Sigma^{-1} mu)).
double compute_mve_sr_cpp(const arma::vec& mu,
                          const arma::mat& sigma,
                          const arma::uvec& selection,
                          double epsilon,
                          bool stabilize_sigma,
                          bool do_checks);

// Compute MVE portfolio weights w = Sigma^{-1} mu with optional normalization.
arma::vec compute_mve_weights_cpp(const arma::vec& mu,
                                  const arma::mat& sigma,
                                  const arma::uvec& selection,
                                  bool normalize_w,
                                  double epsilon,
                                  bool stabilize_sigma,
                                  bool do_checks);

#endif // SHARPE_RATIO_H
