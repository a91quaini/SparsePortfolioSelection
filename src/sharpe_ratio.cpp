// sharpe_ratio.cpp

#include "sharpe_ratio.h"
#include <cmath>
#include <limits>

using std::numeric_limits;

inline void check_mu_sigma(const arma::vec& mu, const arma::mat& sigma) {
  if (sigma.n_rows != sigma.n_cols) Rcpp::stop("sigma must be square");
  if (mu.n_elem != sigma.n_rows)   Rcpp::stop("mu and sigma not conformable");
  if (!mu.is_finite() || !sigma.is_finite()) Rcpp::stop("Non-finite inputs.");
}

inline void check_selection(const arma::uvec& sel, arma::uword n) {
  if (sel.n_elem == 0) return;
  if (sel.max() >= n) Rcpp::stop("selection out of bounds 0..n-1.");
}

// -----------------------------------------------------------------------------
// Public API
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
double compute_sr_cpp(const arma::vec& weights,
                      const arma::vec& mu,
                      const arma::mat& sigma,
                      const arma::uvec& selection,
                      double ridge_epsilon,
                      bool do_checks) {

  const arma::uword n = mu.n_elem;
  if (do_checks) {
    check_mu_sigma(mu, sigma);
    if (weights.n_elem != n) Rcpp::stop("weights and mu must have the same length.");
    check_selection(selection, n);
    if (!weights.is_finite()) Rcpp::stop("weights must be finite.");
    if (!std::isfinite(ridge_epsilon)) Rcpp::stop("epsilon must be finite.");
  }

  arma::mat sigma_ = stabilize_sigma_cpp(sigma, ridge_epsilon);

  double num = 0.0;
  double variance = 0.0;

  if (selection.n_elem == 0 || selection.n_elem == n) {
    num = arma::dot(weights, mu);
    variance = arma::dot(weights, sigma_ * weights);
  } else {
    const arma::vec w  = weights.elem(selection);
    const arma::vec mu_s = mu.elem(selection);
    const arma::mat sigma_s = sigma_.submat(selection, selection);
    num = arma::dot(w, mu_s);
    variance = arma::dot(w, sigma_s * w);
  }

  if (!std::isfinite(variance) || !(variance > 0.0)) {
    return numeric_limits<double>::quiet_NaN();
  }

  const double sr = num / std::sqrt(variance);
  return std::isfinite(sr) ? sr : numeric_limits<double>::quiet_NaN();
}

// [[Rcpp::export]]
double compute_mve_sr_cpp(const arma::vec& mu,
                          const arma::mat& sigma,
                          const arma::uvec& selection,
                          double ridge_epsilon,
                          bool do_checks) {

  const arma::uword n = mu.n_elem;
  if (do_checks) {
    if (n == 0) Rcpp::stop("mu must be non-empty.");
    check_mu_sigma(mu, sigma);
    check_selection(selection, n);
    if (!std::isfinite(ridge_epsilon)) Rcpp::stop("epsilon must be finite.");
  }

  arma::mat sigma_ = stabilize_sigma_cpp(sigma, ridge_epsilon);

  arma::vec mu_s;
  arma::mat sigma_s;
  if (selection.n_elem == 0 || selection.n_elem == n) {
    mu_s = mu;
    sigma_s = sigma_;
  } else {
    mu_s = mu.elem(selection);
    sigma_s = sigma_.submat(selection, selection);
  }

  const arma::vec x = solve_sympd(sigma_s, mu_s);
  const double val = arma::dot(mu_s, x);
  const double sr_sq = std::max(val, 0.0);
  const double sr = std::sqrt(sr_sq);
  return std::isfinite(sr) ? sr : numeric_limits<double>::quiet_NaN();
}

// [[Rcpp::export]]
arma::vec compute_mve_weights_cpp(const arma::vec& mu,
                                  const arma::mat& sigma,
                                  const arma::uvec& selection,
                                  double ridge_epsilon,
                                  bool normalize_weights,
                                  int normalization_type,
                                  bool do_checks) {

  const arma::uword n = mu.n_elem;
  if (do_checks) {
    if (n == 0) Rcpp::stop("mu must be non-empty.");
    check_mu_sigma(mu, sigma);
    check_selection(selection, n);
    if (!std::isfinite(ridge_epsilon)) Rcpp::stop("epsilon must be finite.");
  }

  arma::mat sigma_ = stabilize_sigma_cpp(sigma, ridge_epsilon);

  arma::vec w(n, arma::fill::zeros);
  if (selection.n_elem == 0 || selection.n_elem == n) {
    w = solve_sympd(sigma_, mu);
  } else {
    const arma::vec mu_s = mu.elem(selection);
    const arma::mat sigma_s = sigma_.submat(selection, selection);
    arma::vec ws = solve_sympd(sigma_s, mu_s);
    w.elem(selection) = ws;
  }

  if (normalize_weights) {
    w = normalize_weights_cpp(w, 1e-8, normalization_type);
  }

  return w;
}
