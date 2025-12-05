// sharpe_ratio.cpp
// Port of Julia's SharpeRatio utilities to RcppArmadillo.

#include "sharpe_ratio.h"
#include <cmath>
#include <limits>

using std::numeric_limits;

namespace {
inline double sr_eps() { return std::sqrt(arma::datum::eps); }

inline void check_mu_sigma(const arma::vec& mu, const arma::mat& sigma) {
  if (sigma.n_rows != sigma.n_cols) Rcpp::stop("sigma must be square");
  if (mu.n_elem != sigma.n_rows)   Rcpp::stop("mu and sigma not conformable");
  if (!mu.is_finite() || !sigma.is_finite()) Rcpp::stop("Non-finite inputs.");
}

inline void check_selection(const arma::uvec& sel, arma::uword n) {
  if (sel.n_elem == 0) return;
  if (sel.max() >= n) Rcpp::stop("selection out of bounds 0..n-1.");
}

// Solve Sigma * x = b for symmetric Sigma with fallbacks (chol -> solve -> pinv -> ridge+chol).
arma::vec sym_solve(const arma::mat& Sigma, const arma::vec& b) {
  arma::mat A = Sigma; // local copy for potential ridge
  arma::mat L;

  // 1) Try Cholesky
  if (arma::chol(L, A, "lower")) {
    arma::vec y = arma::solve(arma::trimatl(L), b);
    return arma::solve(arma::trimatu(L.t()), y);
  }

  // 2) Try SPD solve (likely_sympd + no_approx)
  arma::vec x;
  bool ok = arma::solve(
    x, A, b,
    arma::solve_opts::likely_sympd + arma::solve_opts::no_approx
  );
  if (ok) return x;

  // 3) Pseudoinverse
  x = arma::pinv(A) * b;
  if (x.is_finite()) return x;

  // 4) Add tiny ridge, then try again with Cholesky
  const double delta = std::max(arma::datum::eps, 1e-10 * arma::mean(A.diag()));
  A.diag() += delta;
  if (arma::chol(L, A, "lower")) {
    arma::vec y = arma::solve(arma::trimatl(L), b);
    return arma::solve(arma::trimatu(L.t()), y);
  }

  // Last resort: fall back to solve_sympd from utils (allows approximation).
  return solve_sympd<arma::vec>(A, b);
}
} // namespace

// -----------------------------------------------------------------------------
// Public API
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
double compute_sr_cpp(const arma::vec& weights,
                      const arma::vec& mu,
                      const arma::mat& sigma,
                      const arma::uvec& selection,
                      double epsilon,
                      bool stabilize_sigma,
                      bool do_checks) {

  const arma::uword n = mu.n_elem;
  if (do_checks) {
    check_mu_sigma(mu, sigma);
    if (weights.n_elem != n) Rcpp::stop("weights and mu must have the same length.");
    check_selection(selection, n);
    if (!weights.is_finite()) Rcpp::stop("weights must be finite.");
    if (!std::isfinite(epsilon)) Rcpp::stop("epsilon must be finite.");
  }

  const arma::mat Sigma_eff = prep_covariance(sigma, epsilon, stabilize_sigma);

  double num = 0.0;
  double variance = 0.0;

  if (selection.n_elem == 0 || selection.n_elem == n) {
    num = arma::dot(weights, mu);
    variance = arma::dot(weights, Sigma_eff * weights);
  } else {
    arma::uvec sel = selection;
    arma::vec w  = weights.elem(sel);
    arma::vec mu_s = mu.elem(sel);
    arma::mat Sigma_s = Sigma_eff.submat(sel, sel);
    num = arma::dot(w, mu_s);
    variance = arma::dot(w, Sigma_s * w);
  }

  if (!std::isfinite(variance) || !(variance > 0.0)) {
    return numeric_limits<double>::quiet_NaN();
  }

  const double den = std::sqrt(variance);
  if (!(den > sr_eps()) || !std::isfinite(den)) {
    return numeric_limits<double>::quiet_NaN();
  }

  const double sr = num / den;
  return std::isfinite(sr) ? sr : numeric_limits<double>::quiet_NaN();
}

// [[Rcpp::export]]
double compute_mve_sr_cpp(const arma::vec& mu,
                          const arma::mat& sigma,
                          const arma::uvec& selection,
                          double epsilon,
                          bool stabilize_sigma,
                          bool do_checks) {

  const arma::uword n = mu.n_elem;
  if (do_checks) {
    if (n == 0) Rcpp::stop("mu must be non-empty.");
    check_mu_sigma(mu, sigma);
    check_selection(selection, n);
    if (!std::isfinite(epsilon)) Rcpp::stop("epsilon must be finite.");
  }

  const arma::mat Sigma_eff = prep_covariance(sigma, epsilon, stabilize_sigma);

  arma::vec mu_s;
  arma::mat Sigma_s;
  if (selection.n_elem == 0 || selection.n_elem == n) {
    mu_s = mu;
    Sigma_s = Sigma_eff;
  } else {
    arma::uvec sel = selection;
    mu_s = mu.elem(sel);
    Sigma_s = Sigma_eff.submat(sel, sel);
  }

  const arma::vec x = sym_solve(Sigma_s, mu_s);
  const double val = arma::dot(mu_s, x);
  const double sr_sq = std::max(val, 0.0);
  const double sr = std::sqrt(sr_sq);
  return std::isfinite(sr) ? sr : numeric_limits<double>::quiet_NaN();
}

// [[Rcpp::export]]
arma::vec compute_mve_weights_cpp(const arma::vec& mu,
                                  const arma::mat& sigma,
                                  const arma::uvec& selection,
                                  bool normalize_w,
                                  double epsilon,
                                  bool stabilize_sigma,
                                  bool do_checks) {

  const arma::uword n = mu.n_elem;
  if (do_checks) {
    if (n == 0) Rcpp::stop("mu must be non-empty.");
    check_mu_sigma(mu, sigma);
    check_selection(selection, n);
    if (!std::isfinite(epsilon)) Rcpp::stop("epsilon must be finite.");
  }

  const arma::mat Sigma_eff = prep_covariance(sigma, epsilon, stabilize_sigma);

  arma::vec w(n, arma::fill::zeros);
  if (selection.n_elem == 0 || selection.n_elem == n) {
    w = sym_solve(Sigma_eff, mu);
  } else {
    arma::uvec sel = selection;
    arma::vec mu_s = mu.elem(sel);
    arma::mat Sigma_s = Sigma_eff.submat(sel, sel);
    arma::vec ws = sym_solve(Sigma_s, mu_s);
    w.elem(sel) = ws;
  }

  if (normalize_w) {
    const double norm1 = arma::norm(w, 1);
    const double s = arma::accu(w);
    if (!(norm1 > 1e-12) && !(std::abs(s) > 1e-12)) {
      return arma::vec(w.n_elem, arma::fill::zeros);
    }
    w = normalize_weights(w, "relative", 1e-6, false);
  }

  return w;
}
