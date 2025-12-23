#include "utils.h"
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

// ---- 1) Linear solves ----------------------------------------------------

template<typename T>
inline T solve_sympd(const arma::mat& A, const T& b) {
  T x;

  // 1) Strict SPD attempt: likely_sympd + no_approx
  bool ok = arma::solve(
    x, A, b,
    arma::solve_opts::likely_sympd + arma::solve_opts::no_approx
  );
  if (ok) return x;

  // 2) Fallback: allow approximation (SVD / pinv-like)
  ok = arma::solve(
    x, A, b,
    arma::solve_opts::force_approx
  );
  if (ok) return x;

  // 3) Final fallback: explicit pseudoinverse (guarantees a return)
  return arma::pinv(A) * b;
}
// Explicit instantiations (optional but can speed up linking)
template arma::vec solve_sympd<arma::vec>(const arma::mat&, const arma::vec&);
template arma::mat solve_sympd<arma::mat>(const arma::mat&, const arma::mat&);

// Numerical knobs / utilities ---------------------------------------
// [[Rcpp::export]]
arma::mat stabilize_sigma_cpp(const arma::mat& sigma,
                              double ridge_epsilon) {

  if (ridge_epsilon <= 1e-14) {
    return sigma;
  }

  arma::mat sigma_ = sigma;
  for (arma::uword i = 0; i < sigma.n_rows; ++i) {
    sigma_(i, i) += ridge_epsilon;
  }

  return sigma_;
}

// [[Rcpp::export]]
arma::vec normalize_weights_cpp(const arma::vec& w,
                                double epsilon,
                                int type) {
  // light guard
  if (!std::isfinite(epsilon) || epsilon <= 0.0) epsilon = 1e-8;

  arma::vec out = w;

  // support inferred from w
  const arma::uvec sel = arma::find(arma::abs(out) > epsilon);
  if (sel.n_elem == 0) {
    return out;
  }

  if (type == 0) {  // sum-to-one
    const double sum_w = arma::accu(out(sel));
    if (!std::isfinite(sum_w) || std::abs(sum_w) <= 0.0) {
      return out;
    }
    out(sel) /= sum_w;

  } else {  // default: L1 norm = 1
    const double l1norm = arma::accu(arma::abs(out(sel)));
    if (!std::isfinite(l1norm)) {
      return out;
    }
    out(sel) /= l1norm;

  }
  return out;
}

// [[Rcpp::export]]
arma::mat safe_chol_cpp(const arma::mat& Q,
                        double ridge_epsilon) {

  // Symmetrize defensively (Cholesky expects symmetric PD)
  arma::mat Qs = arma::symmatu(Q); // equivalent to 0.5*(Q + Q.t())

  // Optional diagonal stabilization via your helper
  Qs = stabilize_sigma_cpp(Qs, ridge_epsilon);

  arma::mat U;
  if (!arma::chol(U, Qs, "upper")) {
    Rcpp::stop("safe_chol_cpp: Cholesky failed.");
  }
  return U; // Qs = U.t() * U
}

// [[Rcpp::export]]
Rcpp::List design_from_moments_cpp(const arma::vec& mu,
                                   const arma::mat& sigma,
                                   double n_obs,
                                   double ridge_epsilon) {
  // Markowitz (covariance-only) quadratic form:
  // Q = T * Sigma  (no mu*mu' term)
  const arma::mat Q = n_obs * sigma;

  // Cholesky of symmetrized (and optionally epsilon-stabilized) Q:
  // returns upper-triangular U such that Qs = U.t() * U
  const arma::mat U = safe_chol_cpp(Q, ridge_epsilon);

  // Design:
  // X'X = U' U = Q
  // Choose y so that X'y = U' y = T mu  => y solves (U') y = T mu
  const arma::mat X = U;
  const arma::vec y = arma::solve(arma::trimatl(U.t()), n_obs * mu,
                                  arma::solve_opts::fast);

  return Rcpp::List::create(
    Rcpp::Named("X") = X,
    Rcpp::Named("y") = y
  );
}

// ---- Internal helpers (anonymous namespace) -----------------------------
namespace {
// Recursive helper that pushes all combinations into 'out'. Indices 0-based.
void gen_all_rec(unsigned int n, unsigned int k, unsigned int offset,
                 arma::uvec& cur, std::vector<arma::uvec>& out) {
  if (cur.n_elem == k) { out.emplace_back(cur); return; }
  unsigned int need = k - cur.n_elem;
  for (unsigned int i = offset; i <= n - need; ++i) {
    cur.insert_rows(cur.n_elem, 1);
    cur(cur.n_elem - 1) = i;
    gen_all_rec(n, k, i + 1, cur, out);
    cur.shed_row(cur.n_elem - 1);
  }
}

// Recursive helper that streams each combination to a callback.
void gen_cb_rec(unsigned int n, unsigned int k, unsigned int offset,
                arma::uvec& cur, const CombCallback& cb) {
  if (cur.n_elem == k) { cb(cur); return; }
  unsigned int need = k - cur.n_elem;
  for (unsigned int i = offset; i <= n - need; ++i) {
    cur.insert_rows(cur.n_elem, 1);
    cur(cur.n_elem - 1) = i;
    gen_cb_rec(n, k, i + 1, cur, cb);
    cur.shed_row(cur.n_elem - 1);
  }
}
} // anonymous namespace

// ---- Public API ----------------------------------------------------------

void generate_combinations(unsigned int n, unsigned int k,
                           std::vector<arma::uvec>& out) {
  out.clear();
  if (k == 0) { out.emplace_back(arma::uvec()); return; }
  if (k > n) return;

  arma::uvec cur;
  gen_all_rec(n, k, 0u, cur, out);
}

void for_each_combination(unsigned int n, unsigned int k,
                          const CombCallback& cb) {
  if (k == 0) { arma::uvec empty; cb(empty); return; }
  if (k > n) return;

  arma::uvec cur;
  gen_cb_rec(n, k, 0u, cur, cb);
}
