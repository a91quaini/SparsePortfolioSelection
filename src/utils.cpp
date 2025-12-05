#include "utils.h"
#include <algorithm>
#include <cmath>

// ---- 0) Numerical knobs / utilities ---------------------------------------

arma::mat prep_covariance(const arma::mat& sigma,
                          double epsilon,
                          bool stabilize) {
  if (sigma.n_rows != sigma.n_cols) {
    Rcpp::stop("sigma must be square");
  }

  arma::mat A = arma::mat(sigma); // copy to ensure mutability/double
  const arma::uword n = A.n_rows;

  // Symmetrize in-place: A := (A + A.t()) / 2
  for (arma::uword j = 0; j < n; ++j) {
    A(j, j) = 0.5 * (A(j, j) + A(j, j));
    for (arma::uword i = j + 1; i < n; ++i) {
      const double s = 0.5 * (A(i, j) + A(j, i));
      A(i, j) = s;
      A(j, i) = s;
    }
  }

  if (stabilize && epsilon > 0.0) {
    const double ss = epsilon * (arma::trace(A) / static_cast<double>(n));
    for (arma::uword i = 0; i < n; ++i) {
      A(i, i) += ss;
    }
  }

  return A;
}

arma::vec normalize_weights(const arma::vec& w,
                            const std::string& mode,
                            double tol,
                            bool do_checks) {
  if (do_checks) {
    if (!(mode == "absolute" || mode == "relative")) {
      Rcpp::stop("mode must be \"absolute\" or \"relative\"");
    }
    if (!w.is_finite()) Rcpp::stop("weights must be finite.");
    if (!std::isfinite(tol) || tol <= 0.0) {
      Rcpp::stop("tol must be positive and finite.");
    }
  }

  const double s  = arma::accu(w);
  const double d1 = arma::norm(w, 1);

  const double denom = (mode == "absolute")
    ? std::max(1e-10, std::max(std::abs(s), tol))
    : std::max(1e-10, std::max(std::abs(s), tol * d1));

  if (!(denom > 0.0) || !std::isfinite(denom)) {
    return arma::vec(w.n_elem, arma::fill::zeros);
  }

  return w / denom;
}

// [[Rcpp::export]]
arma::mat prep_covariance_cpp(const arma::mat& sigma,
                              double epsilon,
                              bool stabilize) {
  return prep_covariance(sigma, epsilon, stabilize);
}

// [[Rcpp::export]]
arma::vec normalize_weights_cpp(const arma::vec& w,
                                const std::string& mode,
                                double tol,
                                bool do_checks) {
  return normalize_weights(w, mode, tol, do_checks);
}

// [[Rcpp::export]]
double eps_ridge_cpp() {
  return EPS_RIDGE;
}

// [[Rcpp::export]]
arma::mat safe_chol_cpp(const arma::mat& Q,
                        double base_bump,
                        unsigned int max_tries) {
  if (Q.n_rows != Q.n_cols) {
    Rcpp::stop("Matrix must be square for safe_chol_cpp.");
  }
  arma::mat A = Q;
  const arma::uword n = A.n_rows;
  const double trA = arma::trace(A);
  const double tau0 = base_bump * std::max(trA / std::max<double>(n, 1.0), arma::datum::eps);
  for (unsigned int t = 0; t <= max_tries; ++t) {
    const double bump = tau0 * std::pow(2.0, static_cast<int>(t));
    A.diag() += bump;
    arma::mat U;
    if (arma::chol(U, A, "upper")) {
      return U;
    }
  }
  Rcpp::stop("safe_chol_cpp failed after retries.");
}

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

// [[Rcpp::export]]
Rcpp::List design_from_moments_cpp(const arma::vec& mu,
                                   const arma::mat& sigma,
                                   double n_obs,
                                   double epsilon,
                                   bool stabilize_sigma) {
  const arma::uword N = mu.n_elem;
  if (sigma.n_rows != sigma.n_cols || sigma.n_rows != N) {
    Rcpp::stop("sigma must be square and match length(mu).");
  }
  if (!(n_obs > 0.0) || !std::isfinite(n_obs)) {
    Rcpp::stop("n_obs must be positive and finite.");
  }

  const arma::mat sigma_s = prep_covariance(sigma, epsilon, stabilize_sigma);
  arma::mat Q = n_obs * (sigma_s + mu * mu.t());

  // Small ridge to keep chol happy if needed
  double muQ = arma::mean(Q.diag());
  double tau = std::numeric_limits<double>::epsilon() * (std::isfinite(muQ) && muQ > 0.0 ? muQ : 1.0);
  Q.diag() += tau;

  arma::mat U;
  if (!arma::chol(U, Q)) {
    U = safe_chol_cpp(Q, 1e-10, 8);
  }

  arma::mat X = U.t();
  arma::vec y = arma::solve(U, n_obs * mu);

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
