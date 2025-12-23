// exhaustive_search.cpp

#include "exhaustive_search.h"
#include "utils.h"
#include "sharpe_ratio.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <unordered_set>
#include <vector>

namespace {

constexpr double SR_TOL = 0.0; // deterministic tie-break only

  inline void check_inputs(const arma::vec& mu,
                           const arma::mat& sigma,
                           unsigned int k,
                           bool enumerate_all,
                           unsigned int max_samples,
                           double ridge_epsilon) {
    const arma::uword n = mu.n_elem;
    if (n == 0) Rcpp::stop("mve_exhaustive_search_cpp: mu must be non-empty.");
    if (sigma.n_rows != sigma.n_cols) Rcpp::stop("mve_exhaustive_search_cpp: sigma must be square.");
    if (sigma.n_rows != n) Rcpp::stop("mve_exhaustive_search_cpp: sigma must be n x n.");
    if (k < 1 || k > n) Rcpp::stop("mve_exhaustive_search_cpp: k must be between 1 and n.");
    if (!mu.is_finite() || !sigma.is_finite()) Rcpp::stop("mve_exhaustive_search_cpp: non-finite inputs.");
    if (!std::isfinite(ridge_epsilon) || ridge_epsilon < 0.0)
      Rcpp::stop("mve_exhaustive_search_cpp: ridge_epsilon must be finite and nonnegative.");
    if (!enumerate_all && max_samples == 0) {
      Rcpp::stop("mve_exhaustive_search_cpp: when enumerate_all=false, max_samples must be > 0.");
    }
  }

  inline arma::uword draw_index(arma::uword n) {
    return static_cast<arma::uword>(std::floor(R::runif(0.0, 1.0) * static_cast<double>(n)));
  }

  // Sample a sorted k-subset from 0..n-1 without replacement.
  inline void rand_k_subset(arma::uword n, arma::uword k, arma::uvec& out) {
    std::unordered_set<arma::uword> seen;
    seen.reserve(static_cast<std::size_t>(k) * 2);

    arma::uword i = 0;
    while (i < k) {
      const arma::uword x = draw_index(n);
      if (seen.insert(x).second) {
        out[i++] = x;
      }
    }
    std::sort(out.begin(), out.end());
  }

  // Convert a subset to a UInt64 bitmask (only valid if n <= 64).
  inline uint64_t to_mask(const arma::uvec& sel) {
    uint64_t key = 0;
    for (arma::uword v : sel) key |= (uint64_t(1) << v);
    return key;
  }

  // Fast scorer: variance-based MVE Sharpe on subset A: sqrt(mu_A' Sigma_A^{-1} mu_A).
  // Assumes sigma_ has already been stabilized (or not) upstream; does not re-stabilize.
  class MveSharpeScorer {
  public:
    MveSharpeScorer(const arma::vec& mu, const arma::mat& sigma_)
      : mu_(mu), sigma_(sigma_) {}

    double operator()(const arma::uvec& sel) const {
      if (sel.n_elem == 0 || sel.n_elem == mu_.n_elem) {
        const arma::vec x = solve_sympd(sigma_, mu_);
        const double v = arma::dot(mu_, x);
        return std::sqrt(std::max(v, 0.0));
      }

      const arma::vec mu_s = mu_.elem(sel);
      const arma::mat sigma_s = sigma_.submat(sel, sel);
      const arma::vec x = solve_sympd(sigma_s, mu_s);
      const double v = arma::dot(mu_s, x);
      return std::sqrt(std::max(v, 0.0));
    }

  private:
    const arma::vec& mu_;
    const arma::mat& sigma_;
  };

  // Exhaustive enumeration over all combinations.
  std::pair<double, arma::uvec> enumerate_best(arma::uword n,
                                               arma::uword k,
                                               const MveSharpeScorer& scorer) {
    double best_sr = -std::numeric_limits<double>::infinity();
    arma::uvec best_sel;

    for_each_combination(static_cast<unsigned int>(n),
                         static_cast<unsigned int>(k),
                         [&](const arma::uvec& sel) {
                           const double sr = scorer(sel);
                           if (sr > best_sr + SR_TOL) {
                             best_sr = sr;
                             best_sel = sel;
                           }
                         });

    return {best_sr, best_sel};
  }

  // Sampling-based search with optional deduplication.
  std::pair<double, arma::uvec> sample_best(arma::uword n,
                                            arma::uword k,
                                            unsigned int m,
                                            bool dedup,
                                            const MveSharpeScorer& scorer) {
    double best_sr = -std::numeric_limits<double>::infinity();
    arma::uvec best_sel;
    arma::uvec buf(k, arma::fill::zeros);

    if (!dedup) {
      for (unsigned int i = 0; i < m; ++i) {
        rand_k_subset(n, k, buf);
        const double sr = scorer(buf);
        if (sr > best_sr + SR_TOL) {
          best_sr = sr;
          best_sel = buf;
        }
      }
      return {best_sr, best_sel};
    }

    if (n <= 64) {
      std::unordered_set<uint64_t> seen;
      seen.reserve(static_cast<std::size_t>(m) * 2);

      while (seen.size() < m) {
        rand_k_subset(n, k, buf);
        const uint64_t key = to_mask(buf);
        if (seen.insert(key).second) {
          const double sr = scorer(buf);
          if (sr > best_sr + SR_TOL) {
            best_sr = sr;
            best_sel = buf;
          }
        }
      }
    } else {
      struct VecHash {
        std::size_t operator()(const std::vector<arma::uword>& v) const noexcept {
          std::size_t h = 0;
          for (auto x : v) {
            // standard hash combine
            h ^= std::hash<arma::uword>()(x) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
          }
          return h;
        }
      };

      std::unordered_set<std::vector<arma::uword>, VecHash> seen;
      seen.reserve(static_cast<std::size_t>(m) * 2);

      while (seen.size() < m) {
        rand_k_subset(n, k, buf);
        const std::vector<arma::uword> key = arma::conv_to<std::vector<arma::uword>>::from(buf);
        if (seen.insert(key).second) {
          const double sr = scorer(buf);
          if (sr > best_sr + SR_TOL) {
            best_sr = sr;
            best_sel = buf;
          }
        }
      }
    }

    return {best_sr, best_sel};
  }

} // anonymous namespace

// -----------------------------------------------------------------------------
// Public API
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List mve_exhaustive_search_cpp(const arma::vec& mu,
                                     const arma::mat& sigma,
                                     unsigned int k,
                                     double ridge_epsilon,
                                     bool enumerate_all,
                                     unsigned int max_samples,
                                     bool dedup_samples,
                                     bool compute_weights,
                                     bool normalize_weights,
                                     int normalization_type,
                                     bool do_checks) {
  Rcpp::RNGScope scope; // reproducible sampling

  if (do_checks) {
    check_inputs(mu, sigma, k, enumerate_all, max_samples, ridge_epsilon);
  }

  const arma::uword n = mu.n_elem;

  // Stabilize once (or not at all, if ridge_epsilon==0).
  const arma::mat sigma_ = stabilize_sigma_cpp(sigma, ridge_epsilon);

  // Scorer reuses sigma_ without further stabilization.
  const MveSharpeScorer scorer(mu, sigma_);

  double best_sr = -std::numeric_limits<double>::infinity();
  arma::uvec best_sel;
  std::string status;

  if (enumerate_all) {
    std::tie(best_sr, best_sel) = enumerate_best(n, k, scorer);
    status = "EXHAUSTIVE";
  } else {
    std::tie(best_sr, best_sel) = sample_best(n, k, max_samples, dedup_samples, scorer);
    status = "SAMPLED";
  }

  arma::vec weights(n, arma::fill::zeros);
  if (compute_weights && best_sel.n_elem > 0) {
    // Use variance-based MVE weights: w = Sigma^{-1} mu on the selected subset.
    // sigma_ already includes any stabilization; pass ridge_epsilon=0 to avoid double-ridge.
    weights = compute_mve_weights_cpp(mu, sigma_, best_sel,
                                      /*ridge_epsilon=*/0.0,
                                      /*normalize_weights=*/normalize_weights,
                                      /*normalization_type=*/normalization_type,
                                      /*do_checks=*/false);
  }

  return Rcpp::List::create(
    Rcpp::Named("selection") = best_sel,
    Rcpp::Named("weights")   = weights,
    Rcpp::Named("sr")        = best_sr,
    Rcpp::Named("status")    = status
  );
}
