// exhaustive_search.cpp
// Port of Julia's ExhaustiveSearch module (SparseMaxSR) to RcppArmadillo.

#include "exhaustive_search.h"
#include "utils.h"
#include "sharpe_ratio.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <unordered_set>
#include <vector>

using std::numeric_limits;

namespace {
constexpr double SR_TOL = 0.0;  // deterministic tie-break

inline void check_inputs(const arma::vec& mu,
                         const arma::mat& sigma,
                         unsigned int k,
                         bool enumerate_all,
                         unsigned int max_samples,
                         double epsilon) {
  const arma::uword n = mu.n_elem;
  if (n == 0) Rcpp::stop("mu must be non-empty.");
  if (sigma.n_rows != sigma.n_cols) Rcpp::stop("sigma must be square.");
  if (sigma.n_rows != n) Rcpp::stop("sigma must be n x n.");
  if (k < 1 || k > n) Rcpp::stop("k must be between 1 and n.");
  if (!mu.is_finite() || !sigma.is_finite()) Rcpp::stop("Non-finite entries in mu or sigma.");
  if (!std::isfinite(epsilon)) Rcpp::stop("epsilon must be finite.");
  if (!enumerate_all && max_samples == 0) {
    Rcpp::stop("When enumerate_all=false, max_samples must be > 0.");
  }
}

inline arma::uword draw_index(arma::uword n) {
  // Use R RNG for reproducibility.
  return static_cast<arma::uword>(std::floor(R::runif(0.0, 1.0) * static_cast<double>(n)));
}

// Sample a sorted k-subset from 0..n-1 without replacement.
void rand_k_subset(arma::uword n, arma::uword k, arma::uvec& out) {
  std::unordered_set<arma::uword> seen;
  seen.reserve(k * 2);
  arma::uword i = 0;
  while (i < k) {
    const arma::uword x = draw_index(n);
    if (seen.insert(x).second) {
      out[i++] = x;
    }
  }
  std::sort(out.begin(), out.end());
}

// Convert a subset to a UInt64 bitmask (n <= 64).
inline uint64_t to_mask(const arma::uvec& sel) {
  uint64_t key = 0;
  for (arma::uword v : sel) {
    key |= (uint64_t(1) << v);
  }
  return key;
}

// Score closure: assumes Sigma_eff is already stabilized; skips further prep.
class Scorer {
 public:
  Scorer(const arma::vec& mu, const arma::mat& Sigma_eff, double epsilon)
      : mu_(mu), Sigma_eff_(Sigma_eff), epsilon_(epsilon) {}

  double operator()(const arma::uvec& sel) const {
    return compute_mve_sr_cpp(mu_, Sigma_eff_, sel,
                              epsilon_,
                              /*stabilize_sigma=*/false,
                              /*do_checks=*/false);
  }

 private:
  const arma::vec& mu_;
  const arma::mat& Sigma_eff_;
  double epsilon_;
};

// Exhaustive enumeration over all combinations.
std::pair<double, arma::uvec> enumerate_best(arma::uword n,
                                             arma::uword k,
                                             const Scorer& scorer) {
  double best_sr = -numeric_limits<double>::infinity();
  arma::uvec best_sel;

  arma::uvec buf(k, arma::fill::zeros);
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
                                          const Scorer& scorer) {
  double best_sr = -numeric_limits<double>::infinity();
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
          h ^= std::hash<arma::uword>()(x + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2));
        }
        return h;
      }
    };
    std::unordered_set<std::vector<arma::uword>, VecHash> seen;
    seen.reserve(m * 2);
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
                                     double epsilon,
                                     bool stabilize_sigma,
                                     bool do_checks,
                                     bool enumerate_all,
                                     unsigned int max_samples,
                                     bool dedup_samples,
                                     bool compute_weights) {
  Rcpp::RNGScope scope; // ensure R RNG is scoped for sampling paths

  if (do_checks) {
    check_inputs(mu, sigma, k, enumerate_all, max_samples, epsilon);
  }

  const arma::uword n = mu.n_elem;

  // One-time stabilization of covariance.
  const arma::mat Sigma_eff = prep_covariance(sigma, epsilon, stabilize_sigma);

  // Build scorer that reuses Sigma_eff (stabilize_sigma=false inside).
  const Scorer scorer(mu, Sigma_eff, epsilon);

  double best_sr = -numeric_limits<double>::infinity();
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
    weights = compute_mve_weights_cpp(mu, Sigma_eff, best_sel,
                                      /*normalize_w=*/false,
                                      epsilon,
                                      /*stabilize_sigma=*/false,
                                      /*do_checks=*/false);
  }

  return Rcpp::List::create(
    Rcpp::Named("selection") = best_sel,
    Rcpp::Named("weights")   = weights,
    Rcpp::Named("sr")        = best_sr,
    Rcpp::Named("status")    = status
  );
}
