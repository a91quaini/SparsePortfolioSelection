// miqp_heuristic_search.cpp
// Port of MIQPHeuristicSearch.jl using the Gurobi C++ API.

#include "miqp_heuristic_search.h"
#include "utils.h"
#include "sharpe_ratio.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include <gurobi_c++.h>

namespace {
struct SolveResult {
  arma::vec x;
  arma::uvec v;
  double sr;
  std::string status;
  double obj;
};

inline SolveResult fallback(std::size_t n, const std::string& status = "ERROR") {
  return {arma::vec(n, arma::fill::zeros),
          arma::uvec(),
          0.0,
          status,
          std::numeric_limits<double>::quiet_NaN()};
}

void sanitize_bounds(arma::vec& fmin, arma::vec& fmax) {
  const double dmin = 0.0;
  const double dmax = 1.0;
  for (arma::uword i = 0; i < fmin.n_elem; ++i) {
    double a = fmin[i];
    double b = fmax[i];
    if (!std::isfinite(a)) a = dmin;
    if (!std::isfinite(b)) b = dmax;
    if (a > b) std::swap(a, b);
    fmin[i] = a;
    fmax[i] = b;
  }
}

bool expand_bounds(const arma::vec& x,
                   const arma::uvec& v,
                   arma::vec& fmin,
                   arma::vec& fmax,
                   double factor,
                   double tol) {
  bool touched = false;
  for (arma::uword i = 0; i < x.n_elem; ++i) {
    if (v[i] == 1 &&
        (std::abs(x[i] - fmin[i]) <= tol || std::abs(x[i] - fmax[i]) <= tol)) {
      const double range = fmax[i] - fmin[i];
      const double delta = 0.5 * (factor - 1.0) * range;
      fmin[i] -= delta;
      fmax[i] += delta;
      touched = true;
    }
  }
  return touched;
}

std::string status_to_string(int status) {
  switch (status) {
    case GRB_OPTIMAL: return "OPTIMAL";
    case GRB_SUBOPTIMAL: return "SUBOPTIMAL";
    case GRB_INF_OR_UNBD: return "INF_OR_UNBD";
    case GRB_INFEASIBLE: return "INFEASIBLE";
    case GRB_UNBOUNDED: return "UNBOUNDED";
    case GRB_CUTOFF: return "CUTOFF";
    case GRB_ITERATION_LIMIT: return "ITERATION_LIMIT";
    case GRB_NODE_LIMIT: return "NODE_LIMIT";
    case GRB_TIME_LIMIT: return "TIME_LIMIT";
    case GRB_SOLUTION_LIMIT: return "SOLUTION_LIMIT";
    case GRB_INTERRUPTED: return "INTERRUPTED";
    case GRB_NUMERIC: return "NUMERIC";
    default: return "STATUS_" + std::to_string(status);
  }
}

void check_inputs(const arma::vec& mu,
                  const arma::mat& sigma,
                  unsigned int k,
                  unsigned int m_eff,
                  const arma::vec& fmin,
                  const arma::vec& fmax,
                  double gamma,
                  unsigned int expand_rounds,
                  double expand_factor,
                  double expand_tol,
                  double mipgap,
                  int threads,
                  bool normalize_weights,
                  double epsilon) {
  const arma::uword n = mu.n_elem;
  if (sigma.n_rows != sigma.n_cols) Rcpp::stop("sigma must be square.");
  if (sigma.n_rows != n) Rcpp::stop("sigma must be n x n.");
  if (k < 1 || k > n) Rcpp::stop("1 <= k <= n required.");
  if (m_eff > k) Rcpp::stop("m must be <= k.");
  if (fmin.n_elem != n || fmax.n_elem != n) Rcpp::stop("fmin/fmax must have length n.");
  if (gamma <= 0.0) Rcpp::stop("gamma must be positive.");
  if (expand_factor <= 0.0) Rcpp::stop("expand_factor must be > 0.");
  if (expand_tol < 0.0) Rcpp::stop("expand_tol must be >= 0.");
  if (mipgap < 0.0) Rcpp::stop("mipgap must be >= 0.");
  if (threads < 0) Rcpp::stop("threads must be >= 0.");
  if (!mu.is_finite() || !sigma.is_finite()) Rcpp::stop("mu/sigma must be finite.");
  if (!std::isfinite(epsilon)) Rcpp::stop("epsilon must be finite.");

  if (normalize_weights) {
    double smin = 0.0, smax = 0.0;
    for (arma::uword i = 0; i < n; ++i) {
      smin += std::min(0.0, fmin[i]);
      smax += std::max(0.0, fmax[i]);
    }
    if (!(1.0 >= smin - 1e-12 && 1.0 <= smax + 1e-12)) {
      Rcpp::stop("Budget sum(x)=1 incompatible with caps (quick screen).");
    }
  }
}

SolveResult solve_once(const arma::vec& mu,
                       const arma::mat& Sigma_eff,
                       unsigned int k,
                       unsigned int m_eff,
                       double gamma,
                       arma::vec fmin,
                       arma::vec fmax,
                       bool budget_constraint,
                       bool exactly_k,
                       double mipgap,
                       double time_limit,
                       int threads,
                       bool verbose,
                       double epsilon,
                       Rcpp::Nullable<arma::vec> x_start,
                       Rcpp::Nullable<arma::uvec> v_start) {
  const arma::uword n = mu.n_elem;
  sanitize_bounds(fmin, fmax);

  try {
    GRBEnv env = GRBEnv(true);
    env.set(GRB_IntParam_LogToConsole, verbose ? 1 : 0);
    env.start();

    GRBModel model(env);
    std::vector<GRBVar> x(n);
    std::vector<GRBVar> v(n);

    // Decision variables
    for (arma::uword i = 0; i < n; ++i) {
      x[i] = model.addVar(fmin[i], fmax[i], 0.0, GRB_CONTINUOUS);
      v[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
    }

    // Budget constraint
    if (budget_constraint) {
      GRBLinExpr sum = 0.0;
      for (arma::uword i = 0; i < n; ++i) sum += x[i];
      model.addConstr(sum == 1.0);
    }

    // Cardinality
    GRBLinExpr card = 0.0;
    for (arma::uword i = 0; i < n; ++i) card += v[i];
    if (exactly_k) {
      model.addConstr(card == static_cast<int>(k));
    } else {
      model.addConstr(card >= static_cast<int>(m_eff));
      model.addConstr(card <= static_cast<int>(k));
    }

    // Linking constraints
    for (arma::uword i = 0; i < n; ++i) {
      model.addConstr(x[i] <= fmax[i] * v[i]);
      model.addConstr(x[i] >= fmin[i] * v[i]);
    }

    // Objective: 0.5 * gamma * x' Sigma x - mu' x
    GRBQuadExpr obj = 0.0;
    for (arma::uword i = 0; i < n; ++i) {
      obj -= mu[i] * x[i];
    }
    for (arma::uword i = 0; i < n; ++i) {
      for (arma::uword j = 0; j < n; ++j) {
        const double cij = 0.5 * gamma * Sigma_eff(i, j);
        if (std::abs(cij) > 0.0) {
          obj += cij * x[i] * x[j];
        }
      }
    }
    model.setObjective(obj, GRB_MINIMIZE);

    // Parameters
    model.set(GRB_DoubleParam_MIPGap, mipgap);
    model.set(GRB_IntParam_NumericFocus, 1);
    if (threads > 0) model.set(GRB_IntParam_Threads, threads);
    model.set(GRB_DoubleParam_TimeLimit, time_limit);
    model.set(GRB_DoubleParam_FeasibilityTol, 1e-9);
    model.set(GRB_DoubleParam_OptimalityTol, 1e-9);
    model.set(GRB_DoubleParam_IntFeasTol, 1e-9);

    // Warm starts
    if (x_start.isNotNull() || v_start.isNotNull()) {
      if (x_start.isNotNull()) {
        arma::vec xs = Rcpp::as<arma::vec>(x_start);
        if (xs.n_elem == n) {
          for (arma::uword i = 0; i < n; ++i) {
            x[i].set(GRB_DoubleAttr_Start, xs[i]);
          }
        }
      }
      if (v_start.isNotNull()) {
        arma::uvec vs = Rcpp::as<arma::uvec>(v_start);
        if (vs.n_elem == n) {
          for (arma::uword i = 0; i < n; ++i) {
            v[i].set(GRB_DoubleAttr_Start, static_cast<int>(vs[i]));
          }
        }
      }
    }

    model.optimize();

    const int status = model.get(GRB_IntAttr_Status);
    const int sol_count = model.get(GRB_IntAttr_SolCount);
    if (sol_count == 0) return fallback(n, status_to_string(status));

    arma::vec x_out(n, arma::fill::zeros);
    arma::uvec v_out(n, arma::fill::zeros);
    for (arma::uword i = 0; i < n; ++i) {
      double xi = x[i].get(GRB_DoubleAttr_X);
      double vi = v[i].get(GRB_DoubleAttr_X);
      const unsigned int vi_bin = (std::isfinite(vi) && vi >= 0.5) ? 1u : 0u;
      v_out[i] = vi_bin;
      if (vi_bin == 0 && std::abs(xi) <= 1e-12) xi = 0.0;
      x_out[i] = std::isfinite(xi) ? xi : 0.0;
    }

    const double sr = compute_sr_cpp(x_out, mu, Sigma_eff,
                                     arma::uvec(), // all assets
                                     epsilon,
                                     /*stabilize_sigma=*/false,
                                     /*do_checks=*/false);
    const double obj_val = (status == GRB_INFEASIBLE || status == GRB_INF_OR_UNBD)
      ? std::numeric_limits<double>::quiet_NaN()
      : model.get(GRB_DoubleAttr_ObjVal);
    return {x_out, v_out, sr, status_to_string(status), obj_val};
  } catch (const GRBException& e) {
    return fallback(mu.n_elem, "ERROR:" + e.getMessage());
  } catch (...) {
    return fallback(mu.n_elem, "ERROR");
  }
}

} // anonymous namespace

// -----------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List mve_miqp_search_cpp(const arma::vec& mu,
                                         const arma::mat& sigma,
                                         unsigned int k,
                                         bool exactly_k,
                                         int m,
                                         double gamma,
                                         const arma::vec& fmin,
                                         const arma::vec& fmax,
                                         unsigned int expand_rounds,
                                         double expand_factor,
                                         double expand_tol,
                                         double mipgap,
                                         double time_limit,
                                         int threads,
                                         Rcpp::Nullable<arma::vec> x_start,
                                         Rcpp::Nullable<arma::uvec> v_start,
                                         bool compute_weights,
                                         bool normalize_wts,
                                         bool use_refit,
                                         bool verbose,
                                         double epsilon,
                                         bool stabilize_sigma,
                                         bool do_checks) {
  const arma::uword n = mu.n_elem;
  const unsigned int m_eff = exactly_k ? k :
    (m >= 0 ? static_cast<unsigned int>(m) : static_cast<unsigned int>(k > 0 ? k - 1 : 0));

  if (do_checks) {
    check_inputs(mu, sigma, k, m_eff, fmin, fmax, gamma,
                 expand_rounds, expand_factor, expand_tol,
                 mipgap, threads, normalize_wts, epsilon);
  }

  arma::vec fmin_work = fmin;
  arma::vec fmax_work = fmax;
  sanitize_bounds(fmin_work, fmax_work);

  // Prepare covariance once
  const arma::mat Sigma_eff = prep_covariance(sigma, epsilon, stabilize_sigma);

  SolveResult sol = solve_once(mu, Sigma_eff, k, m_eff, gamma,
                               fmin_work, fmax_work,
                               /*budget_constraint=*/normalize_wts,
                               exactly_k, mipgap, time_limit, threads,
                               verbose, epsilon, x_start, v_start);

  // Progressive bound expansion if needed
  for (unsigned int it = 0; it < expand_rounds; ++it) {
    if (!expand_bounds(sol.x, sol.v, fmin_work, fmax_work, expand_factor, expand_tol)) {
      break;
    }
    sol = solve_once(mu, Sigma_eff, k, m_eff, gamma,
                     fmin_work, fmax_work,
                     /*budget_constraint=*/normalize_wts,
                     exactly_k, mipgap, time_limit, threads,
                     verbose, epsilon, x_start, v_start);
  }

  arma::uvec sel = arma::find(sol.v == 1u);

  // Non-refit branch: return MIQP portfolio
  if (!use_refit) {
    arma::vec w = compute_weights ? sol.x : arma::vec(n, arma::fill::zeros);
    if (compute_weights && normalize_wts) {
      w = ::normalize_weights(w, "relative", 1e-6, false);
    }
    const double sr_out = (std::isfinite(sol.sr) ? sol.sr : 0.0);
    return Rcpp::List::create(
      Rcpp::Named("selection") = sel,
      Rcpp::Named("weights")   = w,
      Rcpp::Named("sr")        = sr_out,
      Rcpp::Named("status")    = sol.status
    );
  }

  // Refit branch: compute exact MVE on support
  if (sel.n_elem == 0) {
    return Rcpp::List::create(
      Rcpp::Named("selection") = sel,
      Rcpp::Named("weights")   = arma::vec(n, arma::fill::zeros),
      Rcpp::Named("sr")        = 0.0,
      Rcpp::Named("status")    = sol.status
    );
  }

  const double sr_refit = compute_mve_sr_cpp(mu, Sigma_eff, sel,
                                             epsilon,
                                             /*stabilize_sigma=*/false,
                                             /*do_checks=*/false);
  arma::vec w_refit = compute_weights
    ? compute_mve_weights_cpp(mu, Sigma_eff, sel,
                              normalize_wts,
                              epsilon,
                              /*stabilize_sigma=*/false,
                              /*do_checks=*/false)
    : arma::vec(n, arma::fill::zeros);
  if (compute_weights && !(arma::norm(w_refit, 1) > 1e-12) && !(std::abs(arma::accu(w_refit)) > 1e-12)) {
    w_refit.zeros();
  }

  return Rcpp::List::create(
    Rcpp::Named("selection") = sel,
    Rcpp::Named("weights")   = w_refit,
    Rcpp::Named("sr")        = std::isfinite(sr_refit) ? sr_refit : 0.0,
    Rcpp::Named("status")    = sol.status
  );
}
