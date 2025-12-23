#include "miqp_heuristic_search.h"

#include <RcppArmadillo.h>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include <gurobi_c++.h>

#include "utils.h"
#include "sharpe_ratio.h" // for compute_mve_weights_cpp if you want refit weights

namespace {

struct SolveResult {
  arma::vec x;        // continuous weights
  arma::uvec v;       // 0/1 selection
  std::string status; // Gurobi status string
  double obj;         // objective value (min form)
};

  inline std::string grb_status_to_string(int status) {
    switch (status) {
    case GRB_OPTIMAL:          return "OPTIMAL";
    case GRB_SUBOPTIMAL:       return "SUBOPTIMAL";
    case GRB_INFEASIBLE:       return "INFEASIBLE";
    case GRB_UNBOUNDED:        return "UNBOUNDED";
    case GRB_CUTOFF:           return "CUTOFF";
    case GRB_ITERATION_LIMIT:  return "ITERATION_LIMIT";
    case GRB_NODE_LIMIT:       return "NODE_LIMIT";
    case GRB_TIME_LIMIT:       return "TIME_LIMIT";
    case GRB_SOLUTION_LIMIT:   return "SOLUTION_LIMIT";
    case GRB_INTERRUPTED:      return "INTERRUPTED";
    case GRB_NUMERIC:          return "NUMERIC";
    default:                   return "STATUS_" + std::to_string(status);
    }
  }

  inline SolveResult fallback(std::size_t n, const std::string& status) {
    SolveResult out;
    out.x = arma::vec(n, arma::fill::zeros);
    out.v = arma::uvec(n, arma::fill::zeros);
    out.status = status;
    out.obj = std::numeric_limits<double>::quiet_NaN();
    return out;
  }

  inline void sanitize_bounds(arma::vec& fmin, arma::vec& fmax) {
    // Replace non-finite bounds with a conservative default.
    // (You can tighten this if you want to enforce fmin<=0<=fmax.)
    for (arma::uword i = 0; i < fmin.n_elem; ++i) {
      if (!std::isfinite(fmin[i])) fmin[i] = 0.0;
      if (!std::isfinite(fmax[i])) fmax[i] = 1.0;
      if (fmin[i] > fmax[i]) std::swap(fmin[i], fmax[i]);
    }
  }

  inline void check_inputs(const arma::vec& mu,
                           const arma::mat& sigma,
                           unsigned int k,
                           int m,
                           double gamma,
                           const arma::vec& fmin,
                           const arma::vec& fmax,
                           bool normalize_weights,
                           double ridge_epsilon,
                           double expand_factor,
                           double expand_tol,
                           double mipgap,
                           double time_limit,
                           int threads) {
    const arma::uword n = mu.n_elem;
    if (n == 0) Rcpp::stop("mu must be non-empty.");
    if (sigma.n_rows != sigma.n_cols) Rcpp::stop("sigma must be square.");
    if (sigma.n_rows != n) Rcpp::stop("sigma must be n x n.");
    if (k < 1 || k > n) Rcpp::stop("k must satisfy 1 <= k <= n.");
    if (!mu.is_finite() || !sigma.is_finite()) Rcpp::stop("mu/sigma must be finite.");
    if (fmin.n_elem != n || fmax.n_elem != n) Rcpp::stop("fmin/fmax must have length n.");
    if (!(gamma > 0.0) || !std::isfinite(gamma)) Rcpp::stop("gamma must be positive and finite.");
    if (!std::isfinite(ridge_epsilon) || ridge_epsilon < 0.0) Rcpp::stop("ridge_epsilon must be >= 0 and finite.");
    if (!std::isfinite(expand_factor) || expand_factor <= 0.0) Rcpp::stop("expand_factor must be > 0 and finite.");
    if (!std::isfinite(expand_tol) || expand_tol < 0.0) Rcpp::stop("expand_tol must be >= 0 and finite.");
    if (!std::isfinite(mipgap) || mipgap < 0.0) Rcpp::stop("mipgap must be >= 0 and finite.");
    if (!std::isfinite(time_limit) || time_limit <= 0.0) Rcpp::stop("time_limit must be > 0 and finite.");
    if (threads < 0) Rcpp::stop("threads must be >= 0.");

    if (m > static_cast<int>(k)) Rcpp::stop("m must be <= k.");
    if (m < 0) Rcpp::stop("m must be >= 0.");

    // Quick feasibility screen if we impose sum(w)=1.
    // This ignores selection/binaries, but catches obvious inconsistencies.
    if (normalize_weights) {
      double smin = 0.0, smax = 0.0;
      for (arma::uword i = 0; i < n; ++i) {
        smin += std::min(0.0, fmin[i]);
        smax += std::max(0.0, fmax[i]);
      }
      if (!(1.0 >= smin - 1e-10 && 1.0 <= smax + 1e-10)) {
        Rcpp::stop("Budget constraint sum(w)=1 appears incompatible with fmin/fmax (quick screen).");
      }
    }
  }

  inline arma::mat prepare_sigma(const arma::mat& sigma, double ridge_epsilon) {
    // No hidden regularization:
    // - always symmetrize
    // - only add ridge if user requests ridge_epsilon > 0
    arma::mat S = arma::symmatu(sigma);
    if (ridge_epsilon > 0.0) {
      S.diag() += ridge_epsilon;
    }
    return S;
  }

  inline double portfolio_sr(const arma::vec& w,
                             const arma::vec& mu,
                             const arma::mat& Sigma) {
    // SR(w) = (w'mu)/sqrt(w'Sigma w)
    const double num = arma::dot(w, mu);
    const double den2 = arma::as_scalar(w.t() * Sigma * w);
    if (!(den2 > 0.0) || !std::isfinite(den2)) return std::numeric_limits<double>::quiet_NaN();
    const double den = std::sqrt(den2);
    const double sr = num / den;
    return std::isfinite(sr) ? sr : std::numeric_limits<double>::quiet_NaN();
  }

  inline bool expand_bounds_multiplicative(const arma::vec& x,
                                           const arma::uvec& v,
                                           arma::vec& fmin,
                                           arma::vec& fmax,
                                           double factor,
                                           double tol) {
    // Expand only bounds that are active, without changing sign conventions.
    bool touched = false;
    for (arma::uword i = 0; i < x.n_elem; ++i) {
      if (v[i] != 1u) continue;

      if (std::abs(x[i] - fmin[i]) <= tol) {
        const double a = fmin[i];
        if (a < 0.0) fmin[i] = a * factor;          // more negative
        else if (a > 0.0) fmin[i] = a / factor;     // closer to 0
        // if a == 0: stays 0
        touched = true;
      }

      if (std::abs(x[i] - fmax[i]) <= tol) {
        const double b = fmax[i];
        if (b > 0.0) fmax[i] = b * factor;          // more positive
        else if (b < 0.0) fmax[i] = b / factor;     // closer to 0 (larger max)
        // if b == 0: stays 0
        touched = true;
      }

      if (fmin[i] > fmax[i]) std::swap(fmin[i], fmax[i]);
    }
    return touched;
  }

  inline SolveResult solve_once(const arma::vec& mu,
                                const arma::mat& Sigma_eff,
                                unsigned int k,
                                unsigned int m_eff,
                                double gamma,
                                const arma::vec& fmin,
                                const arma::vec& fmax,
                                bool budget_constraint,
                                bool exactly_k,
                                double mipgap,
                                double time_limit,
                                int threads,
                                bool verbose,
                                const arma::vec* x_start,
                                const arma::uvec* v_start) {
    const arma::uword n = mu.n_elem;

    try {
      GRBEnv env(true);
      env.set(GRB_IntParam_LogToConsole, verbose ? 1 : 0);
      env.start();

      GRBModel model(env);

      std::vector<GRBVar> x(n);
      std::vector<GRBVar> v(n);

      // Variables
      for (arma::uword i = 0; i < n; ++i) {
        x[i] = model.addVar(fmin[i], fmax[i], 0.0, GRB_CONTINUOUS);
        v[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
      }

      // Warm start (optional)
      if (x_start != nullptr && x_start->n_elem == n) {
        for (arma::uword i = 0; i < n; ++i) {
          x[i].set(GRB_DoubleAttr_Start, (*x_start)[i]);
        }
      }
      if (v_start != nullptr && v_start->n_elem == n) {
        for (arma::uword i = 0; i < n; ++i) {
          v[i].set(GRB_DoubleAttr_Start, ((*v_start)[i] ? 1.0 : 0.0));
        }
      }

      // Budget constraint (only when requested)
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

      // Objective (variance-based):
      // maximize  mu'x - (gamma/2) x'Sigma x
      // => minimize -mu'x + (gamma/2) x'Sigma x
      GRBQuadExpr obj = 0.0;

      // Linear part: -mu'x
      for (arma::uword i = 0; i < n; ++i) {
        obj += (-mu[i]) * x[i];
      }

      // Quadratic part: (gamma/2) x' Sigma x
      // Build without double-counting: diag uses (gamma/2)*S_ii; off-diag uses gamma*S_ij for i<j.
      for (arma::uword i = 0; i < n; ++i) {
        const double cii = 0.5 * gamma * Sigma_eff(i, i);
        if (cii != 0.0) obj += cii * x[i] * x[i];
        for (arma::uword j = i + 1; j < n; ++j) {
          const double cij = gamma * Sigma_eff(i, j); // because i<j accounts for both symmetric terms
          if (cij != 0.0) obj += cij * x[i] * x[j];
        }
      }

      model.setObjective(obj, GRB_MINIMIZE);

      // Params
      model.set(GRB_DoubleParam_MIPGap, mipgap);
      model.set(GRB_DoubleParam_TimeLimit, time_limit);
      if (threads > 0) model.set(GRB_IntParam_Threads, threads);

      // Keep tolerances reasonable (too tight can slow MIQP a lot)
      model.set(GRB_DoubleParam_FeasibilityTol, 1e-8);
      model.set(GRB_DoubleParam_OptimalityTol, 1e-8);
      model.set(GRB_DoubleParam_IntFeasTol,   1e-8);

      model.optimize();

      const int status = model.get(GRB_IntAttr_Status);
      const std::string status_s = grb_status_to_string(status);

      // If no incumbent solution, return fallback
      if (model.get(GRB_IntAttr_SolCount) <= 0) {
        return fallback(n, status_s);
      }

      SolveResult out;
      out.x = arma::vec(n, arma::fill::zeros);
      out.v = arma::uvec(n, arma::fill::zeros);
      for (arma::uword i = 0; i < n; ++i) {
        out.x[i] = x[i].get(GRB_DoubleAttr_X);
        out.v[i] = (v[i].get(GRB_DoubleAttr_X) > 0.5) ? 1u : 0u;
      }
      out.status = status_s;
      out.obj = model.get(GRB_DoubleAttr_ObjVal);
      return out;

    } catch (GRBException& e) {
      return fallback(n, std::string("GUROBI_ERROR_") + std::to_string(e.getErrorCode()));
    } catch (...) {
      return fallback(n, "GUROBI_ERROR");
    }
  }

} // anonymous namespace

// [[Rcpp::export]]
Rcpp::List mve_miqp_search_cpp(const arma::vec& mu,
                               const arma::mat& sigma,
                               unsigned int k,
                               const arma::vec& fmin,
                               const arma::vec& fmax,
                               Rcpp::Nullable<arma::vec> x_start,
                               Rcpp::Nullable<arma::uvec> v_start,
                               int m,
                               double gamma,
                               bool exactly_k,
                               unsigned int expand_rounds,
                               double expand_factor,
                               double expand_tol,
                               double mipgap,
                               double time_limit,
                               int threads,
                               double ridge_epsilon,
                               bool normalize_weights,
                               bool use_refit,
                               bool verbose,
                               bool do_checks) {
  if (do_checks) {
    check_inputs(mu, sigma, k, m, gamma, fmin, fmax,
                 normalize_weights, ridge_epsilon,
                 expand_factor, expand_tol, mipgap, time_limit, threads);
  }

  const arma::uword n = mu.n_elem;

  arma::vec fmin_work = fmin;
  arma::vec fmax_work = fmax;
  sanitize_bounds(fmin_work, fmax_work);

  const arma::mat Sigma_eff = prepare_sigma(sigma, ridge_epsilon);

  const unsigned int m_eff = (m <= 0 ? 1u : static_cast<unsigned int>(m));

  // If normalize_weights==true, we enforce the budget constraint in the MIQP,
  // and we also post-normalize with normalization_type=0.
  const bool budget_constraint = normalize_weights;

  // Prepare warm-start pointers once.
  arma::vec x_start_local;
  arma::uvec v_start_local;
  const arma::vec* x_start_ptr = nullptr;
  const arma::uvec* v_start_ptr = nullptr;
  if (x_start.isNotNull()) {
    x_start_local = Rcpp::as<arma::vec>(x_start);
    if (x_start_local.n_elem == n) x_start_ptr = &x_start_local;
  }
  if (v_start.isNotNull()) {
    v_start_local = Rcpp::as<arma::uvec>(v_start);
    if (v_start_local.n_elem == n) v_start_ptr = &v_start_local;
  }

  // First solve
  Rcpp::checkUserInterrupt();
  SolveResult sol = solve_once(mu, Sigma_eff, k, m_eff, gamma,
                               fmin_work, fmax_work,
                               budget_constraint, exactly_k,
                               mipgap, time_limit, threads,
                               verbose,
                               x_start_ptr, v_start_ptr);

  // Progressive bound expansion (warm start with previous solution)
  for (unsigned int it = 0; it < expand_rounds; ++it) {
    Rcpp::checkUserInterrupt();
    if (!expand_bounds_multiplicative(sol.x, sol.v, fmin_work, fmax_work,
                                      expand_factor, expand_tol)) {
      break;
    }
    sol = solve_once(mu, Sigma_eff, k, m_eff, gamma,
                     fmin_work, fmax_work,
                     budget_constraint, exactly_k,
                     mipgap, time_limit, threads,
                     verbose,
                     &sol.x, &sol.v);
  }

  const arma::uvec sel = arma::find(sol.v == 1u);

  // -----------------------------
  // Non-refit branch: return MIQP x
  // -----------------------------
  if (!use_refit) {
    arma::vec w_used = sol.x;

    if (normalize_weights) {
      // MUST be type=0 by your requirement
      w_used = normalize_weights_cpp(w_used, 1e-6, /*normalization_type=*/0);
    }

    const double sr_used = portfolio_sr(w_used, mu, Sigma_eff);

    arma::vec w_out = w_used;

    return Rcpp::List::create(
      Rcpp::Named("selection") = sel,
      Rcpp::Named("weights")   = w_out,
      Rcpp::Named("sr")        = std::isfinite(sr_used) ? sr_used : 0.0,
      Rcpp::Named("status")    = sol.status,
      Rcpp::Named("obj")       = sol.obj
    );
  }

  // -----------------------------
  // Refit branch: exact MVE on selected support
  // -----------------------------
  if (sel.n_elem == 0) {
    return Rcpp::List::create(
      Rcpp::Named("selection") = sel,
      Rcpp::Named("weights")   = arma::vec(n, arma::fill::zeros),
      Rcpp::Named("sr")        = 0.0,
      Rcpp::Named("status")    = sol.status,
      Rcpp::Named("obj")       = sol.obj
    );
  }

  // compute_mve_weights_cpp should solve Sigma_A w_A = mu_A and place into full vector.
  // If normalize_weights==true: MUST be normalization_type=0.
  arma::vec w_refit = compute_mve_weights_cpp(mu, Sigma_eff, sel,
                            /*ridge_epsilon=*/0.0,
                            /*normalize_weights=*/normalize_weights,
                            /*normalization_type=*/0,
                            /*do_checks=*/false);

  // Even if not returning weights, compute SR on the refit allocation actually used:
  arma::vec w_used = w_refit;
  w_used = compute_mve_weights_cpp(mu, Sigma_eff, sel,
                                   /*ridge_epsilon=*/0.0,
                                   /*normalize_weights=*/normalize_weights,
                                   /*normalization_type=*/0,
                                   /*do_checks=*/false);

  const double sr_refit = portfolio_sr(w_used, mu, Sigma_eff);

  return Rcpp::List::create(
    Rcpp::Named("selection") = sel,
    Rcpp::Named("weights")   = w_refit,
    Rcpp::Named("sr")        = sr_refit,
    Rcpp::Named("status")    = sol.status,
    Rcpp::Named("obj")       = sol.obj
  );
}
