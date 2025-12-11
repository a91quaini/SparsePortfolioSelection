#' Exhaustive search over k-subsets maximizing MVE Sharpe ratio
#'
#' Calls the C++ exhaustive (or sampled) search to find the k-asset subset with
#' highest in-sample MVE Sharpe ratio. Covariance is symmetrized/ridged once via
#' `prep_covariance_cpp` and reused in the loop.
#'
#' @param mu Numeric mean vector (length \eqn{N}).
#' @param sigma Numeric covariance matrix (\eqn{N \times N}).
#' @param k Target cardinality (1..N).
#' @param epsilon Ridge factor for covariance stabilization.
#' @param stabilize_sigma Logical; if `TRUE`, symmetrize/ridge-stabilize `sigma`.
#' @param do_checks Logical; if `TRUE`, validate inputs.
#' @param enumerate_all Logical; `TRUE` to enumerate all \eqn{N \choose k};
#'   `FALSE` to sample supports.
#' @param max_samples Number of sampled supports when `enumerate_all = FALSE`.
#' @param dedup_samples Logical; ensure sampled supports are distinct.
#' @param compute_weights Logical; if `TRUE`, return MVE weights for the best subset.
#'
#' @return List with `selection` (1-based indices), `weights`, `sr`, `status`.
#' @export
mve_exhaustive_search <- function(mu, sigma, k,
                                  epsilon = eps_ridge_cpp(),
                                  stabilize_sigma = TRUE,
                                  do_checks = FALSE,
                                  enumerate_all = TRUE,
                                  max_samples = 0L,
                                  dedup_samples = TRUE,
                                  compute_weights = TRUE) {
  if (!enumerate_all && max_samples <= 0) {
    stop("When enumerate_all=FALSE, max_samples must be > 0.")
  }
  res <- mve_exhaustive_search_cpp(mu, sigma, k, epsilon, stabilize_sigma, do_checks,
                                   enumerate_all, max_samples, dedup_samples, compute_weights)
  res$selection <- as.integer(res$selection) + 1L
  res
}

#' MIQP heuristic search for sparse MVE selection (Gurobi backend)
#'
#' Solves a quadratic MIP heuristic to select at most `k` assets (or exactly `k`
#' if `exactly_k=TRUE`). Objective is \eqn{0.5 * gamma * x' Sigma x - mu' x} with
#' linking constraints \eqn{fmin_i v_i \le x_i \le fmax_i v_i}, cardinality
#' constraints, optional budget \eqn{\sum x_i = 1} via `normalize_weights`, and
#' progressive bound expansion when solutions sit on caps. Uses the C++ Gurobi
#' backend; supports optional refit of exact MVE on the selected support.
#'
#' @param mu Numeric mean vector.
#' @param sigma Numeric covariance matrix.
#' @param k Target cardinality.
#' @param exactly_k Logical; if `TRUE`, enforce \eqn{\sum v_i = k} else banded.
#' @param m Lower bound on \eqn{\sum v_i} when `exactly_k=FALSE`.
#' @param gamma Risk-aversion weight on the quadratic term.
#' @param fmin,fmax Lower/upper bounds for each weight (scalar recycled to length N).
#' @param expand_rounds Number of progressive bound-expansion re-solves.
#' @param expand_factor Multiplicative factor for expanding active bounds.
#' @param expand_tol Tolerance for detecting bound hits.
#' @param mipgap Relative MIP gap.
#' @param time_limit Time limit in seconds.
#' @param threads Number of Gurobi threads (0 lets Gurobi decide).
#' @param x_start Optional warm-start weights.
#' @param v_start Optional warm-start binaries.
#' @param compute_weights Logical; return weights (else zeros).
#' @param normalize_weights Logical; if `TRUE`, add budget constraint and
#'   post-scale refit weights.
#' @param use_refit Logical; if `TRUE`, refit exact MVE on the selected support.
#' @param verbose Logical; emit Gurobi output.
#' @param epsilon Ridge for covariance stabilization.
#' @param stabilize_sigma Logical; symmetrize/ridge-stabilize sigma.
#' @param do_checks Logical; validate inputs.
#'
#' @return List with `selection` (1-based), `weights`, `sr`, `status`.
#' @export
mve_miqp_search <- function(mu, sigma, k,
                            exactly_k = FALSE,
                            m = 1L,
                            gamma = 1.0,
                            fmin = -0.25,
                            fmax = 0.25,
                            expand_rounds = 20L,
                            expand_factor = 3.0,
                            expand_tol = 1e-2,
                            mipgap = 1e-4,
                            time_limit = 200,
                            threads = 1L,
                            x_start = NULL,
                            v_start = NULL,
                            compute_weights = TRUE,
                            normalize_weights = FALSE,
                            use_refit = FALSE,
                            verbose = FALSE,
                            epsilon = eps_ridge_cpp(),
                            stabilize_sigma = TRUE,
                            do_checks = FALSE) {
  fmin_vec <- if (length(fmin) == 1L) rep_len(fmin, length(mu)) else fmin
  fmax_vec <- if (length(fmax) == 1L) rep_len(fmax, length(mu)) else fmax
  res <- mve_miqp_search_cpp(mu, sigma, k, exactly_k, m, gamma, fmin_vec, fmax_vec,
                             expand_rounds, expand_factor, expand_tol,
                             mipgap, time_limit, threads,
                             x_start, v_start,
                             compute_weights, normalize_weights,
                             use_refit, verbose,
                             epsilon, stabilize_sigma, do_checks)
  # Convert 0-based selection from C++ to 1-based R
  res$selection <- as.integer(res$selection) + 1L
  res
}

#' LASSO relaxation heuristic for sparse MVE selection (moment-based)
#'
#' Variant that takes precomputed moments (`mu`, `sigma`) and an assumed sample
#' size `n_obs`, builds the synthetic design directly, and runs the same glmnet
#' selection as the return-based version.
#'
#' @param mu Numeric mean vector.
#' @param sigma Numeric covariance matrix.
#' @param n_obs Sample size underlying the moments.
#' @inheritParams mve_lasso_search_from_returns
#'
#' @return List with `selection` (1-based), `weights`, `sr`, `status`, `lambda`, `alpha`.
#' @export
mve_lasso_search <- function(
    mu,
    sigma,
    n_obs,
    k,
    nlambda = 100L,
    lambda_min_ratio = 1e-3,
    lambda = NULL,
    alpha = 1,
    n_folds = 5L,
    nadd = 80L,
    nnested = 2L,
    standardize = FALSE,
    epsilon = eps_ridge_cpp(),
    stabilize_sigma = TRUE,
    compute_weights = TRUE,
    normalize_weights = FALSE,
    use_refit = FALSE,
    do_checks = FALSE,
    R_cv = NULL
) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is required for mve_lasso_search().")
  }

  if (do_checks) {
    if (!is.numeric(mu) || length(mu) == 0) stop("mu must be a numeric vector.")
    if (!is.matrix(sigma) || !is.numeric(sigma)) stop("sigma must be a numeric matrix.")
    if (nrow(sigma) != ncol(sigma) || nrow(sigma) != length(mu)) stop("sigma must be square and match length(mu).")
    if (!is.numeric(n_obs) || length(n_obs) != 1 || n_obs <= 0) stop("n_obs must be a positive scalar.")
    if (!is.numeric(k) || k < 0 || k > length(mu)) stop("k must be between 0 and length(mu).")
    if (!is.null(lambda)) {
      if (any(!is.finite(lambda)) || any(lambda <= 0)) stop("lambda must be positive/finite.")
    }
    if (!is.finite(epsilon)) stop("epsilon must be finite.")
  }

  N <- length(mu)

  alpha_info <- .validate_alpha_grid(alpha)
  alpha_grid <- alpha_info$grid
  alpha_used <- alpha_grid[1]
  if (alpha_info$is_grid && length(alpha_grid) > 1) {
    if (!is.null(R_cv)) {
      alpha_used <- .cv_alpha_return(R_cv, k, alpha_grid,
                                     nlambda, lambda_min_ratio,
                                     nadd, nnested,
                                     standardize, epsilon, stabilize_sigma,
                                     compute_weights, normalize_weights, use_refit,
                                     n_folds)
    }
  }

  if (k == 0) {
    return(list(selection = integer(), weights = numeric(N), sr = 0, status = "LASSO_EMPTY", lambda = NA_real_, alpha = alpha))
  }
  if (k >= N) {
    sigma_s <- prep_covariance_cpp(sigma, epsilon, stabilize_sigma)
    sr_full <- compute_mve_sr_cpp(mu, sigma_s, selection = integer(), epsilon = epsilon,
                                  stabilize_sigma = FALSE, do_checks = FALSE)
    w_full <- if (compute_weights) {
      compute_mve_weights_cpp(mu, sigma_s, selection = integer(), normalize_w = normalize_weights,
                              epsilon = epsilon, stabilize_sigma = FALSE, do_checks = FALSE)
    } else {
      numeric(N)
    }
    return(list(selection = seq_len(N), weights = w_full, sr = sr_full, status = "LASSO_FULL", lambda = NA_real_, alpha = alpha))
  }

  sigma_s <- prep_covariance_cpp(sigma, epsilon, stabilize_sigma)
  design <- .design_from_moments(mu, sigma_s, n_obs)
  X <- design$X
  y <- design$y

  sel <- .select_lasso_support(X, y, k,
                               nlambda = nlambda,
                               lambda_min_ratio = lambda_min_ratio,
                               lambda_override = lambda,
                               alpha = alpha_used,
                               nadd = nadd,
                               nnested = nnested,
                               standardize = standardize)

  beta_sel <- sel$beta
  support <- sel$support
  lambda_sel <- sel$lambda
  status <- sel$status

  weights <- numeric(N)
  if (compute_weights) {
    if (use_refit && length(support) > 0) {
      weights <- compute_mve_weights_cpp(mu, sigma_s, selection = as.integer(support - 1),
                                         normalize_w = normalize_weights,
                                         epsilon = epsilon, stabilize_sigma = FALSE, do_checks = FALSE)
    } else {
      weights[support] <- beta_sel[support]
      if (normalize_weights) {
        weights <- normalize_weights_cpp(weights, mode = "relative", tol = 1e-6, do_checks = FALSE)
      }
    }
  }

  sr <- if (length(weights) == 0) 0 else {
    compute_sr_cpp(weights, mu, sigma_s, selection = integer(), epsilon = epsilon,
                   stabilize_sigma = FALSE, do_checks = FALSE)
  }

  list(
    selection = as.integer(support),
    weights = as.numeric(weights),
    sr = as.numeric(sr),
    status = status,
    lambda = lambda_sel,
    alpha = alpha_used
  )
}

#' LASSO relaxation heuristic for sparse MVE selection (glmnet backend, return-based)
#'
#' Pragmatic LASSO-based selector: builds a synthetic design from moments
#' \eqn{(mu, Sigma, T)}, fits a glmnet path, and chooses the support whose
#' cardinality is closest to `k` (preferring \eqn{\le k}). Optionally refits
#' exact MVE weights on the chosen support. Uses the same stabilized covariance
#' as the MIQP/exhaustive routines for consistency.
#'
#' @param R Numeric matrix of returns (rows = observations, cols = assets).
#' @param k Target cardinality.
#' @param nlambda Number of lambdas for glmnet path (ignored if `lambda` given).
#' @param lambda_min_ratio Smallest lambda as a fraction of lambda_max.
#' @param lambda Optional descending lambda grid.
#' @param alpha Elastic-net mixing parameter (1 = LASSO).
#' @param standardize Logical; whether glmnet standardizes the synthetic design.
#' @param epsilon Ridge for covariance stabilization.
#' @param stabilize_sigma Logical; symmetrize/ridge-stabilize sigma.
#' @param compute_weights Logical; return weights (else zeros).
#' @param normalize_weights Logical; post-scale weights via
#'   `normalize_weights_cpp(mode="relative")`.
#' @param use_refit Logical; refit exact MVE on selected support.
#' @param do_checks Logical; validate inputs.
#'
#' @return List with `selection` (1-based), `weights`, `sr`, `status`, `lambda`, `alpha`.
#' @export
mve_lasso_search_from_returns <- function(
    R,
    k,
    nlambda = 100L,
    lambda_min_ratio = 1e-3,
    lambda = NULL,
    alpha = 1,
    n_folds = 5L,
    nadd = 80L,
    nnested = 2L,
    standardize = FALSE,
    epsilon = eps_ridge_cpp(),
    stabilize_sigma = TRUE,
    compute_weights = TRUE,
    normalize_weights = FALSE,
    use_refit = FALSE,
    do_checks = FALSE
) {
  if (do_checks) {
    if (!is.matrix(R) || !is.numeric(R)) stop("R must be a numeric matrix.")
    if (nrow(R) < 2L) stop("R must have at least 2 rows.")
    if (ncol(R) < 1L) stop("R must have at least 1 column.")
    if (!is.numeric(k) || k < 0 || k > ncol(R)) stop("k must be between 0 and ncol(R).")
    if (!is.null(lambda)) {
      if (any(!is.finite(lambda)) || any(lambda <= 0)) stop("lambda must be positive/finite.")
    }
    if (!is.finite(epsilon)) stop("epsilon must be finite.")
  }

  Tobs <- nrow(R)
  mu <- colMeans(R)
  sigma <- stats::cov(R)

  mve_lasso_search(mu = mu,
                   sigma = sigma,
                   n_obs = Tobs,
                   k = k,
                   nlambda = nlambda,
                   lambda_min_ratio = lambda_min_ratio,
                   lambda = lambda,
                   alpha = alpha,
                   n_folds = n_folds,
                   nadd = nadd,
                   nnested = nnested,
                   standardize = standardize,
                   epsilon = epsilon,
                   stabilize_sigma = stabilize_sigma,
                   compute_weights = compute_weights,
                   normalize_weights = normalize_weights,
                   use_refit = use_refit,
                   do_checks = FALSE,
                   R_cv = R)
}
