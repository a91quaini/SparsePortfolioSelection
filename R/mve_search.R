#' Exhaustive search over k-subsets maximizing Sharpe ratio
#'
#' Calls the C++ exhaustive (or sampled) search to find the k-asset subset with
#' highest in-sample Sharpe ratio. A ridge term can be applied once to the
#' covariance matrix before scoring subsets.
#'
#' @param mu Numeric mean vector (length \eqn{N}).
#' @param sigma Numeric covariance matrix (\eqn{N \times N}).
#' @param k Target cardinality (1..N).
#' @param ridge_epsilon Ridge factor for covariance stabilization (0 disables).
#' @param enumerate_all Logical; `TRUE` to enumerate all \eqn{N \choose k};
#'   `FALSE` to sample supports.
#' @param max_samples Number of sampled supports when `enumerate_all = FALSE`.
#' @param dedup_samples Logical; ensure sampled supports are distinct.
#' @param compute_weights Logical; if `TRUE`, return MVE weights for the best subset.
#' @param normalize_weights Logical; if `TRUE`, normalize weights.
#' @param normalization_type Integer; 1 = L1 (default), 0 = sum-to-one.
#' @param do_checks Logical; if `TRUE`, validate inputs.
#'
#' @return List with `selection` (1-based indices), `weights`, `sr`, `status`.
#' @export
mve_exhaustive_search <- function(mu, sigma, k,
                                  ridge_epsilon = 0.0,
                                  enumerate_all = TRUE,
                                  max_samples = 0L,
                                  dedup_samples = FALSE,
                                  compute_weights = FALSE,
                                  normalize_weights = FALSE,
                                  normalization_type = 1L,
                                  do_checks = FALSE) {
  if (!enumerate_all && max_samples <= 0) {
    stop("When enumerate_all=FALSE, max_samples must be > 0.")
  }
  res <- mve_exhaustive_search_cpp(mu, sigma, k, ridge_epsilon,
                                   enumerate_all, max_samples, dedup_samples,
                                   compute_weights, normalize_weights,
                                   normalization_type, do_checks)
  res$selection <- as.integer(res$selection) + 1L
  res
}

#' MIQP heuristic search for sparse selection (Gurobi backend)
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
#' @param normalize_weights Logical; if `TRUE`, add budget constraint and
#'   post-scale refit weights.
#' @param use_refit Logical; if `TRUE`, refit exact MVE on the selected support.
#' @param verbose Logical; emit Gurobi output.
#' @param ridge_epsilon Ridge for covariance stabilization.
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
                            normalize_weights = FALSE,
                            use_refit = FALSE,
                            verbose = FALSE,
                            ridge_epsilon = 0.0,
                            do_checks = FALSE) {
  fmin_vec <- if (length(fmin) == 1L) rep_len(fmin, length(mu)) else fmin
  fmax_vec <- if (length(fmax) == 1L) rep_len(fmax, length(mu)) else fmax
  res <- mve_miqp_search_cpp(mu, sigma, k, fmin_vec, fmax_vec,
                             x_start, v_start,
                             m, gamma, exactly_k,
                             expand_rounds, expand_factor, expand_tol,
                             mipgap, time_limit, threads,
                             ridge_epsilon,
                             normalize_weights,
                             use_refit, verbose,
                             do_checks)
  # Convert 0-based selection from C++ to 1-based R
  res$selection <- as.integer(res$selection) + 1L
  res
}

# ---- small helpers (top-level, not nested) --------------------------------

.sps_pick_path_idx <- function(nnz, k) {
  idx_exact <- which(nnz == k)
  if (length(idx_exact) > 0L) return(list(idx = idx_exact[1L], status = "LARS_EXACT_K"))

  idx_le <- which(nnz <= k)
  if (length(idx_le) == 0L) return(list(idx = 1L, status = "LARS_BELOW_K_FALLBACK"))

  max_nnz <- max(nnz[idx_le])
  idx_tie <- idx_le[nnz[idx_le] == max_nnz]
  idx <- max(idx_tie)
  status <- if (max_nnz < k) "LARS_BELOW_K" else "LARS_OK"
  list(idx = idx, status = status)
}

.sps_maybe_normalize <- function(w, normalize_weights, normalization_type) {
  if (!isTRUE(normalize_weights)) return(w)
  # If your normalize_weights_cpp signature is (w, epsilon, type), this works.
  # If it is older (w, tol), adapt here once, not everywhere else.
  normalize_weights_cpp(w, epsilon = 1e-8, type = as.integer(normalization_type))
}

#' LARS heuristic for sparse MVE selection (covariance-based)
#'
#' Uses the `lars` package to compute a LASSO solution path on the synthetic
#' regression design implied by (mu, sigma, n_obs) under the mean--variance criterion.
#'
#' Selection rule: (i) first solution with exactly k nonzeros (within tol_nnl),
#' else (ii) solution with the largest number of nonzeros not exceeding k (ties
#' broken in favor of the later point).
#'
#' @param mu Numeric mean vector.
#' @param sigma Numeric covariance matrix.
#' @param n_obs Number of observations underlying the moments.
#' @param k Target sparsity level.
#' @param ridge_epsilon Ridge factor for stabilization.
#' @param tol_nnl Threshold for treating coefficients as non-zero.
#' @param normalize_weights Logical; if TRUE, post-normalize weights via
#'   normalize_weights_cpp(). Default normalization_type is L1.
#' @param normalization_type Integer; 1 = L1 (default), 0 = sum-to-one.
#' @param use_refit Logical; if TRUE, refit exact MVE weights on selected support.
#' @param do_checks Logical; if TRUE, validate inputs.
#'
#' @return List with selection (1-based indices), weights, sr, status.
#' @export
mve_lars_search <- function(
    mu,
    sigma,
    n_obs,
    k,                     # <- scalar or vector
    ridge_epsilon = 0.0,
    tol_nnl = 1e-10,
    normalize_weights = FALSE,
    normalization_type = 1L,
    use_refit = FALSE,
    compute_sr = FALSE,
    do_checks = FALSE
) {

  mu <- as.numeric(mu)
  sigma <- as.matrix(sigma)
  N <- length(mu)

  k_vec <- as.integer(k)
  if (!length(k_vec)) stop("k must be non-empty.")
  if (any(!is.finite(k_vec))) stop("k must be finite.")
  k_vec[k_vec > N] <- N

  if (do_checks) {
    if (N == 0) stop("mu must be non-empty.")
    if (!is.matrix(sigma) || nrow(sigma) != ncol(sigma) || nrow(sigma) != N) {
      stop("sigma must be square and match length(mu).")
    }
    if (!is.numeric(n_obs) || length(n_obs) != 1 || !is.finite(n_obs) || n_obs <= 0) {
      stop("n_obs must be a positive finite scalar.")
    }
    if (any(k_vec < 0L)) stop("k must be >= 0.")
    if (!is.finite(ridge_epsilon) || ridge_epsilon < 0) stop("ridge_epsilon must be finite and nonnegative.")
    if (!is.finite(tol_nnl) || tol_nnl <= 0) stop("tol_nnl must be positive and finite.")
    if (!is.numeric(normalization_type) || length(normalization_type) != 1) stop("normalization_type must be scalar.")
  }

  normalization_type <- as.integer(normalization_type)

  # handle k <= 0 (return zero weights for those entries)
  is_zero_k <- (k_vec <= 0L)

  # FULL support shortcut only if ALL requests are full support and scalar-ish
  # (in practice your k_grid won’t include N, so this is mostly irrelevant)
  # We keep it correct anyway:
  if (all(k_vec >= N)) {
    w_full <- compute_mve_weights_cpp(
      mu, sigma, selection = integer(),
      ridge_epsilon = ridge_epsilon,
      normalize_weights = normalize_weights,
      normalization_type = normalization_type,
      do_checks = FALSE
    )
    sr_full <- if (isTRUE(compute_sr)) {
      compute_sr_cpp(w_full, mu, sigma, selection = integer(),
                     ridge_epsilon = ridge_epsilon, do_checks = FALSE)
    } else NA_real_

    if (length(k_vec) == 1L) {
      return(list(
        selection = seq_len(N),
        weights   = as.numeric(w_full),
        sr        = as.numeric(sr_full),
        status    = "LARS_FULL"
      ))
    }

    W <- matrix(rep(as.numeric(w_full), each = length(k_vec)), nrow = length(k_vec), byrow = TRUE)
    Sel <- replicate(length(k_vec), seq_len(N), simplify = FALSE)
    Status <- rep("LARS_FULL", length(k_vec))
    SR <- rep(as.numeric(sr_full), length(k_vec))
    return(list(k = k_vec, selection = Sel, weights = W, sr = SR, status = Status))
  }

  # Build design from covariance moments (C++ does stabilization via ridge_epsilon)
  design <- tryCatch(
    design_from_moments_cpp(mu, sigma, n_obs, ridge_epsilon),
    error = function(e) NULL
  )
  if (is.null(design) || is.null(design$X) || is.null(design$y)) {
    if (length(k_vec) == 1L) {
      return(list(selection = integer(), weights = numeric(N), sr = NA_real_, status = "LARS_DESIGN_FAIL"))
    }
    return(list(
      k = k_vec,
      selection = replicate(length(k_vec), integer(), simplify = FALSE),
      weights = matrix(0, nrow = length(k_vec), ncol = N),
      sr = rep(NA_real_, length(k_vec)),
      status = rep("LARS_DESIGN_FAIL", length(k_vec))
    ))
  }

  fit <- tryCatch(
    lars::lars(x = design$X, y = design$y, type = "lasso", intercept = FALSE, normalize = FALSE),
    error = function(e) NULL
  )
  if (is.null(fit) || is.null(fit$beta)) {
    if (length(k_vec) == 1L) {
      return(list(selection = integer(), weights = numeric(N), sr = NA_real_, status = "LARS_FIT_FAIL"))
    }
    return(list(
      k = k_vec,
      selection = replicate(length(k_vec), integer(), simplify = FALSE),
      weights = matrix(0, nrow = length(k_vec), ncol = N),
      sr = rep(NA_real_, length(k_vec)),
      status = rep("LARS_FIT_FAIL", length(k_vec))
    ))
  }

  B <- fit$beta
  if (is.null(dim(B))) B <- matrix(B, nrow = 1L)
  nnz <- rowSums(abs(B) > tol_nnl)

  Kreq <- length(k_vec)
  W_out <- matrix(0, nrow = Kreq, ncol = N)
  Sel_out <- vector("list", Kreq)
  Status_out <- character(Kreq)
  SR_out <- rep(NA_real_, Kreq)

  for (ii in seq_len(Kreq)) {
    kk <- k_vec[ii]
    if (kk <= 0L) {
      Sel_out[[ii]] <- integer()
      Status_out[ii] <- "LARS_EMPTY"
      if (isTRUE(compute_sr)) SR_out[ii] <- 0
      next
    }

    pick <- .sps_pick_path_idx(nnz, kk)
    idx <- pick$idx
    status <- pick$status

    w_path <- as.numeric(B[idx, ])
    sel <- which(abs(w_path) > tol_nnl)

    if (!length(sel)) {
      Sel_out[[ii]] <- integer()
      Status_out[ii] <- "LARS_ZERO"
      if (isTRUE(compute_sr)) SR_out[ii] <- 0
      next
    }

    if (isTRUE(use_refit)) {
      sel0 <- as.integer(sel - 1L)
      w_refit <- compute_mve_weights_cpp(
        mu, sigma, selection = sel0,
        ridge_epsilon = ridge_epsilon,
        normalize_weights = normalize_weights,
        normalization_type = normalization_type,
        do_checks = FALSE
      )
      w_final <- as.numeric(w_refit)
    } else {
      w_final <- .sps_maybe_normalize(w_path, normalize_weights, normalization_type)
    }

    W_out[ii, ] <- w_final
    Sel_out[[ii]] <- as.integer(sel)
    Status_out[ii] <- status

    if (isTRUE(compute_sr)) {
      SR_out[ii] <- compute_sr_cpp(
        w_final, mu, sigma, selection = integer(),
        ridge_epsilon = ridge_epsilon, do_checks = FALSE
      )
      if (!is.finite(SR_out[ii])) SR_out[ii] <- 0
    }
  }

  if (length(k_vec) == 1L) {
    return(list(
      selection = Sel_out[[1L]],
      weights   = as.numeric(W_out[1L, ]),
      sr        = as.numeric(SR_out[1L]),
      status    = Status_out[1L]
    ))
  }

  list(k = k_vec, selection = Sel_out, weights = W_out, sr = SR_out, status = Status_out)
}

