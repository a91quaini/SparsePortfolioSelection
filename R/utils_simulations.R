#' FF3-style parameter builder (fixed constants)
#'
#' Returns the fixed factor model parameters used in the MATLAB simulations.
#' Values are monthly and match the hard-coded constants in the MATLAB code.
#'
#' @param factor_premium Length-3 numeric vector of factor premia.
#' @param factor_vol Length-3 numeric vector of factor volatilities.
#' @param alpha_mean Scalar mean for asset-specific intercepts.
#' @param alpha_vol Scalar std dev for asset-specific intercepts.
#' @param beta_mean Length-3 numeric vector of mean betas.
#' @param beta_vol Scalar or length-3 std dev for betas.
#' @param eps_vol Scalar idiosyncratic volatility.
#'
#' @return List of parameters.
#' @export
calibrate_ff3 <- function(
    n_assets,
    factor_premium = c(0.08, 0.03, 0.04) / 12,
    factor_vol = c(0.16, 0.14, 0.14) / sqrt(12),
    alpha_mean = 0.02 / 12,
    alpha_vol = 0.01 / sqrt(12),
    beta_mean = c(1.0, 0.0, 0.0) / 12,
    beta_vol = 0.3 / sqrt(12),
    eps_vol = 0.20 / sqrt(12),
    do_checks = TRUE
) {
  if (do_checks) {
    if (!is.numeric(n_assets) || n_assets < 1) stop("n_assets must be positive.")
  }
  K <- length(factor_premium)
  beta_vol_vec <- if (length(beta_vol) == 1L) rep(beta_vol, K) else beta_vol
  alpha_true <- alpha_mean + alpha_vol * stats::rnorm(n_assets)
  beta_true <- matrix(0, nrow = n_assets, ncol = K)
  for (k in seq_len(K)) {
    beta_true[, k] <- beta_mean[k] + beta_vol_vec[k] * stats::rnorm(n_assets)
  }
  mu_true <- alpha_true + beta_true %*% factor_premium
  Sigma_sys <- beta_true %*% diag(factor_vol^2) %*% t(beta_true)
  Sigma_true <- Sigma_sys + diag((eps_vol^2) * rep(1, n_assets))

  list(
    factor_premium = factor_premium,
    factor_vol = factor_vol,
    alpha_mean = alpha_mean,
    alpha_vol = alpha_vol,
    beta_mean = beta_mean,
    beta_vol = beta_vol,
    eps_vol = eps_vol,
    alpha_true = alpha_true,
    beta_true = beta_true,
    mu_true = as.numeric(mu_true),
    sigma_true = Sigma_true
  )
}

#' Simulate returns from fixed FF3-style parameters
#'
#' Generates asset-specific alphas/betas around fixed means, simulates factor
#' returns and idiosyncratic noise, and returns population moments plus a
#' T x N returns matrix (rows = observations).
#'
#' @param N Number of assets.
#' @param Tobs Number of observations.
#' @param params Parameter list from `calibrate_ff3()`.
#' @param do_checks Logical; validate inputs.
#'
#' @return List with `returns` and `factors`.
#' @export
simulate_ff3 <- function(Tobs, params, do_checks = FALSE) {
  if (do_checks) {
    if (!is.list(params)) stop("params must be a list from calibrate_ff3().")
    req <- c("factor_premium", "factor_vol", "alpha_mean", "alpha_vol", "beta_mean", "beta_vol", "eps_vol")
    req_full <- c(req, "alpha_true", "beta_true", "mu_true", "sigma_true")
    if (!all(req_full %in% names(params))) stop("params missing required fields (run calibrate_ff3()).")
    if (!is.numeric(Tobs) || Tobs < 1) stop("Tobs must be positive.")
  }

  N <- length(params$alpha_true)
  K <- length(params$factor_premium)
  alpha_true <- params$alpha_true
  beta_true <- params$beta_true

  # Simulate factor returns and idiosyncratic noise
  factors <- matrix(params$factor_premium, nrow = Tobs, ncol = K, byrow = TRUE) +
    matrix(params$factor_vol, nrow = Tobs, ncol = K, byrow = TRUE) * matrix(stats::rnorm(Tobs * K), nrow = Tobs, ncol = K)
  eps_noise <- params$eps_vol * matrix(stats::rnorm(Tobs * N), nrow = Tobs, ncol = N)

  returns <- matrix(0, nrow = Tobs, ncol = N)
  for (t in seq_len(Tobs)) {
    returns[t, ] <- alpha_true + factors[t, ] %*% t(beta_true) + eps_noise[t, ]
  }

  list(
    returns = returns,
    factors = factors
  )
}

#' Decompose card-k MVE Sharpe into estimation and selection terms
#'
#' Uses a sparse-search function on sample moments, then evaluates the resulting
#' selection/weights on population moments to split the Sharpe ratio into
#' estimation vs selection components.
#'
#' @param mu_pop Population mean vector.
#' @param sigma_pop Population covariance matrix.
#' @param mu_sample Sample mean vector.
#' @param sigma_sample Sample covariance matrix.
#' @param k Target cardinality.
#' @param mve_search_fn Function used to pick the subset on sample moments
#'   (e.g., `mve_exhaustive_search`, `mve_miqp_search`, `mve_lars_search`).
#' @param mve_search_fn_params List of extra args passed to `mve_search_fn`
#'   (excluding mu/sigma/k).
#' @param do_checks Logical; validate inputs.
#' @param return_selection Logical; also return selection/weights if TRUE.
#'
#' @return List with `mve_sr_cardk_est_term`, `mve_sr_cardk_sel_term`, and
#'   optionally `selection` (1-based) and `weights`.
#' @export
compute_mve_sr_decomposition <- function(mu_pop,
                                         sigma_pop,
                                         mu_sample,
                                         sigma_sample,
                                         k,
                                         mve_search_fn = mve_lars_search,
                                         mve_search_fn_params = list(),
                                         do_checks = FALSE,
                                         return_selection = FALSE) {
  n <- length(mu_pop)
  if (do_checks) {
    if (length(mu_sample) != n) stop("mu_sample length must equal mu_pop.")
    if (!is.matrix(sigma_pop) || any(dim(sigma_pop) != c(n, n))) stop("sigma_pop must be n x n.")
    if (!is.matrix(sigma_sample) || any(dim(sigma_sample) != c(n, n))) stop("sigma_sample must be n x n.")
    if (k < 1 || k > n) stop("k must be in 1..n.")
  }

  search_args <- c(list(mu_sample, sigma_sample, k), mve_search_fn_params)
  res <- tryCatch(do.call(mve_search_fn, search_args),
                  error = function(e) stop("mve_search_fn failed: ", conditionMessage(e)))

  sel <- res$selection
  w <- res$weights

  # Ensure full-length weights
  if (length(w) != n) {
    w_full <- numeric(n)
    if (length(sel) > 0 && length(w) == length(sel)) {
      w_full[sel] <- w
    } else {
      w_full <- rep_len(w, n)
    }
    w <- w_full
  }

  sel0 <- if (length(sel) == 0) integer() else as.integer(sel - 1L)

  est_term <- compute_sr(w, mu_pop, sigma_pop, selection = sel0, do_checks = FALSE)
  sel_term <- compute_mve_sr(mu_pop, sigma_pop, selection = sel0, do_checks = FALSE)

  out <- list(
    mve_sr_cardk_est_term = est_term,
    mve_sr_cardk_sel_term = sel_term
  )
  if (return_selection) {
    out$selection <- sel
    out$weights <- w
  }
  out
}

#' Compute population MVE SR across k-grid using multiple solvers
#'
#' For each k, compute the population SR using the max of:
#' 1) MIQP (boosted time_limit/mipgap, CPLEX-style bounds),
#' 2) LARS heuristic on the synthetic design.
#' Enforces non-decreasing SR in k (SR_k >= SR_{k-1}).
#'
#' @param mu_pop Population mean vector.
#' @param sigma_pop Population covariance matrix.
#' @param k_grid Vector of cardinalities.
#' @param n_obs Synthetic sample size for moment-based LASSO design.
#' @param miqp_time_limit Time limit for MIQP (seconds).
#' @param miqp_mipgap Relative MIP gap for MIQP.
#' @param verbose Logical; pass to MIQP.
#' @return Numeric vector of population SRs aligned with k_grid.
#' @export
compute_population_mve_sr <- function(mu_pop,
                                      sigma_pop,
                                      k_grid,
                                      n_obs,
                                      miqp_time_limit = 200,
                                      miqp_mipgap = 1e-5,
                                      verbose = FALSE) {
  pop_sr <- numeric(length(k_grid))
  prev_sr <- -Inf

  for (i in seq_along(k_grid)) {
    k <- k_grid[i]

    # 1) MIQP boosted
    sr_miqp <- tryCatch({
      res <- mve_miqp_search(
        mu_pop, sigma_pop, k,
        exactly_k = FALSE,
        gamma = 1.0,
        fmin = -0.25,
        fmax = 0.25,
        expand_rounds = 25L,
        expand_factor = 2.0,
        expand_tol = 1e-4,
        mipgap = miqp_mipgap,
        time_limit = miqp_time_limit,
        threads = 0,
        compute_weights = TRUE,
        normalize_weights = FALSE,
        use_refit = FALSE,
        verbose = verbose
      )
      as.numeric(res$sr)
    }, error = function(e) -Inf)

    # 2) LARS heuristic
    sr_lars <- tryCatch({
      res <- mve_lars_search(
        mu = mu_pop,
        sigma = sigma_pop,
        n_obs = n_obs,
        k = k,
        stabilize_sigma = FALSE,
        epsilon = eps_ridge_cpp(),
        tol_nnl = 1e-10,
        compute_weights = TRUE,
        normalize_weights = FALSE,
        use_refit = FALSE,
        do_checks = FALSE
      )
      as.numeric(res$sr)
    }, error = function(e) -Inf)

    best_sr <- max(sr_miqp, sr_lars, na.rm = TRUE)
    if (!is.finite(best_sr)) best_sr <- -Inf
    if (best_sr < prev_sr) best_sr <- prev_sr
    pop_sr[i] <- best_sr
    prev_sr <- best_sr
  }

  pop_sr
}

#' Sharpe-based systematic spanning index for a selected subset
#'
#' Computes a ratio-based "SpanSR" index for a selected subset A of assets:
#'   SpanSR(A) = theta^2(A) / theta^2(A augmented with systematic proxies),
#' where theta^2(X) = mu_X' Sigma_X^{-1} mu_X is the *population* max squared Sharpe
#' over the opportunity set X. By construction SpanSR(A) is in (0,1] (up to numerical
#' issues), and approaches 1 when A already spans the systematic proxies in the
#' mean--variance sense.
#'
#' Two choices of systematic proxies:
#' - sys = "factors": uses true FF factors from params (factor premia/vol + betas).
#' - sys = "pca": uses PCA portfolios g_t = M' r_t with user-supplied loadings M.
#'
#' @param mu_pop Population mean vector of the N assets.
#' @param sigma_pop Population covariance matrix of the N assets.
#' @param selection Integer vector of 1-based indices for the selected subset A.
#' @param sys Character; "factors" or "pca".
#' @param params List from calibrate_ff3(), required if sys = "factors".
#' @param pca_loadings N x m matrix of PCA portfolio weights (columns), required if sys = "pca".
#' @param m Number of systematic proxies to use (default: all columns / all factors).
#' @param ridge_epsilon Optional ridge added to the augmented covariance for stability.
#' @param do_checks Logical; validate inputs.
#' @param return_details Logical; include theta^2 components if TRUE.
#'
#' @return Either scalar SpanSR, or list with SpanSR and components if return_details=TRUE.
#' @export
compute_span_sr_index <- function(mu_pop,
                                  sigma_pop,
                                  selection,
                                  sys = c("factors", "pca"),
                                  params = NULL,
                                  pca_loadings = NULL,
                                  m = NULL,
                                  ridge_epsilon = 0.0,
                                  do_checks = FALSE,
                                  return_details = FALSE) {
  sys <- match.arg(sys)

  n <- length(mu_pop)
  A <- as.integer(unique(selection))
  A <- A[is.finite(A) & A >= 1L & A <= n]
  k <- length(A)

  if (do_checks) {
    if (!is.numeric(mu_pop) || length(mu_pop) != n) stop("mu_pop must be a numeric vector.")
    if (!is.matrix(sigma_pop) || any(dim(sigma_pop) != c(n, n))) stop("sigma_pop must be n x n.")
    if (k < 1L) stop("selection must contain at least one valid index in 1..N.")
    if (!is.numeric(ridge_epsilon) || ridge_epsilon < 0) stop("ridge_epsilon must be nonnegative.")
  }

  # subset moments
  mu_A <- as.numeric(mu_pop[A])
  Sigma_A <- sigma_pop[A, A, drop = FALSE]
  if (ridge_epsilon > 0) Sigma_A <- Sigma_A + diag(ridge_epsilon, k)

  # theta^2(A) = mu_A' Sigma_A^{-1} mu_A
  theta_A_sq <- tryCatch({
    as.numeric(crossprod(mu_A, solve(Sigma_A, mu_A)))
  }, error = function(e) NA_real_)

  # systematic proxy moments + cross-cov with selected assets
  if (sys == "factors") {
    if (is.null(params)) stop("params (from calibrate_ff3()) is required when sys = 'factors'.")
    if (is.null(m)) m <- length(params$factor_premium)
    m <- as.integer(m)

    # factor moments
    mu_g <- as.numeric(params$factor_premium[seq_len(m)])
    Sigma_g <- diag(as.numeric(params$factor_vol[seq_len(m)]^2), nrow = m, ncol = m)

    # cross-cov Cov(r_A, f) = B_A Sigma_f
    B_A <- params$beta_true[A, seq_len(m), drop = FALSE]
    Sigma_Ag <- B_A %*% Sigma_g  # k x m

  } else { # sys == "pca"
    if (is.null(pca_loadings)) stop("pca_loadings is required when sys = 'pca'.")
    if (is.null(m)) m <- ncol(pca_loadings)
    m <- as.integer(m)

    M <- pca_loadings[, seq_len(m), drop = FALSE]
    if (do_checks) {
      if (!is.matrix(M) || nrow(M) != n) stop("pca_loadings must be N x m.")
    }

    # g = M' r => mu_g = M' mu, Sigma_g = M' Sigma M, Cov(r_A, g) = Sigma[A, ] M
    mu_g <- as.numeric(crossprod(M, mu_pop))                        # m x 1
    Sigma_g <- crossprod(M, sigma_pop %*% M)                        # m x m
    Sigma_Ag <- sigma_pop[A, , drop = FALSE] %*% M                  # k x m
  }

  # augmented moments for [r_A ; g]
  mu_plus <- c(mu_A, mu_g)
  Sigma_plus <- rbind(
    cbind(Sigma_A, Sigma_Ag),
    cbind(t(Sigma_Ag), Sigma_g)
  )
  if (ridge_epsilon > 0) {
    Sigma_plus <- Sigma_plus + diag(ridge_epsilon, k + m)
  }

  theta_plus_sq <- tryCatch({
    as.numeric(crossprod(mu_plus, solve(Sigma_plus, mu_plus)))
  }, error = function(e) NA_real_)

  span_sr <- theta_A_sq / theta_plus_sq

  if (!return_details) return(span_sr)

  list(
    span_sr = span_sr,
    theta_A_sq = theta_A_sq,
    theta_plus_sq = theta_plus_sq,
    delta_sq = theta_plus_sq - theta_A_sq,
    k = k,
    m = m,
    sys = sys
  )
}

#' One Monte Carlo replication for FF3 simulation analysis
#'
#' Simulates returns under the FF3 calibration, computes sample moments, and returns
#' vectors of decomposition terms across k_grid. Optionally computes SpanSR indices.
#'
#' @param mc_id Integer replication id (unused except for reproducibility hooks).
#' @param n_obs Sample size T.
#' @param params Output of calibrate_ff3().
#' @param mu_pop Population mean of assets (length N).
#' @param sigma_pop Population covariance of assets (N x N).
#' @param k_grid Integer vector of k values.
#' @param mve_search_fn Function implementing the k-sparse search (e.g., mve_lars_search).
#' @param mve_search_fn_params Named list of params passed to mve_search_fn.
#' @param compute_span Logical; compute SpanSR if TRUE.
#' @param span_sys Character; "factors" or "pca".
#' @param span_m Integer; number of systematic proxies (m).
#' @param span_ridge_epsilon Ridge used in SpanSR computations.
#' @param do_checks Logical.
#'
#' @return List with est, sel, and optionally span vectors (each length(k_grid)).
#' @export
run_one_ff3_mc <- function(mc_id,
                           n_obs,
                           params,
                           mu_pop,
                           sigma_pop,
                           k_grid,
                           mve_search_fn,
                           mve_search_fn_params,
                           compute_span = FALSE,
                           span_sys = c("factors", "pca"),
                           span_m = 3L,
                           span_ridge_epsilon = 0.0,
                           do_checks = FALSE) {

  span_sys <- match.arg(span_sys)

  sim <- simulate_ff3(n_obs, params)
  R <- sim$returns
  mu_sample <- colMeans(R)
  sigma_sample <- stats::cov(R)

  # PCA loadings from in-sample moments if needed
  pca_loadings <- NULL
  if (compute_span && span_sys == "pca") {
    eig <- eigen(sigma_sample, symmetric = TRUE)
    m <- min(as.integer(span_m), ncol(eig$vectors))
    pca_loadings <- eig$vectors[, seq_len(m), drop = FALSE]
  }

  est_vec <- numeric(length(k_grid))
  sel_vec <- numeric(length(k_grid))
  span_vec <- if (compute_span) numeric(length(k_grid)) else NULL

  for (i in seq_along(k_grid)) {
    k <- k_grid[i]

    decomp <- compute_mve_sr_decomposition(
      mu_pop, sigma_pop,
      mu_sample, sigma_sample,
      k = k,
      mve_search_fn = mve_search_fn,
      mve_search_fn_params = mve_search_fn_params,
      do_checks = do_checks,
      return_selection = compute_span  # we need the support for SpanSR
    )

    est_vec[i] <- decomp$mve_sr_cardk_est_term
    sel_vec[i] <- decomp$mve_sr_cardk_sel_term

    if (compute_span) {
      span_vec[i] <- compute_span_sr_index(
        mu_pop, sigma_pop,
        selection = decomp$selection,
        sys = span_sys,
        params = if (span_sys == "factors") params else NULL,
        pca_loadings = if (span_sys == "pca") pca_loadings else NULL,
        m = as.integer(span_m),
        ridge_epsilon = span_ridge_epsilon,
        do_checks = do_checks,
        return_details = FALSE
      )
    }
  }

  if (!compute_span) return(list(est = est_vec, sel = sel_vec))
  list(est = est_vec, sel = sel_vec, span = span_vec)
}


#' Plot MVE SR decomposition results from MC analysis
#'
#' Expects the list saved by `inst/simulations/mc_analysis_ff3.R` containing
#' `population_sr`, `est_terms`, and `sel_terms`, along with the config.
#' If `span_terms` is present, also plots the average SpanSR index.
#' Additionally plots percentiles of the Sample(k) curve.
#'
#' @param results List with fields `population_sr`, `est_terms`, `sel_terms`,
#'   optionally `span_terms`, and `config` (including `k_grid`, `n_assets`, `n_obs`,
#'   `search_method`).
#' @param save Logical; if TRUE, save plots to `inst/simulations/figures`.
#'
#' @return A list with `sr_plot`, `risk_plot`, `loss_plot`, `sample_sr`,
#'   and optionally `span_plot`.
#' @export
plot_mve_sr_decomposition <- function(results, save = FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotting.")
  }

  cfg <- results$config
  k_grid <- cfg$k_grid
  pop_sr <- as.numeric(results$population_sr)
  opt_pop_sr <- results$opt_pop_sr
  est_mean <- colMeans(results$est_terms, na.rm = TRUE)
  sel_mean <- colMeans(results$sel_terms, na.rm = TRUE)

  # --- stylistic constants to match empirics ---
  base_size <- 14
  pt_size   <- 1.8
  line_w    <- 1.1

  # --- color palette (Okabe–Ito, colorblind-friendly) ---
  curve_cols <- c(
    "True"   = "#006400",  # dark green
    "Oracle" = "#0072B2",  # blue
    "Sample" = "#000000"   # black
  )
  if (!is.null(opt_pop_sr)) {
    curve_cols <- c("Optimal" = "#E41A1C", curve_cols) # red
  }

  # ---- compute sd band for Sample(k) from MC draws ----
  est_sd <- apply(results$est_terms, 2, stats::sd, na.rm = TRUE)
  df_sample_band <- data.frame(
    k = k_grid,
    lower = est_mean - est_sd,
    upper = est_mean + est_sd
  )

  # ---- build df_sr (unchanged logic) ----
  if (!is.null(opt_pop_sr)) {
    df_sr <- data.frame(
      k = rep(k_grid, 4),
      sr = c(rep(opt_pop_sr, length(k_grid)), pop_sr, sel_mean, est_mean),
      curve = rep(c("Optimal", "True", "Oracle", "Sample"), each = length(k_grid))
    )
    df_sr$curve <- factor(df_sr$curve, levels = c("Optimal", "True", "Oracle", "Sample"))
    df_sr_points <- df_sr[df_sr$curve != "Optimal", ]
  } else {
    df_sr <- data.frame(
      k = rep(k_grid, 3),
      sr = c(pop_sr, sel_mean, est_mean),
      curve = rep(c("True", "Oracle", "Sample"), each = length(k_grid))
    )
    df_sr$curve <- factor(df_sr$curve, levels = c("True", "Oracle", "Sample"))
    df_sr_points <- df_sr
  }

  # k^* from Sample(k) mean curve
  max_sample_k <- k_grid[which.max(est_mean)]
  max_sample   <- max(est_mean, na.rm = TRUE)

  # shapes: square True, triangle Oracle, circle Sample
  shape_vals <- c("True" = 15, "Oracle" = 17, "Sample" = 16)

  legend_breaks <- if (!is.null(opt_pop_sr)) {
    c("Optimal", "True", "Oracle", "Sample")
  } else {
    c("True", "Oracle", "Sample")
  }

  # override shapes in the single color legend (Optimal line only)
  shape_override <- rep(NA, length(legend_breaks))
  names(shape_override) <- legend_breaks
  shape_override[names(shape_vals)] <- shape_vals

  # --- SR plot base (so we can read y-min) ---
  sr_plot_base <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = df_sample_band,
      ggplot2::aes(x = k, ymin = lower, ymax = upper),
      fill = unname(curve_cols["Sample"]),
      alpha = 0.18,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_line(
      data = df_sr,
      ggplot2::aes(x = k, y = sr, color = curve),
      linewidth = line_w
    ) +
    ggplot2::geom_point(
      data = df_sr_points,
      ggplot2::aes(x = k, y = sr, color = curve, shape = curve),
      size = pt_size
    ) +
    ggplot2::scale_color_manual(
      values = curve_cols,
      labels = if (!is.null(opt_pop_sr)) {
        c(
          expression("True: " * theta^"*"),
          expression("True(k): " * theta[k]^"*"),
          expression("Oracle(k): " * theta[widehat(A)[k]]),
          expression("Sample(k): " * theta(hat(w)[k]))
        )
      } else {
        c(
          expression("True(k): " * theta[k]^"*"),
          expression("Oracle(k): " * theta[widehat(A)[k]]),
          expression("Sample(k): " * theta(hat(w)[k]))
        )
      },
      breaks = legend_breaks,
      name = NULL
    ) +
    ggplot2::scale_shape_manual(values = shape_vals, guide = "none") +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(
          shape = unname(shape_override),
          linetype = 1,
          linewidth = line_w,
          size = pt_size
        )
      )
    ) +
    ggplot2::labs(x = "Number of holdings k", y = "Sharpe ratio", title = NULL) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.title = ggplot2::element_blank()
    )

  pb <- ggplot2::ggplot_build(sr_plot_base)
  pp <- pb$layout$panel_params[[1]]
  y_min <- if (!is.null(pp$y.range)) pp$y.range[1] else pp$y$range[1]

  sr_plot <- sr_plot_base +
    ggplot2::annotate(
      "segment",
      x = max_sample_k, xend = max_sample_k,
      y = y_min, yend = max_sample,
      linetype = "dashed", linewidth = 0.6, color = "black"
    ) +
    ggplot2::annotate(
      "text", x = max_sample_k + 0.2, y = (y_min + max_sample) / 2,
      label = "k^\"*\"",
      hjust = -0.1, vjust = 0.5, parse = TRUE
    )

  # --- Risk decomposition plot: use colors (no gray) ---
  selection_risk  <- pop_sr - sel_mean
  allocation_risk <- sel_mean - est_mean
  df_risk_long <- data.frame(
    k = rep(k_grid, times = 2),
    component = factor(rep(c("Selection", "Allocation"), each = length(k_grid)),
                       levels = c("Selection", "Allocation")),
    value = c(selection_risk, allocation_risk)
  )

  risk_cols <- c(
    "Selection"  = unname(curve_cols["Oracle"]), # same as Oracle(k)
    "Allocation" = unname(curve_cols["Sample"])  # same as Sample(k)
  )

  risk_plot <- ggplot2::ggplot(df_risk_long, ggplot2::aes(x = k, y = value, fill = component)) +
    ggplot2::geom_area(position = ggplot2::position_stack(reverse = TRUE), alpha = 0.75) +
    ggplot2::scale_fill_manual(values = risk_cols, limits = c("Selection", "Allocation")) +
    ggplot2::labs(x = "Number of holdings k", y = "Estimation loss", fill = "Source of risk", title = NULL) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(legend.position = "bottom")

  # # --- Loss plot: use colors + point size 1.8 + unified legend (color+shape) ---
  # efficiency_loss <- if (!is.null(opt_pop_sr)) opt_pop_sr - pop_sr else rep(NA_real_, length(pop_sr))
  # estimation_loss <- pop_sr - est_mean
  # df_loss_long <- data.frame(
  #   k = rep(k_grid, times = 2),
  #   component = factor(rep(c("Efficiency", "Estimation"), each = length(k_grid)),
  #                      levels = c("Efficiency", "Estimation")),
  #   value = c(efficiency_loss, estimation_loss)
  # )
  #
  # loss_cols <- c(
  #   "Efficiency" = unname(curve_cols["True"]),   # same as True(k)
  #   "Estimation" = unname(curve_cols["Sample"])  # same as Sample(k)
  # )
  # loss_shapes <- c(
  #   "Efficiency" = unname(shape_vals["True"]),   # same shape as True(k)
  #   "Estimation" = unname(shape_vals["Sample"])  # same shape as Sample(k)
  # )
  #
  # loss_plot <- ggplot2::ggplot(df_loss_long, ggplot2::aes(x = k, y = value, color = component)) +
  #   ggplot2::geom_line(linewidth = line_w) +
  #   ggplot2::geom_point(ggplot2::aes(shape = component), size = pt_size) +
  #   ggplot2::scale_color_manual(values = loss_cols, name = NULL) +
  #   ggplot2::scale_shape_manual(values = loss_shapes, guide = "none") +
  #   ggplot2::guides(
  #     color = ggplot2::guide_legend(
  #       override.aes = list(
  #         shape = unname(loss_shapes[levels(df_loss_long$component)]),
  #         linetype = 1,
  #         linewidth = line_w,
  #         size = pt_size
  #       )
  #     )
  #   ) +
  #   ggplot2::labs(x = "Number of holdings k", y = "Sharpe ratio loss", title = NULL) +
  #   ggplot2::theme_minimal(base_size = base_size) +
  #   ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_blank())

  # --- Loss plot: annotate curves, no legend ---
  efficiency_loss <- if (!is.null(opt_pop_sr)) opt_pop_sr - pop_sr else rep(NA_real_, length(pop_sr))
  estimation_loss <- pop_sr - est_mean

  df_loss_long <- data.frame(
    k = rep(k_grid, times = 2),
    component = factor(rep(c("Efficiency", "Estimation"), each = length(k_grid)),
                       levels = c("Efficiency", "Estimation")),
    value = c(efficiency_loss, estimation_loss)
  )

  loss_cols <- c(
    "Efficiency" = unname(curve_cols["True"]),   # same as True(k)
    "Estimation" = unname(curve_cols["Sample"])  # same as Sample(k)
  )

  # pick reasonable k's to place labels near the curves
  k_eff_lab <- min(5L, max(k_grid))
  k_est_lab <- min(ceiling(0.75 * max(k_grid)), max(k_grid))

  y_eff_lab <- df_loss_long$value[df_loss_long$component == "Efficiency" & df_loss_long$k == k_eff_lab]
  y_est_lab <- df_loss_long$value[df_loss_long$component == "Estimation" & df_loss_long$k == k_est_lab]

  loss_plot <- ggplot2::ggplot(df_loss_long, ggplot2::aes(x = k, y = value, color = component)) +
    ggplot2::geom_line(linewidth = line_w) +
    ggplot2::scale_color_manual(values = loss_cols, guide = "none") +   # remove legend

    # in-plot curve labels (same color as curve)
    ggplot2::annotate(
      "text",
      x = k_eff_lab, y = y_eff_lab,
      label = "Efficiency loss",
      color = loss_cols["Efficiency"],
      hjust = 0, vjust = -0.6
    ) +
    ggplot2::annotate(
      "text",
      x = k_est_lab, y = y_est_lab,
      label = "Estimation loss",
      color = loss_cols["Estimation"],
      hjust = 1, vjust = -0.6
    ) +

    ggplot2::labs(x = "Number of holdings k", y = "Sharpe ratio loss", title = NULL) +
    ggplot2::theme_minimal(base_size = base_size)

  # (optionally keep points off; if you still want them, add geom_point back)


  # --- Span plot: colored ribbon + colored line/points + point size 1.8 ---
  span_plot <- NULL
  if (!is.null(results$span_terms)) {
    span_med <- apply(results$span_terms, 2, stats::quantile, probs = 0.50, na.rm = TRUE, names = FALSE)
    span_q05 <- apply(results$span_terms, 2, stats::quantile, probs = 0.05, na.rm = TRUE, names = FALSE)
    span_q95 <- apply(results$span_terms, 2, stats::quantile, probs = 0.95, na.rm = TRUE, names = FALSE)

    df_span <- data.frame(
      k = k_grid,
      med = as.numeric(span_med),
      q05 = as.numeric(span_q05),
      q95 = as.numeric(span_q95)
    )

    span_at_kstar <- df_span$med[df_span$k == max_sample_k]
    if (length(span_at_kstar) != 1L || !is.finite(span_at_kstar)) span_at_kstar <- NA_real_
    ymin_span <- suppressWarnings(min(c(df_span$q05, df_span$med), na.rm = TRUE))

    df_span_ribbon <- data.frame(
      k = df_span$k,
      ymin = df_span$q05,
      ymax = df_span$q95,
      band = "0.05–0.95"
    )

    span_col <- "#0072B2"

    span_plot <- ggplot2::ggplot() +
      ggplot2::geom_ribbon(
        data = df_span_ribbon,
        ggplot2::aes(x = k, ymin = ymin, ymax = ymax, fill = band),
        alpha = 0.18
      ) +
      ggplot2::geom_line(
        data = df_span,
        ggplot2::aes(x = k, y = med),
        linewidth = line_w,
        color = span_col
      ) +
      ggplot2::geom_point(
        data = df_span,
        ggplot2::aes(x = k, y = med),
        color = span_col,
        size = pt_size
      ) +
      ggplot2::scale_fill_manual(
        name = "Percentile",
        values = c("0.05–0.95" = span_col),
        breaks = c("0.05–0.95")
      ) +
      ggplot2::labs(x = "Number of holdings k", y = "Systematic risk spanning", title = NULL) +
      ggplot2::theme_minimal(base_size = base_size) +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::guides(color = "none")

    if (is.finite(max_sample_k) && is.finite(span_at_kstar) && is.finite(ymin_span)) {
      span_plot <- span_plot +
        ggplot2::annotate(
          "segment",
          x = max_sample_k, xend = max_sample_k,
          y = ymin_span, yend = span_at_kstar,
          linetype = "dashed", linewidth = 1.2, color = "black"
        ) +
        ggplot2::annotate(
          "text",
          x = max_sample_k + 0.2, y = (ymin_span + span_at_kstar) / 2,
          label = "k^\"*\"",
          hjust = -0.1, vjust = 0.5, parse = TRUE
        )
    }
  }

  if (save) {
    dir.create("inst/simulations/figures", recursive = TRUE, showWarnings = FALSE)
    base <- sprintf("ff3_%s_N%d_T%d_mc_analysis_ff3",
                    cfg$search_method, cfg$n_assets, cfg$n_obs)
    ggplot2::ggsave(file.path("inst/simulations/figures", paste0(base, "_sr.png")),
                    sr_plot, width = 8, height = 5, dpi = 150)
    ggplot2::ggsave(file.path("inst/simulations/figures", paste0(base, "_risk.png")),
                    risk_plot, width = 8, height = 5, dpi = 150)
    ggplot2::ggsave(file.path("inst/simulations/figures", paste0(base, "_loss.png")),
                    loss_plot, width = 8, height = 5, dpi = 150)
    if (!is.null(span_plot)) {
      ggplot2::ggsave(file.path("inst/simulations/figures", paste0(base, "_span.png")),
                      span_plot, width = 8, height = 5, dpi = 150)
    }
  }

  out <- list(sr_plot = sr_plot, risk_plot = risk_plot, loss_plot = loss_plot)
  if (!is.null(span_plot)) out$span_plot <- span_plot
  out
}



