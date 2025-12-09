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
#' @return List with `mu_true`, `Sigma_true`, `returns`, `alpha_true`, `beta_true`.
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

  returns
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
#'   (e.g., `mve_exhaustive_search`, `mve_miqp_search`, `mve_lasso_search_return_based`).
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
                                         mve_search_fn = mve_exhaustive_search,
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
#' 2) LASSO with boosted path density,
#' 3) LASSO with alpha CV (if alpha grid provided and optionally R_cv).
#' Enforces non-decreasing SR in k (SR_k >= SR_{k-1}).
#'
#' @param mu_pop Population mean vector.
#' @param sigma_pop Population covariance matrix.
#' @param k_grid Vector of cardinalities.
#' @param n_obs Synthetic sample size for moment-based LASSO design.
#' @param alpha_grid Alpha grid for CV (set length>1 to enable).
#' @param lasso_nlambda Number of lambdas for boosted LASSO path.
#' @param lasso_lambda_min_ratio Smallest lambda fraction for LASSO.
#' @param miqp_time_limit Time limit for MIQP (seconds).
#' @param miqp_mipgap Relative MIP gap for MIQP.
#' @param R_cv Optional returns matrix for alpha CV (moments API).
#' @param verbose Logical; pass to MIQP.
#' @return Numeric vector of population SRs aligned with k_grid.
#' @export
compute_population_mve_sr <- function(mu_pop,
                                      sigma_pop,
                                      k_grid,
                                      n_obs,
                                      lasso_nlambda = 300L,
                                      lasso_lambda_min_ratio = 1e-4,
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
        exactly_k = TRUE,
        gamma = 1.0,
        fmin = -0.25,
        fmax = 0.25,
        expand_rounds = 20L,
        expand_factor = 3.0,
        expand_tol = 1e-2,
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

    # 2) LASSO boosted (alpha fixed)
    sr_lasso <- tryCatch({
      res <- mve_lasso_search(
        mu = mu_pop,
        sigma = sigma_pop,
        n_obs = n_obs,
        k = k,
        nlambda = lasso_nlambda,
        lambda_min_ratio = lasso_lambda_min_ratio,
        alpha = 1,
        nadd = 50L,
        nnested = 2L,
        standardize = FALSE,
        compute_weights = TRUE,
        normalize_weights = FALSE,
        use_refit = FALSE
      )
      as.numeric(res$sr)
    }, error = function(e) -Inf)

    best_sr <- max(sr_miqp, sr_lasso, na.rm = TRUE)
    if (!is.finite(best_sr)) best_sr <- -Inf
    if (best_sr < prev_sr) best_sr <- prev_sr
    pop_sr[i] <- best_sr
    prev_sr <- best_sr
  }

  pop_sr
}

#' Plot MVE SR decomposition results from MC analysis
#'
#' Expects the list saved by `inst/simulations/mc_analysis_ff3.R` containing
#' `population_sr`, `est_terms`, and `sel_terms`, along with the config.
#'
#' @param results List with fields `population_sr`, `est_terms`, `sel_terms`,
#'   and `config` (including `k_grid`, `n_assets`, `n_obs`, `search_method`).
#' @param save Logical; if TRUE, save plots to `inst/simulations/figures`.
#'
#' @return A list with `sr_plot` and `risk_plot` ggplot objects.
#' @export
plot_mve_sr_decomposition <- function(results, save = FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotting.")
  }

  cfg <- results$config
  k_grid <- cfg$k_grid
  pop_sr <- as.numeric(results$population_sr)
  est_mean <- colMeans(results$est_terms, na.rm = TRUE)
  sel_mean <- colMeans(results$sel_terms, na.rm = TRUE)

  df_sr <- data.frame(
    k = rep(k_grid, 3),
    sr = c(pop_sr, sel_mean, est_mean),
    curve = rep(c("True", "Oracle", "Sample"), each = length(k_grid))
  )
  df_sr$curve <- factor(df_sr$curve, levels = c("True", "Oracle", "Sample"))

  sr_plot <- ggplot2::ggplot(df_sr, ggplot2::aes(x = k, y = sr, color = curve)) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::geom_point() +
    ggplot2::labs(x = "k (allocation cardinality)", y = "Population Sharpe ratio", color = "MVE weights", title = NULL) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      legend.title = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 12),
      axis.title = ggplot2::element_text(size = 14),
      axis.text = ggplot2::element_text(size = 12)
    )

  allocation_risk <- pop_sr - sel_mean
  selection_risk <- sel_mean - est_mean
  df_risk <- data.frame(
    k = k_grid,
    Selection = selection_risk,
    Allocation = allocation_risk
  )
  df_risk_long <- tidyr::pivot_longer(df_risk, cols = c("Selection", "Allocation"),
                                      names_to = "component", values_to = "value")
  df_risk_long$component <- factor(df_risk_long$component, levels = c("Selection", "Allocation"))

  risk_plot <- ggplot2::ggplot(df_risk_long, ggplot2::aes(x = k, y = value, fill = component)) +
    ggplot2::geom_area(position = "stack", alpha = 0.8) +
    ggplot2::scale_x_continuous() +
    ggplot2::labs(x = "k (allocation cardinality)", y = "Population Sharpe ratio gap", fill = "Source of risk", title = NULL) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      legend.title = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 12),
      axis.title = ggplot2::element_text(size = 14),
      axis.text = ggplot2::element_text(size = 12)
    )

  if (save) {
    dir.create("inst/simulations/figures", recursive = TRUE, showWarnings = FALSE)
    base <- sprintf("ff3_%s_N%d_T%d_mc_analysis_ff3",
                    cfg$search_method, cfg$n_assets, cfg$n_obs)
    ggplot2::ggsave(file.path("inst/simulations/figures", paste0(base, "_sr.png")), sr_plot, width = 7, height = 5, dpi = 150)
    ggplot2::ggsave(file.path("inst/simulations/figures", paste0(base, "_risk.png")), risk_plot, width = 7, height = 5, dpi = 150)
  }

  list(sr_plot = sr_plot, risk_plot = risk_plot)
}
