# Out-of-sample analysis for ME×BE/ME (100) + 49 industry monthly portfolios.
# Loads the mebeme_ind panel via load_data(..., type = "mebeme_ind", frequency = "monthly").

## ---- thread control: must be at the very top ------------------------------
# Nn = 1L
Nn = 14L
Nn = min(Nn, parallel::detectCores(logical = TRUE) - 1L)
suppressPackageStartupMessages({
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(Nn)
    RhpcBLASctl::omp_set_num_threads(Nn)
    cat("BLAS threads set to:",
        RhpcBLASctl::blas_get_num_procs(), "\n")
    cat("OpenMP max threads set to:",
        RhpcBLASctl::omp_get_max_threads(), "\n")
  } else {
    warning("Package 'RhpcBLASctl' not available; falling back to env vars.")
    Sys.setenv(
      OMP_NUM_THREADS       = Nn,
      OPENBLAS_NUM_THREADS  = Nn,
      MKL_NUM_THREADS       = Nn,
      BLIS_NUM_THREADS      = Nn
    )
  }
})

# Cap Gurobi threads (default to BLAS cap; override via SPS_GUROBI_THREADS)
MIQP_THREADS <- as.integer(Sys.getenv("SPS_GUROBI_THREADS", Nn))
if (is.na(MIQP_THREADS) || MIQP_THREADS < 1L) MIQP_THREADS <- 1L

library(SparsePortfolioSelection)

# Configuration: 658 observations
PANEL_TYPE <- "mebeme"
MISSINGS <- "median"    # how to treat missing values
N_ASSETS <- 200         # subset of assets to use: total = 152 (100 + 49 + 3)
RNG_SEED <- 12345
W_IN_GRID <- c(360L)  # in-sample lengths (months)
W_OUT <- 1              # OOS block length (months)
OOS_TYPE <- "rolling"   # "rolling" or "expanding"
ADD_MKT <- TRUE         # append MKT-RF
ADD_FACTORS <- TRUE    # append FF3 (MKT, SMB, HML)
CHECK_K <- TRUE         # warn if solver returns sparsity different from k
K_TOL <- 1e-9           # tolerance for nonzero weights when checking sparsity
K_MIN <- 3
K_STEP <- 5
K_CAP <- N_ASSETS - 1
METHOD <- "lasso"        # "lasso" | "elnet" | "miqp"
REFIT <- FALSE
PARALLEL <- TRUE
COMPLETE_ANALYSIS <- FALSE  # if TRUE, run complete analysis (turnover/instability)

# Decide filename/label stems (append _refit if refit enabled)
refit_suffix <- if ((METHOD %in% c("lasso", "elnet") && REFIT) ||
                    (METHOD == "miqp" && REFIT)) "_refit" else ""
METHOD_LABEL <- paste0(METHOD, refit_suffix)
METHOD_STEM <- paste0(METHOD, refit_suffix)

OUT_DIR <- file.path("inst", "empirics", "results", "mebeme_ind_monthly")
FIG_DIR <- file.path("inst", "empirics", "figures", "mebeme_ind_monthly")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# Fix RNG before any shuffling inside load_data
if (!is.null(RNG_SEED)) set.seed(RNG_SEED)

ld <- load_data(type = PANEL_TYPE, missing = MISSINGS, path = "data", frequency = "monthly",
                add_mkt = ADD_MKT, add_factors = ADD_FACTORS)
R_all <- ld$returns
rf_vec <- if (is.null(ld$rf)) rep(0, nrow(R_all)) else ld$rf
T_full <- nrow(R_all); N_full <- ncol(R_all)

# Select a subset of assets for speed/reproducibility
N <- min(N_ASSETS, N_full)
R_all <- R_all[, 1:N, drop = FALSE]
if (!is.null(rf_vec)) rf_vec <- rf_vec[seq_len(nrow(R_all))]

k_grid <- NULL  # will be set after we know N

# Define parameter lists
alpha_grid = seq(0.30, 1.00, by = 0.05)
if (METHOD == "lasso") {
  alpha_grid = 1.00
}
lasso_params <- list(
  nlambda = 400L,
  lambda_min_ratio = 1e-4,
  alpha = alpha_grid,
  n_folds = 5L,
  nadd = 100L,
  nnested = 4L,
  standardize = FALSE,
  stabilize_sigma = FALSE,
  compute_weights = TRUE,
  normalize_weights = FALSE,
  use_refit = REFIT
)

miqp_params <- list(
  exactly_k = TRUE,
  m = 1L,
  gamma = 1.0,
  fmin = -0.25,
  fmax = 0.25,
  expand_rounds = 6L,
  expand_factor = 3.0,
  expand_tol = 1e-2,
  mipgap = 1e-4,
  time_limit = 100,
  threads = MIQP_THREADS,
  compute_weights = TRUE,
  normalize_weights = FALSE,
  use_refit = REFIT,
  verbose = FALSE,
  stabilize_sigma = TRUE
)

compute_weights_fn <- if (METHOD == "miqp") {
  function(Rin, k) {
    mu <- colMeans(Rin)
    sigma <- cov_fast(Rin)
    res <- do.call(mve_miqp_search, c(list(mu, sigma, k), miqp_params))
    wopt <- res$weights
    if (CHECK_K && !is.null(wopt)) {
      nz <- sum(abs(wopt) > K_TOL)
      if (nz != k) warning(sprintf("MIQP returned %d nonzero weights but target k=%d", nz, k))
    }
    list(weights = wopt, selection = res$selection, status = res$status)
  }
} else {
  function(Rin, k) {
    mu <- colMeans(Rin)
    sigma <- cov_fast(Rin)
    res <- do.call(mve_lasso_search, c(list(mu = mu, sigma = sigma, n_obs = nrow(Rin), k = k), lasso_params))
    wopt <- res$weights
    if (CHECK_K && !is.null(wopt)) {
      nz <- sum(abs(wopt) > K_TOL)
      if (nz != k) warning(sprintf("LASSO returned %d nonzero weights but target k=%d", nz, k))
    }
    list(weights = wopt, selection = res$selection, status = res$status)
  }
}

for (W_IN in W_IN_GRID) {
  R <- R_all
  Tobs <- nrow(R); N <- ncol(R)
  k_max <- min(K_CAP, N - 1)
  if (k_max < K_MIN) stop("k_max < K_MIN; reduce K_MIN or increase N.")
  k_grid <- seq.int(K_MIN, k_max, by = K_STEP)

  message(sprintf("Starting OOS run: T=%d, N=%d, W_IN=%d, W_OUT=%d, k ∈ [%d..%d]", Tobs, N, W_IN, W_OUT, K_MIN, k_max))

# Optional parallel run: set PARALLEL <- TRUE to enable
  if (COMPLETE_ANALYSIS) {
    if (PARALLEL) {
      n_cores <- Nn
      res <- run_complete_oos_evaluation_parallel(
        R = R,
        size_w_in = W_IN,
        size_w_out = W_OUT,
        k_grid = k_grid,
        oos_type = OOS_TYPE,
        compute_weights_fn = compute_weights_fn,
        compute_weights_fn_params = list(),
        rf = rf_vec,
        sharpe_fn = "median",
        n_cores = n_cores,
        return_details = TRUE
      )
    } else {
      res <- run_complete_oos_evaluation(
        R = R,
        size_w_in = W_IN,
        size_w_out = W_OUT,
        k_grid = k_grid,
        oos_type = OOS_TYPE,
        compute_weights_fn = compute_weights_fn,
        compute_weights_fn_params = list(),
        rf = rf_vec,
        sharpe_fn = "median",
        return_details = TRUE
      )
    }
    SR <- matrix(res$summary$oos_sr, ncol = 1)
    labels <- METHOD_LABEL
    less_than_k <- integer(length(k_grid))
  } else {
    if (PARALLEL) {
      n_cores <- Nn
      res <- run_oos_evaluation_parallel(
        R = R,
        size_w_in = W_IN,
        size_w_out = W_OUT,
        k_grid = k_grid,
        oos_type = OOS_TYPE,
        compute_weights_fn = compute_weights_fn,
        compute_weights_fn_params = list(),
        n_cores = n_cores,
        return_details = TRUE
      )
    } else {
      res <- run_oos_evaluation(
        R = R,
        size_w_in = W_IN,
        size_w_out = W_OUT,
        k_grid = k_grid,
        oos_type = OOS_TYPE,
        compute_weights_fn = compute_weights_fn,
        compute_weights_fn_params = list(),
        return_details = TRUE
      )
    }

    SR <- matrix(res$oos_by_k, ncol = 1)
    labels <- METHOD_LABEL
    less_than_k <- integer(length(k_grid))
    if (!is.null(res$selection)) {
      for (ik in seq_along(k_grid)) {
        k <- k_grid[ik]
        sel_mat <- res$selection[[ik]]
        if (!is.null(sel_mat)) {
          nz_counts <- apply(sel_mat, 1, function(x) sum(abs(x) > K_TOL))
          less_than_k[ik] <- sum(nz_counts < k)
        }
      }
    }
  }

  cat("\n=== Average OOS Sharpe by k ===\n")
  print_results(k_grid, SR, method_labels = labels, digits = 4)

  panel_tag <- if (PANEL_TYPE == "mebeme") "mebeme" else "mebeme_ind"
  factors_tag <- if (ADD_FACTORS) "ff3" else "nofactors"
  mkt_tag <- if (ADD_MKT) "mkt" else "nomkt"
  stem <- sprintf("oos_sr_%s_%s_monthly_%s_%s_N%d_Win%d_Wout%d",
                  METHOD_STEM, panel_tag, factors_tag, mkt_tag, N, W_IN, W_OUT)
  csv_path <- file.path(OUT_DIR, paste0(stem, ".csv"))
  plot_base <- file.path(FIG_DIR, stem)

  # augment results with less_than_k column
  res_table <- data.frame(k = k_grid, SharpeRatio = SR[, 1], less_than_k = less_than_k)
  write.csv(res_table, csv_path, row.names = FALSE)
  message("Saved results to: ", csv_path)

  plot_sr_empirics(k_grid, SR, save_path = plot_base)
  message("Saved figure to: ", plot_base, ".png")

  if (COMPLETE_ANALYSIS) {
    # Plot turnover / instabilities if available
    if (!is.null(res$summary$median_turnover)) {
      plot_turnover_empirics(k_grid, res$summary$median_turnover, method_labels = labels,
                             save_path = paste0(plot_base, "_turnover"))
    }
    if (!is.null(res$summary$median_weight_instability_L1)) {
      plot_weight_instability_empirics(
        k_grid,
        res$summary$median_weight_instability_L1,
        res$summary$median_weight_instability_L2,
        method_labels = labels,
        save_path = paste0(plot_base, "_weight_instability")
      )
    }
    if (!is.null(res$summary$median_selection_instability)) {
      plot_selection_instability_empirics(
        k_grid,
        res$summary$median_selection_instability,
        method_labels = labels,
        save_path = paste0(plot_base, "_selection_instability")
      )
    }
  }
}
