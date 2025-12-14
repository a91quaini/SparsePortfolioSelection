# Out-of-sample analysis for ME×BE/ME (100) + 49 industry monthly portfolios.
# Loads the mebeme_ind panel via load_data(..., type = "mebeme_ind", frequency = "monthly").

## ---- thread control: must be at the very top ------------------------------
# Nn = 1L
Nn = 6L
Nn = min(Nn, parallel::detectCores(logical = TRUE) - 1L)
# suppressPackageStartupMessages({
#   if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
#     RhpcBLASctl::blas_set_num_threads(Nn)
#     RhpcBLASctl::omp_set_num_threads(Nn)
#     cat("BLAS threads set to:",
#         RhpcBLASctl::blas_get_num_procs(), "\n")
#     cat("OpenMP max threads set to:",
#         RhpcBLASctl::omp_get_max_threads(), "\n")
#   } else {
#     warning("Package 'RhpcBLASctl' not available; falling back to env vars.")
#     Sys.setenv(
#       OMP_NUM_THREADS       = Nn,
#       OPENBLAS_NUM_THREADS  = Nn,
#       MKL_NUM_THREADS       = Nn,
#       BLIS_NUM_THREADS      = Nn
#     )
#   }
# })

# Cap Gurobi threads (default to BLAS cap; override via SPS_GUROBI_THREADS)
MIQP_THREADS <- as.integer(Sys.getenv("SPS_GUROBI_THREADS", Nn))
if (is.na(MIQP_THREADS) || MIQP_THREADS < 1L) MIQP_THREADS <- 1L

library(SparsePortfolioSelection)

# Configuration
PANEL_TYPE <- "mebeme_ind"
MISSINGS <- "median"    # how to treat missing values
N_ASSETS <- 150         # subset of assets to use: total = 149 (100 + 49)
RNG_SEED <- 12345
W_IN_GRID <- c(120L, 240L, 360L, 480L, 600L)  # in-sample lengths (months)
W_OUT <- 1              # OOS block length (months)
OOS_TYPE <- "rolling"   # "rolling" or "expanding"
K_MIN <- 3
K_STEP <- 2
K_CAP <- N_ASSETS - 1
METHOD <- "lasso"        # "lasso" | "elnet" | "miqp"
REFIT <- FALSE
PARALLEL <- TRUE

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

R_all <- load_data(type = PANEL_TYPE, missing = MISSINGS, path = "data", frequency = "monthly")
T_full <- nrow(R_all); N_full <- ncol(R_all)

# Select a subset of assets for speed/reproducibility
N <- min(N_ASSETS, N_full)
asset_idx <- sort(sample.int(N_full, N))
R_all <- R_all[, asset_idx, drop = FALSE]

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
  stabilize_sigma = TRUE,
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
  expand_rounds = 10L,
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
    list(weights = res$weights, selection = res$selection, status = res$status)
  }
} else {
  function(Rin, k) {
    mu <- colMeans(Rin)
    sigma <- cov_fast(Rin)
    res <- do.call(mve_lasso_search, c(list(mu = mu, sigma = sigma, n_obs = nrow(Rin), k = k), lasso_params))
    list(weights = res$weights, selection = res$selection, status = res$status)
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
  if (PARALLEL) {
    n_cores <- Nn
    sr_vec <- run_oos_evaluation_parallel(
      R = R,
      size_w_in = W_IN,
      size_w_out = W_OUT,
      k_grid = k_grid,
      oos_type = OOS_TYPE,
      compute_weights_fn = compute_weights_fn,
      compute_weights_fn_params = list(),
      n_cores = n_cores,
      return_details = FALSE
    )
  } else {
    sr_vec <- run_oos_evaluation(
      R = R,
      size_w_in = W_IN,
      size_w_out = W_OUT,
      k_grid = k_grid,
      oos_type = OOS_TYPE,
      compute_weights_fn = compute_weights_fn,
      compute_weights_fn_params = list(),
      return_details = FALSE
    )
  }

  SR <- matrix(sr_vec, ncol = 1)
  labels <- METHOD_LABEL

  cat("\n=== Average OOS Sharpe by k ===\n")
  print_results(k_grid, SR, method_labels = labels, digits = 4)

  stem <- sprintf("oos_sr_%s_mebeme_ind_monthly_N%d_Win%d_Wout%d", METHOD_STEM, N, W_IN, W_OUT)
  csv_path <- file.path(OUT_DIR, paste0(stem, ".csv"))
  plot_base <- file.path(FIG_DIR, stem)

  save_results(csv_path, k_grid, SR, method_labels = labels)
  message("Saved results to: ", csv_path)

  plot_sr_empirics(k_grid, SR, save_path = plot_base)
  message("Saved figure to: ", plot_base, ".png")
}
