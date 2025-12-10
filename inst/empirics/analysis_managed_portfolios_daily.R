# Out-of-sample analysis for daily managed portfolios using run_oos_evaluation.
# This script:
#  - loads processed daily panels (US or International) from data/,
#  - optionally imputes/handles missing values,
#  - randomly subselects a universe of assets,
#  - runs rolling/expanding OOS Sharpe evaluation for a chosen solver (lasso/elnet/miqp),
#  - saves CSV results and plots under inst/empirics/{results,figures}.

## ---- thread control: must be at the very top ------------------------------
Nn = 1L
# Nn = 12L
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

# Configuration
PANEL_TYPE <- "US"      # "US" or "International"
MISSINGS <- "median"    # how to treat missing values
N_ASSETS <- 274         # subset of assets to use -> 274 for "US" and 400 for "International"
RNG_SEED <- 12345
W_IN <- 252             # in-sample length (days)
W_OUT <- 30             # OOS block length (non-overlapping)
OOS_TYPE <- "rolling"   # "rolling" or "expanding"
K_MIN <- 3
K_STEP <- 9
K_CAP <- N_ASSETS - 1
METHOD <- "miqp"    # "lasso" | "elnet" | "miqp"
REFIT <- TRUE
PARALLEL <- FALSE 

# Decide filename/label stems (append _refit if refit enabled)
refit_suffix <- if ((METHOD %in% c("lasso", "elnet") && REFIT) ||
                    (METHOD == "miqp" && REFIT)) "_refit" else ""
METHOD_LABEL <- paste0(METHOD, refit_suffix)
METHOD_STEM <- paste0(METHOD, refit_suffix)

OUT_DIR <- file.path("inst", "empirics", "results", "managed_portfolios_daily")
FIG_DIR <- file.path("inst", "empirics", "figures", "managed_portfolios_daily")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

R <- load_data(type = PANEL_TYPE, missing = MISSINGS, path = "data")
T_full <- nrow(R); N_full <- ncol(R)

# Select a subset of assets for speed/reproducibility
if (!is.null(RNG_SEED)) set.seed(RNG_SEED)
N <- min(N_ASSETS, N_full)
asset_idx <- sort(sample.int(N_full, N))
R <- R[, asset_idx, drop = FALSE]

Tobs <- nrow(R); N <- ncol(R)
k_max <- min(K_CAP, N - 1)
if (k_max < K_MIN) stop("k_max < K_MIN; reduce K_MIN or increase N.")
k_grid <- seq.int(K_MIN, k_max, by = K_STEP)

# Define parameter lists
alpha_grid = seq(0.30, 1.00, by = 0.05)
if (METHOD == "lasso") {
  alpha_grid = 1.00
}
lasso_params <- list(
  nlambda = 100L, # 300L,
  lambda_min_ratio = 1e-3, # 1e-2,
  alpha = alpha_grid,
  n_folds = 5L,
  nadd = 80L,
  nnested = 3L,
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

message(sprintf("Starting OOS run: T=%d, N=%d, W_IN=%d, W_OUT=%d, k ∈ [%d..%d]", Tobs, N, W_IN, W_OUT, K_MIN, k_max))

# Optional parallel run: set PARALLEL <- TRUE to enable
# n_threads = parallel::detectCores(logical = TRUE) - 1L
n_threads = 12
if (PARALLEL) {
  n_cores <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
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

stem <- sprintf("oos_sr_%s_%s_N%d_Win%d_Wout%d", METHOD_STEM, PANEL_TYPE, N, W_IN, W_OUT)
csv_path <- file.path(OUT_DIR, paste0(stem, ".csv"))
plot_base <- file.path(FIG_DIR, stem)

save_results(csv_path, k_grid, SR, method_labels = labels)
message("Saved results to: ", csv_path)

plot_sr_empirics(k_grid, SR, save_path = plot_base)
message("Saved figure to: ", plot_base, ".png")
