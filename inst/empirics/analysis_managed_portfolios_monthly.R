# Out-of-sample analysis for monthly managed portfolios using run_oos_evaluation.
# Mirrors the daily script but loads monthly data via load_data(..., frequency = "monthly").

## ---- thread control: must be at the very top ------------------------------
# Nn = 1L
Nn = 96L
# Nn = min(Nn, parallel::detectCores(logical = TRUE) - 1L)
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
#
# # Cap Gurobi threads (default to BLAS cap; override via SPS_GUROBI_THREADS)
# MIQP_THREADS <- as.integer(Sys.getenv("SPS_GUROBI_THREADS", Nn))
# if (is.na(MIQP_THREADS) || MIQP_THREADS < 1L) MIQP_THREADS <- 1L

library(SparsePortfolioSelection)

# Configuration: 636 observations
PANEL_TYPE <- "US"      # ignored for monthly frequency
MISSINGS <- "median"    # how to treat missing values
N_ASSETS <- 200         # subset of assets to use: total = 355
RNG_SEED <- 12345
W_IN_GRID <- c(240L, 360L, 480L)  # in-sample lengths (months)
W_OUT <- 1             # OOS block length (months): [1]
OOS_TYPE <- "rolling"   # "rolling" or "expanding"
ADD_MKT <- TRUE         # append MKT-RF
ADD_FACTORS <- TRUE    # append FF3 (MKT, SMB, HML)
COMPLETE_ANALYSIS <- TRUE  # if TRUE, run complete analysis (turnover/instability)
CHECK_K <- TRUE         # warn if solver returns sparsity different from k
K_TOL <- 1e-9           # tolerance for nonzero weights when checking sparsity
K_MIN <- 3
K_STEP <- 3
K_CAP <- N_ASSETS - 1
METHOD <- "lasso"        # "lasso" | "elnet" | "miqp"
REFIT <- FALSE
PARALLEL <- TRUE

# Decide filename/label stems (append _refit if refit enabled)
refit_suffix <- if ((METHOD %in% c("lasso", "elnet") && REFIT) ||
                    (METHOD == "miqp" && REFIT)) "_refit" else ""
METHOD_LABEL <- paste0(METHOD, refit_suffix)
METHOD_STEM <- paste0(METHOD, refit_suffix)

OUT_DIR <- file.path("inst", "empirics", "results", "managed_portfolios_monthly")
FIG_DIR <- file.path("inst", "empirics", "figures", "managed_portfolios_monthly")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# Fix RNG before any shuffling inside load_data
if (!is.null(RNG_SEED)) set.seed(RNG_SEED)

ld <- load_data(type = PANEL_TYPE, missing = MISSINGS, path = "data", frequency = "monthly",
                add_mkt = ADD_MKT, add_factors = ADD_FACTORS)
R_all <- ld$returns
rf_vec <- if (is.null(ld$rf)) rep(0, nrow(R_all)) else ld$rf
T_full <- nrow(R_all); N_full <- ncol(R_all)

# Select a subset of assets for speed/reproducibility; keep any appended factors
N <- min(N_ASSETS, N_full)
coln <- colnames(R_all)
factor_idx <- integer(0)
if (!is.null(coln)) {
  factor_idx <- which(grepl("MKT", coln) | grepl("SMB", coln) | grepl("HML", coln))
}
asset_idx <- setdiff(seq_len(ncol(R_all)), factor_idx)
asset_idx <- asset_idx[seq_len(min(N, length(asset_idx)))]
keep_idx <- c(asset_idx, factor_idx)
R_all <- R_all[, keep_idx, drop = FALSE]
if (!is.null(rf_vec)) rf_vec <- rf_vec[seq_len(nrow(R_all))]

k_grid <- NULL  # defined per run once N is known
kg_ref <- NULL  # reference k-grid for combined plots
sr_list <- list()
turn_list <- list()
instab1_list <- list()
instab2_list <- list()
selinst_list <- list()
lessk_list <- list()
lessk_totals <- integer(0)

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

# miqp_params <- list(
#   exactly_k = TRUE,
#   m = 1L,
#   gamma = 1.0,
#   fmin = -0.25,
#   fmax = 0.25,
#   expand_rounds = 10L,
#   expand_factor = 3.0,
#   expand_tol = 1e-2,
#   mipgap = 1e-4,
#   time_limit = 100,
#   threads = MIQP_THREADS,
#   compute_weights = TRUE,
#   normalize_weights = FALSE,
#   use_refit = REFIT,
#   verbose = FALSE,
#   stabilize_sigma = TRUE
# )

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
  if (!is.null(res$selections_path)) {
    for (ik in seq_along(k_grid)) {
      sels <- res$selections_path[[ik]]
      if (length(sels) > 0) {
        nz_counts <- vapply(sels, function(sel) length(sel), integer(1))
        less_than_k[ik] <- sum(nz_counts < k_grid[ik])
      }
    }
  }
  factors_tag <- if (ADD_FACTORS) "ff3" else "nofactors"
  mkt_tag <- if (ADD_MKT) "mkt" else "nomkt"
  stem_base <- sprintf("combined_oos_%s_%s_%s_%s_%s_Wout%d",
                       METHOD,
                       if (REFIT) "refit" else "norefit",
                       "us",
                       factors_tag,
                       mkt_tag,
                       W_OUT)
  res_table <- data.frame(
    k = k_grid,
    SharpeRatio = SR[, 1],
    less_than_k = less_than_k,
    turnover = res$summary$median_turnover,
    selection_instability = res$summary$mean_selection_instability,
    weight_instability_l1 = res$summary$median_weight_instability_L1,
    weight_instability_l2 = res$summary$median_weight_instability_L2
  )
  write.csv(res_table, file.path(OUT_DIR, paste0(stem_base, "_Win", W_IN, "_sr.csv")), row.names = FALSE)
  sr_list[[length(sr_list) + 1L]] <- SR[, 1]
  turn_list[[length(turn_list) + 1L]] <- res$summary$median_turnover
  instab1_list[[length(instab1_list) + 1L]] <- res$summary$median_weight_instability_L1
  instab2_list[[length(instab2_list) + 1L]] <- res$summary$median_weight_instability_L2
  selinst_list[[length(selinst_list) + 1L]] <- res$summary$mean_selection_instability
  lessk_list[[length(lessk_list) + 1L]] <- less_than_k
  lessk_totals[length(lessk_totals) + 1L] <- length(res$oos_returns[[1]])
}

# Combined plots across W_IN (if multiple)
if (COMPLETE_ANALYSIS && length(sr_list) >= 1) {
  comb_labels <- as.character(W_IN_GRID)
  sr_mat <- do.call(cbind, sr_list)
  panel_tag <- "us"
  factors_tag <- if (ADD_FACTORS) "ff3" else "nofactors"
  mkt_tag <- if (ADD_MKT) "mkt" else "nomkt"
  base_stem <- sprintf("combined_oos_%s_%s_%s_%s_%s_Wout%d",
                       METHOD,
                       if (REFIT) "refit" else "norefit",
                       panel_tag,
                       factors_tag,
                       mkt_tag,
                       W_OUT)
  plot_sr_empirics(k_grid, sr_mat, labels = comb_labels,
                   save_path = file.path(FIG_DIR, paste0(base_stem, "_sr")))
  if (length(turn_list) == length(sr_list)) {
    turn_mat <- do.call(cbind, turn_list)
    plot_turnover_empirics(k_grid, turn_mat, labels = comb_labels,
                           save_path = file.path(FIG_DIR, paste0(base_stem, "_turnover")))
  }
  if (length(instab1_list) == length(sr_list)) {
    instab1_mat <- do.call(cbind, instab1_list)
    instab2_mat <- do.call(cbind, instab2_list)
    plot_weight_instability_empirics(k_grid, instab1_mat, instab2_mat,
                                     labels = comb_labels,
                                     save_path_base = file.path(FIG_DIR, paste0(base_stem, "_weight_instability")))
  }
  if (length(selinst_list) == length(sr_list)) {
    sel_mat <- do.call(cbind, selinst_list)
    plot_selection_instability_empirics(k_grid, sel_mat, labels = comb_labels,
                                        save_path = file.path(FIG_DIR, paste0(base_stem, "_selection_instability")))
  }
  if (length(lessk_list) == length(sr_list)) {
    lessk_mat <- do.call(cbind, lessk_list)
    plot_less_than_k(k_grid, lessk_mat, labels = comb_labels,
                     save_path = file.path(FIG_DIR, paste0(base_stem, "_less_than_k")),
                     total_windows = lessk_totals)
  }
}
