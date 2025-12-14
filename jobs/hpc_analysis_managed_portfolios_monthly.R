#!/usr/bin/env Rscript
# hpc_analysis_managed_portfolios_monthly.R
# Out-of-sample analysis for MONTHLY managed portfolios on Snellius (HPC-safe).
#
# Usage (from package root):
#   Rscript inst/empirics/hpc_analysis_managed_portfolios_monthly.R
#
# Recommended via SLURM:
#   #SBATCH --cpus-per-task=12
#   export OMP_NUM_THREADS=12 OPENBLAS_NUM_THREADS=12 MKL_NUM_THREADS=12 BLIS_NUM_THREADS=12
#   export SPS_GUROBI_THREADS=1   # if METHOD=miqp and PARALLEL=TRUE
#   Rscript inst/empirics/hpc_analysis_managed_portfolios_monthly.R

# ------------------------------ helpers -------------------------------------

get_script_path <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", ca, value = TRUE)
  if (length(f) == 0L) return(NA_character_)
  sub("^--file=", "", f[1L])
}

find_up <- function(target, start_dir = getwd(), max_up = 20L) {
  d <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
  for (i in seq_len(max_up)) {
    cand <- file.path(d, target)
    if (file.exists(cand)) return(d)
    parent <- normalizePath(file.path(d, ".."), winslash = "/", mustWork = FALSE)
    if (identical(parent, d)) break
    d <- parent
  }
  NA_character_
}

as_int_env <- function(name, default) {
  x <- Sys.getenv(name, unset = "")
  if (x == "") return(as.integer(default))
  suppressWarnings({
    v <- as.integer(x)
    if (is.na(v)) as.integer(default) else v
  })
}

as_lgl_env <- function(name, default) {
  x <- tolower(Sys.getenv(name, unset = ""))
  if (x == "") return(as.logical(default))
  x %in% c("1", "true", "t", "yes", "y")
}

as_chr_env <- function(name, default) {
  x <- Sys.getenv(name, unset = "")
  if (x == "") default else x
}

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

log_msg <- function(...) {
  msg <- paste0("[", timestamp(), "] ", paste0(..., collapse = ""))
  message(msg)
}

# ------------------------------ locate root --------------------------------

script_path <- get_script_path()
script_dir  <- if (!is.na(script_path)) dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE)) else getwd()

# We want to run relative to package root (where DESCRIPTION lives)
pkg_root <- find_up("DESCRIPTION", start_dir = script_dir)
if (is.na(pkg_root)) {
  stop("Could not locate package root (DESCRIPTION not found) starting from: ", script_dir)
}
setwd(pkg_root)
log_msg("Working directory set to package root: ", pkg_root)

# ------------------------------ thread control ------------------------------

# Target CPU budget:
# - Prefer SLURM_CPUS_PER_TASK if present
# - else fall back to 12
Nn <- as_int_env("SLURM_CPUS_PER_TASK", default = 12L)
if (is.na(Nn) || Nn < 1L) Nn <- 12L

# Always set env vars (inherit into parallel workers)
Sys.setenv(
  OMP_NUM_THREADS       = Nn,
  OPENBLAS_NUM_THREADS  = Nn,
  MKL_NUM_THREADS       = Nn,
  BLIS_NUM_THREADS      = Nn
)

suppressPackageStartupMessages({
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(Nn)
    RhpcBLASctl::omp_set_num_threads(Nn)
    log_msg("BLAS threads set to: ", RhpcBLASctl::blas_get_num_procs())
    log_msg("OpenMP max threads set to: ", RhpcBLASctl::omp_get_max_threads())
  } else {
    log_msg("RhpcBLASctl not available; using env vars only (OK on HPC).")
  }
})

# Cap Gurobi threads:
# - default: if PARALLEL is TRUE and METHOD=miqp, use 1 thread (safer)
# - else default to Nn
# - override with SPS_GUROBI_THREADS
# NOTE: If you run MIQP in parallel, using 1 thread per worker is usually best.
#       If you run MIQP sequentially, you can set SPS_GUROBI_THREADS=Nn.
#       (You can also make it smaller than Nn.)
# We'll set MIQP_THREADS after PARALLEL/METHOD are defined (below).

# ------------------------------ load package --------------------------------

suppressPackageStartupMessages({
  library(SparsePortfolioSelection)
})

# ------------------------------ configuration --------------------------------
# You can override key options via environment variables if desired:
#   SPS_METHOD=lasso|elnet|miqp
#   SPS_REFIT=true|false
#   SPS_PARALLEL=true|false
#   SPS_N_ASSETS=250, SPS_W_IN=480, SPS_W_OUT=1, SPS_K_MIN=3, SPS_K_STEP=2, SPS_K_CAP=...
#   SPS_RNG_SEED=12345
#   SPS_OOS_TYPE=rolling|expanding
#   SPS_MISSINGS=median|...
#   SPS_GUROBI_THREADS=1|...

PANEL_TYPE <- "US"                        # ignored for monthly frequency
MISSINGS   <- as_chr_env("SPS_MISSINGS",  "median")
N_ASSETS   <- as_int_env("SPS_N_ASSETS",  250L)   # total in monthly panel ~353
RNG_SEED   <- as_int_env("SPS_RNG_SEED",  12345L)

W_IN       <- as_int_env("SPS_W_IN",      480L)   # months: [240 360 480]
W_OUT      <- as_int_env("SPS_W_OUT",     1L)     # months: [1]
OOS_TYPE   <- as_chr_env("SPS_OOS_TYPE",  "rolling")

K_MIN      <- as_int_env("SPS_K_MIN",     3L)
K_STEP     <- as_int_env("SPS_K_STEP",    2L)
# default K_CAP = N_ASSETS - 1 (applied after subsampling)
K_CAP_ENV  <- Sys.getenv("SPS_K_CAP", unset = "")
K_CAP      <- if (K_CAP_ENV == "") NA_integer_ else suppressWarnings(as.integer(K_CAP_ENV))

METHOD     <- as_chr_env("SPS_METHOD",    "lasso")  # "lasso" | "elnet" | "miqp"
REFIT      <- as_lgl_env("SPS_REFIT",     FALSE)
PARALLEL   <- as_lgl_env("SPS_PARALLEL",  TRUE)

if (!METHOD %in% c("lasso", "elnet", "miqp")) stop("Unknown METHOD: ", METHOD)
if (!OOS_TYPE %in% c("rolling", "expanding")) stop("Unknown OOS_TYPE: ", OOS_TYPE)

# Decide method label/stem
refit_suffix <- if (REFIT) "_refit" else ""
METHOD_LABEL <- paste0(METHOD, refit_suffix)
METHOD_STEM  <- paste0(METHOD, refit_suffix)

# Output locations (relative to package root)
OUT_DIR <- file.path("inst", "empirics", "results", "managed_portfolios_monthly")
FIG_DIR <- file.path("inst", "empirics", "figures", "managed_portfolios_monthly")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# MIQP thread cap logic (now that METHOD/PARALLEL are known)
default_miqp_threads <- if (METHOD == "miqp" && isTRUE(PARALLEL)) 1L else Nn
MIQP_THREADS <- as_int_env("SPS_GUROBI_THREADS", default = default_miqp_threads)
if (is.na(MIQP_THREADS) || MIQP_THREADS < 1L) MIQP_THREADS <- 1L

log_msg("Config: METHOD=", METHOD_LABEL,
        ", PARALLEL=", PARALLEL,
        ", Nn=", Nn,
        ", MIQP_THREADS=", MIQP_THREADS,
        ", W_IN=", W_IN,
        ", W_OUT=", W_OUT,
        ", OOS_TYPE=", OOS_TYPE)

# ------------------------------ data loading --------------------------------

if (!is.null(RNG_SEED)) set.seed(RNG_SEED)

R <- load_data(type = PANEL_TYPE, missing = MISSINGS, path = "data", frequency = "monthly")

T_full <- nrow(R)
N_full <- ncol(R)
log_msg("Loaded monthly panel: T_full=", T_full, ", N_full=", N_full)

# Subsample assets (reproducible)
N <- min(N_ASSETS, N_full)
asset_idx <- sort(sample.int(N_full, N))
R <- R[, asset_idx, drop = FALSE]

# Persist asset subset so reruns are identical even if you tweak other things
saveRDS(asset_idx, file.path(OUT_DIR, sprintf("asset_idx_monthly_N%d_seed%d.rds", N, RNG_SEED)))

Tobs <- nrow(R)
N    <- ncol(R)

k_cap_default <- N - 1L
k_cap_user <- if (is.na(K_CAP)) k_cap_default else min(K_CAP, k_cap_default)
k_max <- min(k_cap_user, N - 1L)
if (k_max < K_MIN) stop("k_max < K_MIN; reduce K_MIN or increase N.")
k_grid <- seq.int(K_MIN, k_max, by = K_STEP)

log_msg(sprintf("Working sample: T=%d, N=%d, k ∈ [%d..%d] step %d",
                Tobs, N, K_MIN, k_max, K_STEP))

# ------------------------------ parameter lists -----------------------------

alpha_grid <- seq(0.30, 1.00, by = 0.05)
if (METHOD == "lasso") alpha_grid <- 1.00

lasso_params <- list(
  nlambda          = 100L,
  lambda_min_ratio = 1e-3,
  alpha            = alpha_grid,
  n_folds          = 5L,
  nadd             = 80L,
  nnested          = 3L,
  standardize      = FALSE,
  stabilize_sigma  = TRUE,
  compute_weights  = TRUE,
  normalize_weights= FALSE,
  use_refit        = REFIT
)

miqp_params <- list(
  exactly_k        = TRUE,
  m                = 1L,
  gamma            = 1.0,
  fmin             = -0.25,
  fmax             = 0.25,
  expand_rounds    = 10L,
  expand_factor    = 3.0,
  expand_tol       = 1e-2,
  mipgap           = 1e-4,
  time_limit       = 100,
  threads          = MIQP_THREADS,
  compute_weights  = TRUE,
  normalize_weights= FALSE,
  use_refit        = REFIT,
  verbose          = FALSE,
  stabilize_sigma  = TRUE
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
    res <- do.call(
      mve_lasso_search,
      c(list(mu = mu, sigma = sigma, n_obs = nrow(Rin), k = k), lasso_params)
    )
    list(weights = res$weights, selection = res$selection, status = res$status)
  }
}

# ------------------------------ run OOS -------------------------------------

log_msg(sprintf("Starting OOS run: T=%d, N=%d, W_IN=%d, W_OUT=%d",
                Tobs, N, W_IN, W_OUT))

# Conservative parallel core choice:
# - Prefer SLURM_CPUS_PER_TASK (Nn)
# - Never exceed detectCores-1
detect_cap <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
n_cores <- min(Nn, detect_cap)

if (PARALLEL && n_cores >= 2L) {
  log_msg("Parallel enabled with n_cores=", n_cores, " (cap Nn=", Nn, ", detect_cap=", detect_cap, ")")
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
  if (PARALLEL) log_msg("Parallel requested but n_cores<2; falling back to sequential.")
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

stem <- sprintf("oos_sr_%s_monthly_N%d_Win%d_Wout%d", METHOD_STEM, N, W_IN, W_OUT)
csv_path  <- file.path(OUT_DIR, paste0(stem, ".csv"))
plot_base <- file.path(FIG_DIR, stem)

save_results(csv_path, k_grid, SR, method_labels = labels)
log_msg("Saved results to: ", csv_path)

plot_sr_empirics(k_grid, SR, save_path = plot_base)
log_msg("Saved figure to: ", plot_base, ".png")

# Write a small run manifest for reproducibility
manifest <- list(
  timestamp      = timestamp(),
  R_version      = R.version.string,
  method         = METHOD_LABEL,
  parallel       = PARALLEL,
  n_threads_env  = list(
    Nn = Nn,
    OMP_NUM_THREADS = Sys.getenv("OMP_NUM_THREADS"),
    OPENBLAS_NUM_THREADS = Sys.getenv("OPENBLAS_NUM_THREADS"),
    MKL_NUM_THREADS = Sys.getenv("MKL_NUM_THREADS"),
    BLIS_NUM_THREADS = Sys.getenv("BLIS_NUM_THREADS"),
    SPS_GUROBI_THREADS = Sys.getenv("SPS_GUROBI_THREADS", unset = "")
  ),
  data           = list(type = PANEL_TYPE, missing = MISSINGS, frequency = "monthly"),
  dims           = list(T = Tobs, N = N),
  windows        = list(W_IN = W_IN, W_OUT = W_OUT, oos_type = OOS_TYPE),
  k_grid         = k_grid,
  rng_seed       = RNG_SEED,
  asset_idx_rds  = file.path(OUT_DIR, sprintf("asset_idx_monthly_N%d_seed%d.rds", N, RNG_SEED)),
  output_csv     = csv_path,
  output_plot    = paste0(plot_base, ".png")
)
saveRDS(manifest, file.path(OUT_DIR, paste0(stem, "_manifest.rds")))
log_msg("Saved manifest to: ", file.path(OUT_DIR, paste0(stem, "_manifest.rds")))

log_msg("Done.")
