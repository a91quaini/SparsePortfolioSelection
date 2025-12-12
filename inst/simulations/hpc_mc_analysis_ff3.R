#!/usr/bin/env Rscript
# hpc_simulation_ff3_mc.R
# Monte Carlo analysis for FF3-style simulations (HPC/Snellius-safe).
#
# Run from package root:
#   Rscript inst/simulations/hpc_simulation_ff3_mc.R
#
# Recommended via SLURM:
#   #SBATCH --cpus-per-task=192
#   export OMP_NUM_THREADS=192 OPENBLAS_NUM_THREADS=192 MKL_NUM_THREADS=192 BLIS_NUM_THREADS=192
#   export SPS_PARALLEL=false            # this script is sequential by default
#   export SPS_METHOD=miqp               # miqp or lasso
#   export SPS_GUROBI_THREADS=192        # if sequential MIQP; set 1 if you parallelize MC yourself
#   Rscript inst/simulations/hpc_simulation_ff3_mc.R

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
log_msg <- function(...) message(paste0("[", timestamp(), "] ", paste0(..., collapse = "")))

# ------------------------------ locate root --------------------------------

script_path <- get_script_path()
script_dir  <- if (!is.na(script_path)) dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE)) else getwd()

pkg_root <- find_up("DESCRIPTION", start_dir = script_dir)
if (is.na(pkg_root)) {
  stop("Could not locate package root (DESCRIPTION not found) starting from: ", script_dir)
}
setwd(pkg_root)
log_msg("Working directory set to package root: ", pkg_root)

# ------------------------------ thread control ------------------------------

Nn <- as_int_env("SLURM_CPUS_PER_TASK", default = 192L)
if (is.na(Nn) || Nn < 1L) Nn <- 192L

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

# ------------------------------ load package --------------------------------

suppressPackageStartupMessages({
  library(SparsePortfolioSelection)
})

# ------------------------------ configuration --------------------------------
# Override via env vars:
#   SPS_SEED=123
#   SPS_N_ASSETS=100
#   SPS_N_OBS=480
#   SPS_N_MC=200
#   SPS_METHOD=miqp|lasso
#   SPS_K_MIN=2
#   SPS_K_MAX=99          (default: N-1)
#   SPS_K_STEP=1
#   SPS_GUROBI_THREADS=12 (defaults to Nn if METHOD=miqp)
#
# NOTE: This script is sequential. If you parallelize MC externally (e.g., SLURM array),
#       set SPS_GUROBI_THREADS to match your per-task CPU budget (often Nn), or 1 if needed.

SEED         <- as_int_env("SPS_SEED",     123L)

n_assets     <- as_int_env("SPS_N_ASSETS", 100L)   # [25, 50, 100]
n_obs        <- as_int_env("SPS_N_OBS",    480L)   # [120, 240, 480]
n_MC         <- as_int_env("SPS_N_MC",     200L)
search_method<- as_chr_env("SPS_METHOD",  "miqp")  # "lasso" or "miqp"

k_min        <- as_int_env("SPS_K_MIN",    2L)
k_max_env    <- Sys.getenv("SPS_K_MAX", unset = "")
k_max        <- if (k_max_env == "") (n_assets - 1L) else suppressWarnings(as.integer(k_max_env))
k_step       <- as_int_env("SPS_K_STEP",   1L)

if (!search_method %in% c("lasso", "miqp")) stop("SPS_METHOD must be 'lasso' or 'miqp'")
if (n_assets < 3L) stop("n_assets must be >= 3")
if (n_obs < 10L) stop("n_obs too small")
if (n_MC < 1L) stop("n_MC must be >= 1")

k_max <- min(k_max, n_assets - 1L)
if (k_max < k_min) stop("k_max < k_min; adjust SPS_K_MIN/SPS_K_MAX.")
k_grid <- seq.int(k_min, k_max, by = k_step)

# MIQP thread cap: default to Nn (sequential script)
default_miqp_threads <- Nn
MIQP_THREADS <- as_int_env("SPS_GUROBI_THREADS", default = default_miqp_threads)
if (is.na(MIQP_THREADS) || MIQP_THREADS < 1L) MIQP_THREADS <- 1L

set.seed(SEED)

log_msg("Config: method=", search_method,
        ", N=", n_assets,
        ", T=", n_obs,
        ", n_MC=", n_MC,
        ", k=[", k_min, "..", k_max, "] step ", k_step,
        ", CPU cap Nn=", Nn,
        ", MIQP_THREADS=", MIQP_THREADS)

# ------------------------------ search params --------------------------------

lasso_params <- list(
  n_obs             = n_obs,
  nlambda           = 100L,
  lambda_min_ratio  = 1e-3,
  nadd              = 80L,
  nnested           = 2L,
  alpha             = 1,
  standardize       = FALSE,
  compute_weights   = TRUE,
  normalize_weights = FALSE,
  use_refit         = FALSE
)

miqp_params <- list(
  exactly_k         = TRUE,
  gamma             = 1.0,
  fmin              = 0,
  fmax              = 0.25,
  mipgap            = 1e-4,
  time_limit        = 60,
  threads           = MIQP_THREADS,
  compute_weights   = TRUE,
  normalize_weights = FALSE,
  verbose           = FALSE
)

if (search_method == "lasso") {
  mve_search_fn <- mve_lasso_search
  mve_search_fn_params <- lasso_params
} else {
  mve_search_fn <- mve_miqp_search
  mve_search_fn_params <- miqp_params
}

# ------------------------------ calibrate population -------------------------

params <- calibrate_ff3(n_assets)
mu_pop <- params$mu_true
sigma_pop <- params$sigma_true

# Population SR across k (use robust boosted computation inside package)
pop_sr <- compute_population_mve_sr(
  mu_pop = mu_pop,
  sigma_pop = sigma_pop,
  k_grid = k_grid,
  n_obs = n_obs,
  lasso_nlambda = 300L,
  lasso_lambda_min_ratio = 1e-4,
  miqp_time_limit = 200,
  miqp_mipgap = 1e-5,
  verbose = FALSE
)

# ------------------------------ Monte Carlo ----------------------------------

est_terms <- matrix(NA_real_, nrow = n_MC, ncol = length(k_grid))
sel_terms <- matrix(NA_real_, nrow = n_MC, ncol = length(k_grid))
colnames(est_terms) <- paste0("k", k_grid)
colnames(sel_terms) <- paste0("k", k_grid)

for (mc in seq_len(n_MC)) {
  if (mc == 1L || mc %% 10L == 0L) log_msg("MC run ", mc, " / ", n_MC)
  R <- simulate_ff3(n_obs, params)
  mu_sample <- colMeans(R)
  sigma_sample <- stats::cov(R)

  for (i in seq_along(k_grid)) {
    k <- k_grid[i]
    decomp <- compute_mve_sr_decomposition(
      mu_pop, sigma_pop,
      mu_sample, sigma_sample,
      k = k,
      mve_search_fn = mve_search_fn,
      mve_search_fn_params = mve_search_fn_params,
      do_checks = FALSE,
      return_selection = FALSE
    )
    est_terms[mc, i] <- decomp$mve_sr_cardk_est_term
    sel_terms[mc, i] <- decomp$mve_sr_cardk_sel_term
  }
}

# ------------------------------ save outputs ---------------------------------

RES_DIR <- file.path("inst", "simulations", "results")
FIG_DIR <- file.path("inst", "simulations", "figures")
dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

stem <- sprintf("ff3_%s_N%d_T%d_MC%d_seed%d_k%d_%d_step%d",
                search_method, n_assets, n_obs, n_MC, SEED, k_min, k_max, k_step)

fname_rds <- file.path(RES_DIR, paste0("results_", stem, ".rds"))

result_obj <- list(
  config = list(
    n_assets = n_assets,
    n_obs = n_obs,
    k_grid = k_grid,
    n_MC = n_MC,
    search_method = search_method,
    seed = SEED,
    cpu = list(
      Nn = Nn,
      MIQP_THREADS = MIQP_THREADS,
      OMP_NUM_THREADS = Sys.getenv("OMP_NUM_THREADS"),
      OPENBLAS_NUM_THREADS = Sys.getenv("OPENBLAS_NUM_THREADS"),
      MKL_NUM_THREADS = Sys.getenv("MKL_NUM_THREADS"),
      BLIS_NUM_THREADS = Sys.getenv("BLIS_NUM_THREADS")
    ),
    mve_search_fn_params = mve_search_fn_params
  ),
  population_sr = pop_sr,
  est_terms = est_terms,
  sel_terms = sel_terms
)

saveRDS(result_obj, file = fname_rds)
log_msg("Saved results to: ", fname_rds)

# Plot (let the package decide file naming/locations, but ensure working dir is pkg root)
# If plot_mve_sr_decomposition() supports save paths, adapt accordingly; otherwise it saves internally.
plots <- plot_mve_sr_decomposition(result_obj, save = TRUE)
log_msg("Plotting complete (save=TRUE).")

log_msg("Done.")
