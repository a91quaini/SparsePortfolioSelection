# Monte Carlo analysis for FF3-style simulations with high-quality population SR.
# Population SR per k is max of: LASSO, LASSO+refit, and a boosted MIQP.
# Sample estimates use user-chosen method (LASSO or MIQP) with user-tuned params.

library(SparsePortfolioSelection)

set.seed(123)
# HPC-friendly defaults: limit BLAS/OMP to 1 and choose worker count via SPS_CORES (default 64)
Sys.setenv(
  OMP_NUM_THREADS = 1L,
  OPENBLAS_NUM_THREADS = 1L,
  MKL_NUM_THREADS = 1L,
  BLIS_NUM_THREADS = 1L
)
Nn <- as.integer(Sys.getenv("SPS_CORES", "64"))
if (is.na(Nn) || Nn < 1L) Nn <- 64L
det_cores <- parallel::detectCores(logical = TRUE)
n_cores <- max(1L, min(Nn, if (is.na(det_cores)) Nn else det_cores))
use_psock <- identical(Sys.getenv("SPS_FORK"), "0")

# Simulation parameters
n_assets <- 100
n_obs_grid <- c(120L, 240L, 480L, 960L)   # sample sizes
k_grid <- 2:(n_assets - 1L)
n_MC <- 1000

# User-chosen sample method: "lasso" or "miqp"
search_method <- "lasso"

# Sample LASSO params (adjust as desired)
lasso_sample_base <- list(
  nlambda = 300L,
  lambda_min_ratio = 1e-3,
  nadd = 100L,
  nnested = 3L,
  alpha = 1,
  standardize = FALSE,
  compute_weights = TRUE,
  normalize_weights = FALSE,
  use_refit = FALSE
)

# Sample MIQP params (adjust as desired)
miqp_sample <- list(
  exactly_k = TRUE,
  gamma = 1.0,
  fmin = -0.25,
  fmax = 0.25,
  mipgap = 1e-4,
  time_limit = 200,
  threads = max(1L, min(8L, n_cores)),
  compute_weights = TRUE,
  normalize_weights = FALSE,
  verbose = FALSE
)

# Population search params (aim for quality, not speed)
lasso_pop_base <- list(
  nlambda = 400L,
  lambda_min_ratio = 1e-2,
  nadd = 100L,
  nnested = 5L,
  alpha = 1,
  standardize = FALSE,
  compute_weights = TRUE,
  normalize_weights = FALSE,
  use_refit = FALSE
)
lasso_pop_refit <- modifyList(lasso_pop_base, list(use_refit = TRUE))
miqp_pop <- list(
  exactly_k = FALSE,
  gamma = 1.0,
  fmin = -0.25,
  fmax =0.25,
  expand_rounds = 20L,
  expand_factor = 2.0,
  expand_tol = 1e-4,
  mipgap = 5e-6,
  time_limit = 500,
  threads = max(1L, min(8L, n_cores)),
  compute_weights = TRUE,
  normalize_weights = FALSE,
  use_refit = FALSE,
  verbose = FALSE
)

# Calibrate population parameters
params <- calibrate_ff3(n_assets)
mu_pop <- params$mu_true
sigma_pop <- params$sigma_true

# Population SR per k: max over LASSO, LASSO+refit, boosted MIQP
pop_sr <- numeric(length(k_grid))
for (i in seq_along(k_grid)) {
  k <- k_grid[i]

  # LASSO (no refit)
  res_lasso <- mve_lasso_search(
    mu = mu_pop,
    sigma = sigma_pop,
    k = k,
    n_obs = n_obs_grid[[1]],
    nlambda = lasso_pop_base$nlambda,
    lambda_min_ratio = lasso_pop_base$lambda_min_ratio,
    alpha = lasso_pop_base$alpha,
    nadd = lasso_pop_base$nadd,
    nnested = lasso_pop_base$nnested,
    standardize = lasso_pop_base$standardize,
    compute_weights = TRUE,
    normalize_weights = FALSE,
    use_refit = lasso_pop_base$use_refit
  )
  sr_lasso <- as.numeric(res_lasso$sr)

  # LASSO + refit
  res_lasso_refit <- mve_lasso_search(
    mu = mu_pop,
    sigma = sigma_pop,
    k = k,
    n_obs = n_obs_grid[[1]],
    nlambda = lasso_pop_refit$nlambda,
    lambda_min_ratio = lasso_pop_refit$lambda_min_ratio,
    alpha = lasso_pop_refit$alpha,
    nadd = lasso_pop_refit$nadd,
    nnested = lasso_pop_refit$nnested,
    standardize = lasso_pop_refit$standardize,
    compute_weights = TRUE,
    normalize_weights = FALSE,
    use_refit = lasso_pop_refit$use_refit
  )
  sr_lasso_refit <- as.numeric(res_lasso_refit$sr)

  # MIQP boosted
  res_miqp <- mve_miqp_search(
    mu_pop, sigma_pop, k,
    exactly_k = miqp_pop$exactly_k,
    gamma = miqp_pop$gamma,
    fmin = miqp_pop$fmin,
    fmax = miqp_pop$fmax,
    expand_rounds = miqp_pop$expand_rounds,
    expand_factor = miqp_pop$expand_factor,
    expand_tol = miqp_pop$expand_tol,
    mipgap = miqp_pop$mipgap,
    time_limit = miqp_pop$time_limit,
    threads = miqp_pop$threads,
    compute_weights = TRUE,
    normalize_weights = FALSE,
    use_refit = miqp_pop$use_refit,
    verbose = miqp_pop$verbose
  )
  sr_miqp <- as.numeric(res_miqp$sr)

  pop_sr[i] <- max(sr_lasso, sr_lasso_refit, sr_miqp, na.rm = TRUE)
}

dir.create("inst/simulations/results", recursive = TRUE, showWarnings = FALSE)

for (n_obs in n_obs_grid) {
  lasso_sample <- modifyList(lasso_sample_base, list(n_obs = n_obs))

  if (search_method == "lasso") {
    mve_search_fn <- mve_lasso_search
    mve_search_fn_params <- lasso_sample
  } else if (search_method == "miqp") {
    mve_search_fn <- mve_miqp_search
    mve_search_fn_params <- miqp_sample
  } else {
    stop("search_method must be 'lasso' or 'miqp'")
  }

  run_one <- function(mc_id) {
    R <- simulate_ff3(n_obs, params)
    mu_sample <- colMeans(R)
    sigma_sample <- stats::cov(R)
    est_vec <- numeric(length(k_grid))
    sel_vec <- numeric(length(k_grid))
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
      est_vec[i] <- decomp$mve_sr_cardk_est_term
      sel_vec[i] <- decomp$mve_sr_cardk_sel_term
    }
    list(est = est_vec, sel = sel_vec)
  }

  if (n_cores == 1L || use_psock) {
    if (n_cores == 1L) {
      res_list <- lapply(seq_len(n_MC), run_one)
    } else {
      cl <- parallel::makeCluster(n_cores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      parallel::clusterEvalQ(cl, {
        Sys.setenv(OMP_NUM_THREADS = 1L,
                   OPENBLAS_NUM_THREADS = 1L,
                   MKL_NUM_THREADS = 1L,
                   BLIS_NUM_THREADS = 1L)
        library(SparsePortfolioSelection)
        NULL
      })
      parallel::clusterExport(cl, c("params", "k_grid", "miqp_sample", "lasso_sample", "miqp_pop",
                                    "mu_pop", "sigma_pop", "pop_sr", "n_obs", "search_method",
                                    "mve_lasso_search", "mve_miqp_search", "compute_mve_sr_decomposition",
                                    "simulate_ff3", "n_MC", "run_one"),
                              envir = environment())
      res_list <- parallel::parLapply(cl, seq_len(n_MC), run_one)
    }
  } else {
    res_list <- parallel::mclapply(seq_len(n_MC), run_one, mc.cores = n_cores)
  }

  est_terms <- t(vapply(res_list, `[[`, numeric(length(k_grid)), "est"))
  sel_terms <- t(vapply(res_list, `[[`, numeric(length(k_grid)), "sel"))
  colnames(est_terms) <- paste0("k", k_grid)
  colnames(sel_terms) <- paste0("k", k_grid)

  fname <- sprintf("inst/simulations/results/results_ff3_%s_N%d_T%d_hpc_simulation_analysis_ff3.rds",
                   search_method, n_assets, n_obs)
  saveRDS(list(
    config = list(
      n_assets = n_assets,
      n_obs = n_obs,
      k_grid = k_grid,
      n_MC = n_MC,
      search_method = search_method,
      mve_search_fn_params = mve_search_fn_params
    ),
    population_sr = pop_sr,
    est_terms = est_terms,
    sel_terms = sel_terms
  ), file = fname)
  message(sprintf("Saved results to %s", fname))

  plot_mve_sr_decomposition(list(
    config = list(
      n_assets = n_assets,
      n_obs = n_obs,
      k_grid = k_grid,
      n_MC = n_MC,
      search_method = search_method
    ),
    population_sr = pop_sr,
    est_terms = est_terms,
    sel_terms = sel_terms
  ), save = TRUE)
}
