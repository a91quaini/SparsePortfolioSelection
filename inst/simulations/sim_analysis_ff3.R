# Simulation analysis for FF3-style simulations (LARS only)

library(SparsePortfolioSelection)

set.seed(123)
Nn <- 6L
n_cores <- max(1L, Nn - 1L)

# 1) Simulation parameters
n_assets <- 100
n_obs_grid <- c(120L, 240L, 480L, 960L)
k_grid <- 2:(n_assets - 1L)
n_MC <- 1000

# 2) Calibrate population parameters (fixed constants + one draw of alphas/betas)
params <- calibrate_ff3(n_assets)
mu_pop <- params$mu_true
sigma_pop <- params$sigma_true

# 3) Population SR for each k (uses first n_obs value; SR does not change with n_obs)
miqp_params <- list(
  exactly_k = FALSE,
  gamma = 1.0,
  fmin = -0.25,
  fmax = 0.25,
  mipgap = 1e-4,
  time_limit = 150,
  threads = 0,
  compute_weights = TRUE,
  normalize_weights = FALSE,
  verbose = FALSE
)

# 3) Population SR for each k (uses first n_obs value; SR does not change with n_obs)
pop_sr <- vapply(k_grid, function(k) {
  sr_lars <- tryCatch({
    res <- mve_lars_search(
      mu = mu_pop,
      sigma = sigma_pop,
      n_obs = n_obs_grid[[1]],
      k = k,
      stabilize_sigma = FALSE,
      tol_nnl = 1e-10,
      compute_weights = TRUE,
      normalize_weights = FALSE,
      use_refit = FALSE,
      do_checks = FALSE
    )
    as.numeric(res$sr)
  }, error = function(e) -Inf)

  sr_miqp <- tryCatch({
    res <- mve_miqp_search(
      mu_pop, sigma_pop, k,
      exactly_k = miqp_params$exactly_k,
      m = 1L,
      gamma = miqp_params$gamma,
      fmin = miqp_params$fmin,
      fmax = miqp_params$fmax,
      expand_rounds = 25L,
      expand_factor = 2.0,
      expand_tol = 1e-4,
      mipgap = miqp_params$mipgap,
      time_limit = miqp_params$time_limit,
      threads = miqp_params$threads,
      compute_weights = miqp_params$compute_weights,
      normalize_weights = miqp_params$normalize_weights,
      use_refit = FALSE,
      verbose = miqp_params$verbose,
      epsilon = 1e-6,
      stabilize_sigma = FALSE,
      do_checks = FALSE
    )
    as.numeric(res$sr)
  }, error = function(e) -Inf)

  max(sr_lars, sr_miqp, na.rm = TRUE)
}, numeric(1))

dir.create("inst/simulations/results", recursive = TRUE, showWarnings = FALSE)

for (n_obs in n_obs_grid) {
  mve_search_fn <- mve_lars_search
  mve_search_fn_params <- list(
    n_obs = n_obs,
    stabilize_sigma = FALSE,
    tol_nnl = 1e-10,
    compute_weights = TRUE,
    normalize_weights = FALSE,
    use_refit = FALSE,
    do_checks = FALSE
  )

  # Monte Carlo runs: simulate returns, compute decomposition terms
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

  res_list <- parallel::mclapply(seq_len(n_MC), run_one, mc.cores = n_cores)

  est_terms <- t(vapply(res_list, `[[`, numeric(length(k_grid)), "est"))
  sel_terms <- t(vapply(res_list, `[[`, numeric(length(k_grid)), "sel"))
  colnames(est_terms) <- paste0("k", k_grid)
  colnames(sel_terms) <- paste0("k", k_grid)

  fname <- sprintf("inst/simulations/results/results_ff3_lars_N%d_T%d_sim_analysis_ff3.rds",
                   n_assets, n_obs)
  saveRDS(list(
    config = list(
      n_assets = n_assets,
      n_obs = n_obs,
      k_grid = k_grid,
      n_MC = n_MC,
      search_method = "lars",
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
      search_method = "lars"
    ),
    population_sr = pop_sr,
    est_terms = est_terms,
    sel_terms = sel_terms
  ), save = TRUE)
}
