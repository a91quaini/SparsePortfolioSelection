# Simulation analysis for FF3-style simulations (LARS only)

library(SparsePortfolioSelection)

set.seed(123)
Nn <- 6L
n_cores <- max(1L, Nn - 1L)

# 1) Simulation parameters
n_assets <- 100
n_obs_grid <- c(120L, 240L, 480L, 960L)
k_grid <- 2:(n_assets - 1L)
n_MC <- 10000

# 2) Calibrate population parameters (fixed constants + one draw of alphas/betas)
params <- calibrate_ff3(n_assets)
mu_pop <- params$mu_true
sigma_pop <- params$sigma_true
opt_pop_sr <- sqrt(sum(mu_pop * solve(sigma_pop, mu_pop)))

# 3) Population SR for each k (uses first n_obs value; SR does not change with n_obs)
miqp_params <- list(
  exactly_k = TRUE,
  gamma = 1.0,
  fmin = -0.25,
  fmax = 0.25,
  expand_rounds = 20L,
  expand_factor = 3.0,
  expand_tol = 1e-2,
  mipgap = 1e-4,
  time_limit = 150,
  threads = 0,
  normalize_weights = FALSE,
  verbose = FALSE,
  ridge_epsilon = 0.0
)

# 3) Population SR for each k (uses first n_obs value; SR does not change with n_obs)
pop_sr <- vapply(k_grid, function(k) {
  sr_lars <- tryCatch({
    res <- mve_lars_search(
      mu = mu_pop,
      sigma = sigma_pop,
      n_obs = n_obs_grid[[1]],
      k = k,
      ridge_epsilon = 0.0,
      tol_nnl = 1e-10,
      normalize_weights = FALSE,
      normalization_type = 1L,
      use_refit = FALSE,
      compute_sr = TRUE,
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
      expand_rounds = miqp_params$expand_rounds,
      expand_factor = miqp_params$expand_factor,
      expand_tol = miqp_params$expand_tol,
      mipgap = miqp_params$mipgap,
      time_limit = miqp_params$time_limit,
      threads = miqp_params$threads,
      normalize_weights = miqp_params$normalize_weights,
      use_refit = FALSE,
      verbose = miqp_params$verbose,
      ridge_epsilon = miqp_params$ridge_epsilon,
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
    ridge_epsilon = 0.0,
    tol_nnl = 1e-10,
    normalize_weights = FALSE,
    normalization_type = 1L,
    use_refit = FALSE,
    compute_sr = FALSE,
    do_checks = FALSE
  )

  # Monte Carlo runs: simulate returns, compute decomposition terms
  run_one <- function(mc_id) {
    run_one_ff3_mc(
      mc_id = mc_id,
      n_obs = n_obs,
      params = params,
      mu_pop = mu_pop,
      sigma_pop = sigma_pop,
      k_grid = k_grid,
      mve_search_fn = mve_search_fn,
      mve_search_fn_params = mve_search_fn_params,
      compute_span = TRUE,
      span_sys = "factors",   # or "pca"
      span_m = 3L
    )
  }

  res_list <- parallel::mclapply(seq_len(n_MC), run_one, mc.cores = n_cores)

  est_terms  <- t(vapply(res_list, `[[`, numeric(length(k_grid)), "est"))
  sel_terms  <- t(vapply(res_list, `[[`, numeric(length(k_grid)), "sel"))
  span_terms <- t(vapply(res_list, `[[`, numeric(length(k_grid)), "span"))

  colnames(est_terms)  <- paste0("k", k_grid)
  colnames(sel_terms)  <- paste0("k", k_grid)
  colnames(span_terms) <- paste0("k", k_grid)

  fname <- sprintf("inst/simulations/results/results_ff3_lars_N%d_T%d_sim_analysis_ff3.rds",
                   n_assets, n_obs)
  # saveRDS(list(
  #   config = list(
  #     n_assets = n_assets,
  #     n_obs = n_obs,
  #     k_grid = k_grid,
  #     n_MC = n_MC,
  #     search_method = "lars",
  #     mve_search_fn_params = mve_search_fn_params,
  #     span_sys = "pca",
  #     span_m = 3L
  #   ),
  #   population_sr = pop_sr,
  #   opt_pop_sr = opt_pop_sr,
  #   est_terms = est_terms,
  #   sel_terms = sel_terms,
  #   span_terms = span_terms
  # ), file = fname)
  #
  # message(sprintf("Saved results to %s", fname))

  plot_mve_sr_decomposition(list(
    config = list(
      n_assets = n_assets,
      n_obs = n_obs,
      k_grid = k_grid,
      n_MC = n_MC,
      search_method = "lars"
    ),
    population_sr = pop_sr,
    opt_pop_sr = opt_pop_sr,
    est_terms = est_terms,
    sel_terms = sel_terms,
    span_terms = span_terms
  ), save = TRUE)
}
