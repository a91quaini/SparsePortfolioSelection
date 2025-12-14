# Monte Carlo analysis for FF3-style simulations (HPC-safe: LASSO only, no MIQP)

library(SparsePortfolioSelection)

set.seed(123)
Nn <- 6L
n_cores <- max(1L, Nn - 1L)

# 1) Simulation parameters
n_assets <- 100              # [25, 50, 100]
n_obs_grid <- c(120L, 480L)  # allow multiple T (population SR does not depend on T)
k_grid <- 2:(n_assets - 1L)
n_MC <- 200
search_method <- "lasso"     # fixed to lasso to avoid MIQP on HPC

# 2) LASSO search parameters (update n_obs inside loop)
lasso_base <- list(
  nlambda = 300L,
  lambda_min_ratio = 1e-3,
  nadd = 80L,
  nnested = 2L,
  alpha = 1,
  standardize = FALSE,
  compute_weights = TRUE,
  normalize_weights = FALSE,
  use_refit = FALSE
)

# 3) Calibrate population parameters (fixed constants + one draw of alphas/betas)
params <- calibrate_ff3(n_assets)
mu_pop <- params$mu_true
sigma_pop <- params$sigma_true

# 4) Population SR via LASSO only (uses first n_obs value for design; SR does not change with n_obs)
lasso_pop <- modifyList(lasso_base, list(
  n_obs = n_obs_grid[[1]],
  lambda_min_ratio = 1e-4,  # denser path for population calc
  nadd = 50L
))
pop_sr <- numeric(length(k_grid))
for (i in seq_along(k_grid)) {
  k <- k_grid[i]
  res <- mve_lasso_search(
    mu = mu_pop,
    sigma = sigma_pop,
    k = k,
    n_obs = lasso_pop$n_obs,
    nlambda = lasso_pop$nlambda,
    lambda_min_ratio = lasso_pop$lambda_min_ratio,
    alpha = lasso_pop$alpha,
    nadd = lasso_pop$nadd,
    nnested = lasso_pop$nnested,
    standardize = lasso_pop$standardize,
    compute_weights = TRUE,
    normalize_weights = FALSE,
    use_refit = lasso_pop$use_refit
  )
  pop_sr[i] <- as.numeric(res$sr)
}

dir.create("inst/simulations/results", recursive = TRUE, showWarnings = FALSE)

for (n_obs in n_obs_grid) {
  lasso_params <- modifyList(lasso_base, list(n_obs = n_obs))
  mve_search_fn <- mve_lasso_search
  mve_search_fn_params <- lasso_params

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

  fname <- sprintf("inst/simulations/results/results_ff3_lasso_N%d_T%d_hpc_mc_analysis_ff3.rds",
                   n_assets, n_obs)
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
