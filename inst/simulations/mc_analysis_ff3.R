# Monte Carlo analysis for FF3-style simulations (mirrors MATLAB scripts)

library(SparsePortfolioSelection)

set.seed(123)
Nn <- 6L
n_cores <- max(1L, Nn - 1L)

# 1) Simulation parameters
n_assets <- 100              # [25, 50, 100]
n_obs_grid <- c(120L, 480L)  # allow multiple T (population SR does not depend on T)
k_grid <- 2:(n_assets-1)
n_MC <- 1000                  # 200
search_method <- "lasso"     # set to "lasso" or "miqp"

# 2) Search function parameters (static across n_obs unless n_obs is needed)
miqp_params <- list(
  exactly_k = TRUE,
  gamma = 1.0,
  fmin = 0,
  fmax = 0.25,
  mipgap = 1e-4,
  time_limit = 100,
  threads = 0,
  compute_weights = TRUE,
  normalize_weights = FALSE,
  verbose = FALSE
)

# 3) Calibrate population parameters (fixed constants + one draw of alphas/betas)
params <- calibrate_ff3(n_assets)
mu_pop <- params$mu_true
sigma_pop <- params$sigma_true

# 4) Population SR for each k (uses first n_obs value; SR does not change with n_obs)
pop_sr <- compute_population_mve_sr(
  mu_pop = mu_pop,
  sigma_pop = sigma_pop,
  k_grid = k_grid,
  n_obs = n_obs_grid[[1]],
  lasso_nlambda = 300L,
  lasso_lambda_min_ratio = 1e-4,
  miqp_time_limit = 200,
  miqp_mipgap = 1e-5,
  verbose = FALSE
)

dir.create("inst/simulations/results", recursive = TRUE, showWarnings = FALSE)

for (n_obs in n_obs_grid) {
  lasso_params <- list(
    n_obs = n_obs,
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

  if (search_method == "lasso") {
    mve_search_fn <- mve_lasso_search
    mve_search_fn_params <- lasso_params
  } else if (search_method == "miqp") {
    mve_search_fn <- mve_miqp_search
    mve_search_fn_params <- miqp_params
  } else {
    stop("search_method must be 'lasso' or 'miqp'")
  }

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

  fname <- sprintf("inst/simulations/results/results_ff3_%s_N%d_T%d_mc_analysis_ff3.rds",
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
