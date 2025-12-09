# Monte Carlo analysis for FF3-style simulations (mirrors MATLAB scripts)

library(SparsePortfolioSelection)

# Ensure Gurobi env vars are set when needed
if (Sys.getenv("GUROBI_HOME") == "") {
  Sys.setenv(GUROBI_HOME = "/Library/gurobi1201/macos_universal2")
}
if (Sys.getenv("DYLD_LIBRARY_PATH") == "" && Sys.getenv("GUROBI_HOME") != "") {
  Sys.setenv(DYLD_LIBRARY_PATH = file.path(Sys.getenv("GUROBI_HOME"), "lib"))
}

set.seed(123)

# 1) Simulation parameters
n_assets <- 100              # [25, 50, 100]
n_obs <- 480                 # [120 240 480]
k_grid <- 2:(n_assets-1)
n_MC <- 200                  # 200
search_method <- "miqp"     # set to "lasso" or "miqp"

# 2) Search function parameters
lasso_params <- list(
  n_obs = n_obs,
  nlambda = 100L,
  lambda_min_ratio = 1e-3,
  nadd = 80L,
  nnested = 2L,
  alpha = 1,
  standardize = FALSE,
  compute_weights = TRUE,
  normalize_weights = FALSE,
  use_refit = FALSE
)

miqp_params <- list(
  exactly_k = TRUE,
  gamma = 1.0,
  fmin = 0,
  fmax = 0.25,
  mipgap = 1e-4,
  time_limit = 60,
  threads = 0,
  compute_weights = TRUE,
  normalize_weights = FALSE,
  verbose = FALSE
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

# 3) Calibrate population parameters (fixed constants + one draw of alphas/betas)
params <- calibrate_ff3(n_assets)
mu_pop <- params$mu_true
sigma_pop <- params$sigma_true

# 4) Population SR for each k (use boosted MIQP/LASSO + monotone)
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

# 5) Monte Carlo runs: simulate returns, compute decomposition terms
est_terms <- matrix(NA_real_, nrow = n_MC, ncol = length(k_grid))
sel_terms <- matrix(NA_real_, nrow = n_MC, ncol = length(k_grid))
colnames(est_terms) <- paste0("k", k_grid)
colnames(sel_terms) <- paste0("k", k_grid)

for (mc in seq_len(n_MC)) {
  message(sprintf("MC run %d / %d", mc, n_MC))
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

# Save results
dir.create("inst/simulations/results", recursive = TRUE, showWarnings = FALSE)
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

# Plot and save figures
plots <- plot_mve_sr_decomposition(list(
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
