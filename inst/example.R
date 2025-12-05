# Simple example demonstrating MIQP and LASSO searches on simulated data

library(SparsePortfolioSelection)

set.seed(123)

# Simulated means/covariance for 5 assets; make two assets clearly better
mu <- c(0.15, 0.14, 0.05, 0.04, 0.03)  # assets 1 & 2 have higher expected return
A <- matrix(rnorm(25), 5, 5)
base_sigma <- crossprod(A) + 0.05 * diag(5)  # SPD baseline
# Lower variance for the first two assets to boost their SR
scale_vec <- c(0.5, 0.6, 1.0, 1.0, 1.0)
sigma <- diag(scale_vec) %*% base_sigma %*% diag(scale_vec)

# MIQP heuristic (requires Gurobi)
# Set GUROBI_HOME and load the package with Gurobi installed before running.
miqp_res <- mve_miqp_search(
  mu, sigma,
  k = 2,
  exactly_k = TRUE,
  gamma = 1.0,
  fmin = 0,
  fmax = 0.8,
  mipgap = 1e-4,
  time_limit = 10,
  threads = 0,
  compute_weights = TRUE,
  normalize_weights = TRUE,
  verbose = FALSE
)
print(miqp_res)

# LASSO heuristic (glmnet backend)
lasso_res <- mve_lasso_search(
  mu = mu,
  sigma = sigma,
  n_obs = 200,
  k = 2,
  nlambda = 50,
  alpha = 1,
  standardize = FALSE,
  compute_weights = TRUE,
  use_refit = TRUE
)
print(lasso_res)

# LASSO heuristic (glmnet backend)
lasso_res <- mve_lasso_search(
  mu = mu,
  sigma = sigma,
  n_obs = 200,
  k = 2,
  nlambda = 50,
  alpha = c(0.5, 0.6, 0.7, 0.8, 0.9),
  n_folds = 3,
  standardize = TRUE,
  compute_weights = TRUE,
  use_refit = TRUE
)
print(lasso_res)

# Exhaustive search on population moments
exhaustive_res <- mve_exhaustive_search(
  mu,
  sigma,
  k = 2,
  stabilize_sigma = TRUE,
  do_checks = TRUE,
  enumerate_all = TRUE,
  max_samples = 0L,
  dedup_samples = TRUE,
  compute_weights = TRUE
)
print(exhaustive_res)
