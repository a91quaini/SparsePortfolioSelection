######################################################
#### Portfolios: Three Factor Model
######################################################

######################################################
#### Compute population mu and sigma
# We compute the mean vector (mu) and the covariance matrix (sigma) for a set
# of portfolios using the three-factor model:
#   mu = beta * mu_f
#   sigma = beta %*% sigma_f %*% t(beta) + sigma0
# where:
#   mu_f = mean of the factors
#   sigma_f = variance of the factors
#   sigma_0 = variance of the idiosyncratic errors

# Load portfolio returns and factor data.
returns <- load_data(type = "p")
ff5 <- load_data(type = "f")

# shuffle returns
set.seed(123)  # for reproducibility
shuffled_returns <- returns[, sample(ncol(returns))]

# Extract the 3ff factors from the FF5 factors
r_f <- ff5[, 1:3]
mu_f <- colMeans(r_f)
sigma_f <- var(r_f)

# For each portfolio, compute the regression (with intercept) coefficient beta
# and the mean vector.
beta <- t(solve(sigma_f, stats::cov(r_f, returns)))
mu <- beta %*% mu_f

# Compute the residuals for each portfolio:
residuals <- returns - matrix(1, nrow(returns), 1) %*% colMeans(returns) - r_f %*% t(beta)

# Compute the residual variance for each portfolio.
sigma0 <- diag(apply(residuals, 2, var))

# Construct the covariance matrix:
#   sigma = sigma_m2 * beta %*% t(beta) + sigma0_2 * I,
# where I is the identity matrix.
sigma <- beta %*% sigma_f %*% t(beta) + sigma0


######################################################
#### Simulation

# Set simulation parameters
n_sim <- 1000
n_returns <- 20

# Select returns (using the first n_returns from population parameters mu and sigma)
mu_ <- mu[1:n_returns]
sigma_ <- sigma[1:n_returns, 1:n_returns]

# Compute the population mean-variance efficient Sharpe ratio (for all assets)
mve_sr <- compute_mve_sr(mu_, sigma_, 1:n_returns)

# Create the parameter grid: different n_obs and max_card values.
# param_grid <- expand.grid(n_obs = c(10, 20, 30, 50, 100),
#                           max_card = c(5, 10, 15))
param_grid <- expand.grid(n_obs = c(10, 20, 50, 100),
                          max_card = c(5, 10))

# Choose number of cores to use: all but one (at least one)
n_cores <- max(parallel::detectCores() - 1, 1)

# Run the simulation in parallel over each parameter combination.
# Each simulation run:
#   1. Computes the sparse MVE SR (mve_sr_sparse) for the given max_card.
#   2. Runs n_sim simulations using simulate_sr_loss with the given n_obs and max_card.
#   3. Stores the resulting n_sim x 3 simulation matrix.
results_grid <- parallel::mclapply(1:nrow(param_grid), function(i) {
  n_obs_val <- param_grid$n_obs[i]
  max_card_val <- param_grid$max_card[i]

  # Compute the sparse population MVE Sharpe ratio for this max_card.
  mve_sr_sparse <- compute_sparse_mve_sr(mu_, sigma_, max_card_val)

  # Run n_sim simulations.
  sim_results <- replicate(n_sim, {
    output <- simulate_sr_loss(
      mve_sr = mve_sr_sparse$sr,
      mu = mu_,
      sigma = sigma_,
      n_obs = n_obs_val,
      max_card = max_card_val
    )
    c(sr_loss = output$sr_loss,
      sr_loss_selection = output$sr_loss_selection,
      sr_loss_estimation = output$sr_loss_estimation)
  })

  # replicate returns a 3 x n_sim matrix; transpose it to get n_sim x 3.
  sim_matrix <- t(sim_results)

  # Return a list with the parameters and the simulation matrix.
  list(n_obs = n_obs_val, max_card = max_card_val, sim_matrix = sim_matrix)

}, mc.cores = n_cores)

# Organize the results in a nested list keyed by n_obs and max_card.
results <- list()
for (res in results_grid) {
  n_obs_key <- as.character(res$n_obs)
  max_card_key <- as.character(res$max_card)
  if (is.null(results_nested[[n_obs_key]]))
    results_nested[[n_obs_key]] <- list()
  results_nested[[n_obs_key]][[max_card_key]] <- res$sim_matrix
}

output_path <- file.path("inst", "simulations", "results", "results_portfolios_three_factor_model.rds")
saveRDS(results, file = output_path)

