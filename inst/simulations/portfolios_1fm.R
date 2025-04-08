######################################################
#### Portfolios: One Factor Model
######################################################

######################################################
#### Compute population mu and sigma
# We compute the mean vector (mu) and the covariance matrix (sigma) for a set
# of portfolios using the one-factor model:
#   mu = mu_m * beta
#   sigma = sigma_m^2 * beta %*% t(beta) + sigma_0^2 * I
# where:
#   mu_m = mean(market returns)
#   sigma_m^2 = variance(market returns)
#   sigma_0^2 = average residual variance across portfolios

source("inst/simulations/utils.R")

# Load portfolio returns and factor data.
returns <- load_data(type = "p")
ff5 <- load_data(type = "f")

# shuffle returns
set.seed(123)  # for reproducibility
returns <- returns[, sample(ncol(returns))]

# Select returns
n_returns <- 20
returns <- returns[, 1:n_returns]

# Calibrate the model
model <- calibrate_factor_model(returns, ff5[, 1, drop=FALSE], idiosy_vol_type = 0)
mu <- model$mu
sigma <- model$sigma


######################################################
#### Simulation

# Compute the population mean-variance efficient Sharpe ratio (for all assets)
mve_sr <- compute_mve_sr(mu, sigma, 1:n_returns)

# Set simulation parameters
n_sim <- 1000
# Set the various sample sizes
n_obs <- c(10, 15, 20, 50, 100)
# Set the various cardinality constraints
max_card <- c(5, 10, 15)

# Compute the population mean-variance efficient Sharpe ratio with cardinality constraints
mve_sr_sparse <- list()
for (card in max_card) {
  # Compute the sparse MVE Sharpe ratio for the current cardinality
  mve_sr_sparse[[as.character(card)]] <- compute_sparse_mve_sr(mu, sigma, card)$sr
}

################################################################################
##### LONG COMPUTATION TIME, YOU CAN SKIP THIS #################################
################################################################################
#### Compute simulations
#### This part of the code takes time. Instead of running it you can source
#### the results from the file "results_portfolios_1fm.rds" in the "inst/simulations/results" folder
#### This is done below at ####*

# Create the parameter grid: different n_obs and max_card values.
# param_grid <- expand.grid(n_obs = c(10, 20, 30, 50, 100),
#                           max_card = c(5, 10, 15))
param_grid <- expand.grid(n_obs = n_obs,
                          max_card = max_card)

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
  mve_sr_sparse <- compute_sparse_mve_sr(mu, sigma, max_card_val)

  # Run n_sim simulations.
  sim_results <- replicate(n_sim, {
    output <- simulate_sr_loss(
      mve_sr = mve_sr_sparse$sr,
      mu = mu,
      sigma = sigma,
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
  if (is.null(results[[n_obs_key]]))
    results[[n_obs_key]] <- list()
  results[[n_obs_key]][[max_card_key]] <- res$sim_matrix
}

# Save the results to a file.
output_path <- file.path("inst", "simulations", "results", "results_portfolios_1fm_n20.rds")
saveRDS(results, file = output_path)

################################################################################
#### ----> End of the computation part
################################################################################


####* Load the results from the file
evaluate_simulation_results(f_name = "portfolios_1fm_n20",  N = n_returns, mve_sr = mve_sr, mve_sr_sparse = mve_sr_sparse)




######################################################
### Weak factor
model <- calibrate_factor_model(returns, ff5[, 1, drop=FALSE], weak_coeff = 0.5, idiosy_vol_type = 0)
mu <- model$mu
sigma <- model$sigma

######################################################
#### Simulation

# Set simulation parameters
n_sim <- 1000

# Compute the population mean-variance efficient Sharpe ratio (for all assets)
mve_sr <- compute_mve_sr(mu, sigma, 1:n_returns)

# Create the parameter grid: different n_obs and max_card values.
# param_grid <- expand.grid(n_obs = c(10, 20, 30, 50, 100),
#                           max_card = c(5, 10, 15))
param_grid <- expand.grid(n_obs = c(10, 15, 20, 50, 100),
                          max_card = c(5, 10, 15))

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
  mve_sr_sparse <- compute_sparse_mve_sr(mu, sigma, max_card_val)

  # Run n_sim simulations.
  sim_results <- replicate(n_sim, {
    output <- simulate_sr_loss(
      mve_sr = mve_sr_sparse$sr,
      mu = mu,
      sigma = sigma,
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
  if (is.null(results[[n_obs_key]]))
    results[[n_obs_key]] <- list()
  results[[n_obs_key]][[max_card_key]] <- res$sim_matrix
}

# Save the results to a file.
output_path <- file.path("inst", "simulations", "results", "results_portfolios_1fm_n20_weak05.rds")
saveRDS(results, file = output_path)
