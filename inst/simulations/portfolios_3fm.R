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

source("inst/simulations/utils.R")

# Load portfolio returns and factor data.
returns <- load_data(type = "p")
ff5 <- load_data(type = "f")

# shuffle returns
set.seed(123)  # for reproducibility
returns <- returns[, sample(ncol(returns))]

# Select returns and factors
n_returns <- 20
returns <- returns[, 1:n_returns]
factors <- ff5[, 3, drop=FALSE]

# Calibrate the model
model <- calibrate_factor_model(returns, factors, idiosy_vol_type = 1)
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
#### the results from the file "results_portfolios_3fm_n20.rds" in the "inst/simulations/results" folder
#### This is done below at ####*

compute_simulation_results(n_obs = n_obs,
                           max_card = max_card,,
                           n_sim = n_sim,
                           mve_sr_sparse = mve_sr_sparse,
                           mu = mu,
                           sigma = sigma,
                           simulate_sr_loss = simulate_sr_loss,
                           seed = 123,
                           file_name = "results_portfolios_3fm_n20.rds")

################################################################################
#### ----> End of the computation part



####* Evaluate the results
results <- evaluate_simulation_results(f_name = "portfolios_3fm_n20",
                                       N = n_returns,
                                       mve_sr = mve_sr,
                                       mve_sr_sparse = mve_sr_sparse)











######################################################
######################################################
######################################################
### Weak factor
######################################################

model <- calibrate_factor_model(returns, factors, weak_coeff = 0.5, idiosy_vol_type = 1)
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

compute_simulation_results(n_obs = n_obs,
                           max_card = max_card,,
                           n_sim = n_sim,
                           mve_sr_sparse = mve_sr_sparse,
                           mu = mu,
                           sigma = sigma,
                           simulate_sr_loss = simulate_sr_loss,
                           seed = 123,
                           file_name = "results_portfolios_3fm_n20_weak05.rds")

################################################################################
#### ----> End of the computation part
################################################################################


####* Evaluate the results
results <- evaluate_simulation_results(f_name = "portfolios_3fm_n20_weak05",
                                       N = n_returns,
                                       mve_sr = mve_sr,
                                       mve_sr_sparse = mve_sr_sparse)

