######################################################
#### Portfolios: One Factor Model
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

# Load portfolio returns and factor data.
returns <- load_data(type = "p")
ff5 <- load_data(type = "f")

# Extract the market factor from the FF5 factors (first column)
r_m <- ff5[, 1]
mu_m <- mean(r_m)
sigma_m2 <- var(r_m)

# For each portfolio, compute the regression (with intercept) coefficient beta
# and the mean vector.
beta <- cov(returns, r_m) / sigma_m2
mu <- mu_m * beta

# Compute the residuals for each portfolio:
residuals <- returns - matrix(1, nrow(returns), 1) %*% colMeans(returns) - r_m %*% t(beta)

# Compute the residual variance for each portfolio.
# Assume a common idiosyncratic variance: average of the portfolio residual variances.
sigma0_2 <- mean(apply(residuals, 2, var))

# Construct the covariance matrix:
#   sigma = sigma_m2 * beta %*% t(beta) + sigma0_2 * I,
# where I is the identity matrix.
sigma <- sigma_m2 * (beta %*% t(beta)) + sigma0_2 * diag(length(beta))

#### Simulation

# Compute the population mean-variance efficient Sharpe ratio
mve_sr <- compute_mve_sr(mu, sigma, 0:(ncol(returns)-1))
