##################################################################
#### load_data

# --- Simulate the SparsePortfolioSelection namespace for testing ---
if (!exists("SparsePortfolioSelection")) {
  SparsePortfolioSelection <- new.env()
}

# Dummy data dimensions:
# We'll use 5 rows for all dummy matrices.
# For portfolio objects, each dummy is a 5 x 2 matrix (first column: Date, second: return).
dummy_portfolios <- list(
  returns_mebeme25    = matrix(c(rep(200805, 5), runif(5)), ncol = 2),
  returns_meop25      = matrix(c(rep(200805, 5), runif(5)), ncol = 2),
  returns_opinv25     = matrix(c(rep(200805, 5), runif(5)), ncol = 2),
  returns_bemeinv25   = matrix(c(rep(200805, 5), runif(5)), ncol = 2),
  returns_bemeop25    = matrix(c(rep(200805, 5), runif(5)), ncol = 2),
  returns_meac25      = matrix(c(rep(200805, 5), runif(5)), ncol = 2),
  returns_mebeta25    = matrix(c(rep(200805, 5), runif(5)), ncol = 2),
  returns_meinv25     = matrix(c(rep(200805, 5), runif(5)), ncol = 2),
  returns_meprior10   = matrix(c(rep(200805, 5), runif(5)), ncol = 2),
  returns_meprior122  = matrix(c(rep(200805, 5), runif(5)), ncol = 2),
  returns_meprior6013 = matrix(c(rep(200805, 5), runif(5)), ncol = 2),
  returns_mevar25     = matrix(c(rep(200805, 5), runif(5)), ncol = 2),
  returns_ind17       = matrix(c(rep(200805, 5), runif(5)), ncol = 2)
)

# Assign dummy portfolio objects to the simulated namespace
for (nm in names(dummy_portfolios)) {
  assign(nm, dummy_portfolios[[nm]], envir = SparsePortfolioSelection)
}

# Dummy CRSP_returns: a 5 x 6 matrix (first column is Date, 5 return columns)
CRSP_returns <- matrix(c(rep(200805, 5), runif(5 * 5)), ncol = 6)
assign("returns_crsp", CRSP_returns, envir = SparsePortfolioSelection)

# Dummy factors_ff5: a 5 x 7 matrix (first column is Date, 6 factor columns)
factors_ff5 <- matrix(c(rep(200805, 5), runif(5 * 6)), ncol = 7)
assign("factors_ff5", factors_ff5, envir = SparsePortfolioSelection)

# --- Begin tests ---

test_that("load_data rejects invalid type", {
  expect_error(load_data(type = "x"), "type must be either 'US' or 'International'")
})

test_that("load_data('US') returns a matrix with rows/cols > 0", {
  data <- load_data(type = "US")
  expect_true(is.matrix(data))
  expect_gt(nrow(data), 0)
  expect_gt(ncol(data), 0)
})

test_that("load_data('International') returns a matrix with rows/cols > 0", {
  data <- load_data(type = "International")
  expect_true(is.matrix(data))
  expect_gt(nrow(data), 0)
  expect_gt(ncol(data), 0)
})

# ##################################################################
# #### simulate_sr_loss
#
# test_that("simulate_sr_loss errors on invalid inputs", {
#   # Test: empty mu vector.
#   expect_error(
#     simulate_sr_loss(mve_sr = 0.8, mu = numeric(0), sigma = diag(3), n_obs = 100, max_card = 2, greedy_perc = 1.0, do_checks = TRUE),
#     "mu must be provided and be a non-empty numeric vector"
#   )
#
#   # Test: sigma is not square.
#   expect_error(
#     simulate_sr_loss(mve_sr = 0.8, mu = c(0.1, 0.2, 0.15), sigma = matrix(1:6, nrow = 2), n_obs = 100, max_card = 2, greedy_perc = 1.0, do_checks = TRUE),
#     "sigma must be a square matrix"
#   )
#
#   # Test: max_card exceeds number of assets.
#   expect_error(
#     simulate_sr_loss(mve_sr = 0.8, mu = c(0.1, 0.2, 0.15), sigma = diag(3), n_obs = 100, max_card = 4, greedy_perc = 1.0, do_checks = TRUE),
#     "max_card cannot exceed the number of assets"
#   )
#
#   # Test: negative greedy_perc.
#   expect_error(
#     simulate_sr_loss(mve_sr = 0.8, mu = c(0.1, 0.2, 0.15), sigma = diag(3), n_obs = 100, max_card = 2, greedy_perc = -0.5, do_checks = TRUE),
#     "greedy_perc must be a nonnegative numeric scalar"
#   )
#
#   # Test: mve_sr is not a finite scalar.
#   expect_error(
#     simulate_sr_loss(mve_sr = Inf, mu = c(0.1, 0.2, 0.15), sigma = diag(3), n_obs = 100, max_card = 2, greedy_perc = 1.0, do_checks = TRUE),
#     "mve_sr must be a finite numeric scalar"
#   )
# })
#
# test_that("simulate_sr_loss returns expected output structure", {
#   # Valid input for three assets.
#   mu <- c(0.1, 0.2, 0.15)
#   sigma <- diag(3)
#   mve_sr <- 0.8
#   n_obs <- 100
#   max_card <- 2
#   result <- simulate_sr_loss(mve_sr, mu, sigma, n_obs, max_card, greedy_perc = 1.0, do_checks = TRUE)
#
#   expect_true(is.list(result))
#
#   # Check that the returned list contains the expected components.
#   expected_names <- c("sr_loss", "sr_loss_selection", "sr_loss_estimation")
#   expect_equal(sort(names(result)), sort(expected_names))
# })
#
# ####################################################################
# #### calibrate_factor_model
#
# test_that("calibrate_factor_model returns proper output for valid input", {
#   set.seed(123)
#   T <- 100  # number of time periods
#   N <- 5    # number of assets
#   K <- 3    # number of factors
#
#   # Generate dummy factor returns (T x K)
#   factors <- matrix(rnorm(T * K), nrow = T, ncol = K)
#
#   # Generate dummy asset returns using a linear factor model plus noise
#   beta_true <- matrix(runif(N * K), nrow = N, ncol = K)
#   returns <- factors %*% t(beta_true) + matrix(rnorm(T * N, sd = 0.5), nrow = T, ncol = N)
#
#   # Calibrate the model with no factor weakness (weak_coeff = 0)
#   model_no_weak <- calibrate_factor_model(returns, factors, weak_coeff = 0, idiosy_vol_type = 0, do_checks = TRUE)
#
#   # Calibrate the model with moderate factor weakness (weak_coeff = 0.5)
#   model_weak <- calibrate_factor_model(returns, factors, weak_coeff = 0.5, idiosy_vol_type = 0, do_checks = TRUE)
#
#   # Check that the output is a list with components 'mu' and 'sigma'
#   expect_type(model_no_weak, "list")
#   expect_true("mu" %in% names(model_no_weak))
#   expect_true("sigma" %in% names(model_no_weak))
#
#   # Check dimensions:
#   if (is.matrix(model_no_weak$mu)) {
#     expect_equal(nrow(model_no_weak$mu), N)
#     expect_equal(ncol(model_no_weak$mu), 1)
#   } else {
#     expect_equal(length(model_no_weak$mu), N)
#   }
#   expect_true(is.matrix(model_no_weak$sigma))
#   expect_equal(nrow(model_no_weak$sigma), N)
#   expect_equal(ncol(model_no_weak$sigma), N)
#
#   # Compare outputs with different weak_coeff settings.
#   # With weak_coeff > 0, the factor betas should be scaled down,
#   # leading to a (potentially) different mu.
#   expect_false(identical(model_no_weak$mu, model_weak$mu))
#   expect_false(identical(model_no_weak$sigma, model_weak$sigma))
# })
#
# test_that("calibrate_factor_model errors with invalid input", {
#   # Test with non-matrix returns.
#   expect_error(
#     calibrate_factor_model(returns = "not a matrix", factors = matrix(1:9, 3, 3), do_checks = TRUE),
#     "returns must be provided and be a non-empty numeric matrix"
#   )
#
#   # Test with non-numeric factors.
#   expect_error(
#     calibrate_factor_model(returns = matrix(1:9, 3, 3), factors = "not a matrix", do_checks = TRUE),
#     "factors must be provided as a non-empty numeric matrix"
#   )
#
#   # Test with mismatched dimensions.
#   expect_error(
#     calibrate_factor_model(returns = matrix(1:9, 3, 3), factors = matrix(1:16, 4, 4), do_checks = TRUE),
#     "The number of rows in returns must equal the number of rows in factors"
#   )
#
#   # Test with an invalid weak_coeff.
#   expect_error(
#     calibrate_factor_model(returns = matrix(1:9, 3, 3), factors = matrix(1:9, 3, 3), weak_coeff = -0.1, do_checks = TRUE),
#     "weak_coeff must be between 0 and 1"
#   )
#   expect_error(
#     calibrate_factor_model(returns = matrix(1:9, 3, 3), factors = matrix(1:9, 3, 3), weak_coeff = 1.5, do_checks = TRUE),
#     "weak_coeff must be between 0 and 1"
#   )
#
#   # Test with an invalid idiosy_vol_type.
#   expect_error(
#     calibrate_factor_model(returns = matrix(1:9, 3, 3), factors = matrix(1:9, 3, 3), idiosy_vol_type = 2, do_checks = TRUE),
#     "idiosy_vol_type must be a numeric scalar equal to 0 or 1"
#   )
# })
