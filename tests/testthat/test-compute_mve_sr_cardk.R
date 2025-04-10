test_that("Full selection returns full combination when max_card equals number of assets", {
  # For a 3-asset portfolio, with sigma = identity, the square Sharpe ratio is simply sum(mu^2)
  mu <- c(0.1, 0.2, 0.15)
  sigma <- diag(3)
  max_card <- 3
  expected_sqsr <- c(sqrt(t(mu) %*% solve(sigma, mu)))  # 0.1^2 + 0.2^2 + 0.15^2 = 0.01 + 0.04 + 0.0225 = 0.0725

  result <- compute_mve_sr_cardk(mu = mu, sigma = sigma, max_card = max_card, do_checks = TRUE)
  expect_equal(result$sr, expected_sqsr, tolerance = 1e-6)
  # Since we expect the full set to be optimal, the selection should be 0-indexed: 0,1,2
  expect_equal(result$selection, matrix(0:2, ncol = 1))
})

test_that("Subset search with max_card < full set returns the best subset", {
  # For a 3-asset portfolio with sigma = identity, evaluating subsets of size 1 and 2.
  # For one asset:
  #   best is asset 2 (0-indexed asset 1) with sqsr = 0.2^2 = 0.04.
  # For two assets:
  #   combination {1,2} (0-indexed assets 1 and 2) yields sr = sqrt(0.2^2 + 0.15^2) = 0.25.
  # Thus, with max_card = 2, the best should be {1,2} with sr = 0.25.
  mu <- c(0.1, 0.2, 0.15)
  sigma <- diag(3)
  max_card <- 2

  result <- compute_mve_sr_cardk(mu = mu, sigma = sigma, max_card = max_card, do_checks = TRUE)
  expect_equal(result$sr, 0.25, tolerance = 1e-6)
  expect_equal(result$selection, matrix(c(1, 2), ncol=1))
})

test_that("Random search (max_comb = 5) returns valid results", {
  # For a 4-asset portfolio, use max_comb = 5 so that only a fraction of combinations is evaluated.
  mu <- c(0.1, 0.2, 0.15, 0.12)
  sigma <- diag(4)
  max_card <- 3
  max_comb <- 5

  result <- compute_mve_sr_cardk(mu = mu, sigma = sigma, max_card = max_card, max_comb = max_comb, do_checks = TRUE)
  expect_true(is.numeric(result$sr))
  expect_true(length(result$selection) >= 1 && length(result$selection) <= max_card)
})

test_that("do_checks detects non-square sigma", {
  mu <- c(0.1, 0.2, 0.15)
  sigma <- matrix(1:8, nrow = 2, ncol = 4)  # non-square matrix
  expect_error(
    compute_mve_sr_cardk(mu = mu, sigma = sigma, max_card = 2, do_checks = TRUE),
    "sigma must be a square matrix"
  )
})

test_that("do_checks detects dimension mismatch between mu and sigma", {
  mu <- c(0.1, 0.2, 0.15, 0.12)
  sigma <- diag(3)
  expect_error(
    compute_mve_sr_cardk(mu = mu, sigma = sigma, max_card = 2, do_checks = TRUE),
    "The length of mu must equal the number of rows of sigma"
  )
})

# test_that("do_checks detects invalid max_card", {
#   mu <- c(0.1, 0.2, 0.15)
#   sigma <- diag(3)
#   # max_card must be between 1 and the number of assets; here 0 and 4 are invalid.
#   expect_error(
#     compute_mve_sr_cardk(mu = mu, sigma = sigma, max_card = 0, do_checks = TRUE),
#     "max_card must be between 1 and the number of assets (length of mu)"
#   )
# })
