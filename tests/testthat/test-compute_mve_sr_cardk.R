test_that("Full selection returns full combination when max_card equals number of assets", {
  mu <- c(0.1, 0.2, 0.15)
  sigma <- diag(3)
  max_card <- 3
  gamma <- 1

  # expected Sharpe ratio: sqrt(mu^T Sigma^{-1} mu)
  expected_sr <- sqrt(as.numeric(t(mu) %*% solve(sigma, mu)))
  # expected weights = (1/gamma) * (Sigma)^{-1} * mu = mu / gamma
  expected_weights <- solve(sigma, mu) / gamma

  result <- compute_mve_sr_cardk(
    mu         = mu,
    sigma      = sigma,
    max_card   = max_card,
    gamma      = gamma,
    do_checks  = TRUE
  )

  expect_equal(result$sr, expected_sr,       tolerance = 1e-6)
  expect_equal(as.numeric(result$weights), expected_weights, tolerance = 1e-6)
  expect_equal(as.integer(result$selection), 0:2)
})

test_that("Subset search with max_card < full set returns the best subset and weights", {
  mu <- c(0.1, 0.2, 0.15)
  sigma <- diag(3)
  max_card <- 2
  gamma <- 1

  # best subset {1,2} (0-indexed) gives sr = sqrt(0.2^2 + 0.15^2) = 0.25
  expected_sr      <- sqrt(mu[2]^2 + mu[3]^2)
  expected_weights <- c(0, solve(sigma[c(2,3),c(2,3)], mu[c(2,3)]) / gamma)

  result <- compute_mve_sr_cardk(
    mu         = mu,
    sigma      = sigma,
    max_card   = max_card,
    gamma      = gamma,
    do_checks  = TRUE
  )

  expect_equal(result$sr, expected_sr,       tolerance = 1e-6)
  expect_equal(as.integer(result$selection), c(1, 2))
  expect_equal(as.numeric(result$weights),   expected_weights, tolerance = 1e-6)
})

test_that("Gamma parameter scales weights and leaves sr unchanged", {
  mu <- c(0.05, 0.15, 0.10)
  sigma <- diag(3)
  max_card <- 2

  # compute with gamma = 1
  r1 <- compute_mve_sr_cardk(mu, sigma, max_card = max_card, gamma = 1, do_checks = FALSE)
  # compute with gamma = 5
  r2 <- compute_mve_sr_cardk(mu, sigma, max_card = max_card, gamma = 5, do_checks = FALSE)

  expect_equal(r2$sr, r1$sr, tolerance = 1e-8)
  expect_equal(r2$weights, r1$weights / 5, tolerance = 1e-8)
})

test_that("Random search (max_comb = 5) returns valid sr and weights", {
  mu <- c(0.1, 0.2, 0.15, 0.12)
  sigma <- diag(4)
  max_card <- 3
  max_comb <- 5
  gamma    <- 2

  result <- compute_mve_sr_cardk(
    mu        = mu,
    sigma     = sigma,
    max_card  = max_card,
    max_comb  = max_comb,
    gamma     = gamma,
    do_checks = TRUE
  )

  expect_true(is.numeric(result$sr))
  expect_true(is.numeric(result$weights))
  expect_true(length(result$selection) >= 1 && length(result$selection) <= max_card)
  expect_equal(length(result$weights), 4)
})

test_that("do_checks detects non-square sigma", {
  mu <- c(0.1, 0.2, 0.15)
  sigma <- matrix(1:8, nrow = 2, ncol = 4)  # non-square
  expect_error(
    compute_mve_sr_cardk(mu, sigma, max_card = 2, gamma = 1, do_checks = TRUE),
    "sigma must be a square matrix"
  )
})

test_that("do_checks detects dimension mismatch between mu and sigma", {
  mu <- c(0.1, 0.2, 0.15, 0.12)
  sigma <- diag(3)
  expect_error(
    compute_mve_sr_cardk(mu, sigma, max_card = 2, gamma = 1, do_checks = TRUE),
    "The length of mu must equal the number of rows of sigma"
  )
})

test_that("do_checks detects non-numeric gamma", {
  mu <- c(0.1, 0.2, 0.15)
  sigma <- diag(3)
  expect_error(
    compute_mve_sr_cardk(mu, sigma, max_card = 2, gamma = "high", do_checks = TRUE),
    "gamma must be numeric"
  )
})

