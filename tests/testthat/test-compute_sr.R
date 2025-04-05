test_that("compute_sr works for valid input", {
  # Two-asset example
  weights <- c(0.5, 0.5)
  mu <- c(0.1, 0.2)
  sigma <- diag(2)
  # Expected: (0.5*0.1 + 0.5*0.2) / sqrt(0.5^2 + 0.5^2)
  expected <- (0.05 + 0.1) / sqrt(0.25 + 0.25)  # 0.15 / sqrt(0.5) = 0.15 / 0.7071068
  result <- compute_sr(weights, mu, sigma, do_checks = TRUE)
  expect_equal(result, expected, tolerance = 1e-6)
})

test_that("compute_sr errors on mismatched dimensions", {
  weights <- c(0.5, 0.5, 0.0)
  mu <- c(0.1, 0.2)
  sigma <- diag(3)
  expect_error(
    compute_sr(weights, mu, sigma, do_checks = TRUE),
    "weights and mu must be of the same length"
  )
})

test_that("compute_sr errors on non-square sigma", {
  weights <- c(0.5, 0.5)
  mu <- c(0.1, 0.2)
  sigma <- matrix(1:6, nrow = 2, ncol = 3)
  expect_error(
    compute_sr(weights, mu, sigma, do_checks = TRUE),
    "sigma must be a square matrix"
  )
})

test_that("compute_sr errors on sigma dimension mismatch", {
  weights <- c(0.5, 0.5)
  mu <- c(0.1, 0.2)
  sigma <- diag(3)  # 3x3 matrix, but weights/mu length is 2
  expect_error(
    compute_sr(weights, mu, sigma, do_checks = TRUE),
    "The dimensions of sigma must match the length of weights"
  )
})
