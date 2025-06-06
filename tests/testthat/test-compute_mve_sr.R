test_that("Full selection with identity sigma and ones mu works", {
  # Use full selection with mu = (1,1,1) and sigma = identity (3x3).
  mu <- rep(1, 3)
  sigma <- diag(3)
  expected <- c(sqrt(t(mu) %*% solve(sigma, mu)))
  # Explicit full selection: indices 1,2,3
  result <- compute_mve_sr(mu = mu, sigma = sigma, selection = 1:3, do_checks = TRUE)
  expect_equal(result, expected, tolerance = 1e-6)
})

test_that("Empty selection defaults to full selection", {
  # When selection is empty, the function should use the full set.
  mu <- c(1, 2, 3)
  sigma <- diag(3)
  expected <- compute_mve_sr(mu = mu, sigma = sigma, selection = 1:3, do_checks = TRUE)
  result <- compute_mve_sr(mu = mu, sigma = sigma, selection = integer(0), do_checks = TRUE)
  expect_equal(result, expected, tolerance = 1e-6)
})

test_that("Subset selection works correctly", {
  # When a subset of assets is selected, the ratio is computed on the subset.
  mu <- c(1, 2, 3, 4)
  sigma <- diag(4)
  # For selection of indices 1 and 3 (i.e., assets 1 and 3 in R's 1-indexed view),
  # the expected ratio is sqrt(1^2 + 3^2) = sqrt(10).
  expected <- sqrt(1^2 + 3^2)
  result <- compute_mve_sr(mu = mu, sigma = sigma, selection = c(1, 3), do_checks = TRUE)
  expect_equal(result, expected, tolerance = 1e-6)
})

test_that("do_checks errors on non-square sigma", {
  mu <- c(1, 2)
  sigma <- matrix(1:6, nrow = 2, ncol = 3)
  expect_error(
    compute_mve_sr(mu = mu, sigma = sigma, selection = 1:2, do_checks = TRUE),
    "sigma must be a square matrix"
  )
})

test_that("do_checks errors on dimension mismatch between mu and sigma", {
  mu <- c(1, 2, 3)
  sigma <- diag(2)
  expect_error(
    compute_mve_sr(mu = mu, sigma = sigma, selection = 1:2, do_checks = TRUE),
    "The length of mu must equal the number of rows of sigma"
  )
})

test_that("do_checks errors on out-of-bound selection", {
  mu <- c(1, 2, 3)
  sigma <- diag(3)
  # Here, index 4 is out-of-bounds (valid indices are 1,2,3).
  expect_error(
    compute_mve_sr(mu = mu, sigma = sigma, selection = c(1, 4), do_checks = TRUE),
    "Asset selection indices out of bounds"
  )
})
