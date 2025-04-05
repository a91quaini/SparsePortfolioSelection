library(testthat)

test_that("Full selection with identity second_moment and ones mu works", {
  mu <- rep(1, 2)
  second_moment <- diag(2)
  gamma <- 1
  # Expected: (1/gamma)*solve(second_moment, mu)
  expected <- solve(second_moment, mu) / gamma
  result <- compute_mve_weights(gamma = gamma,
                        mu = mu,
                        second_moment = second_moment,
                        selection = 0:1,
                        do_checks = TRUE)
  expect_equal(result, matrix(expected, 2, 1), tolerance = 1e-6)
})

test_that("Full selection with gamma != 1 scales correctly", {
  mu <- rep(1, 2)
  second_moment <- diag(2)
  gamma <- 2
  expected <- solve(second_moment, mu) / gamma
  result <- compute_mve_weights(gamma = gamma,
                        mu = mu,
                        second_moment = second_moment,
                        selection = 0:1,
                        do_checks = TRUE)
  expect_equal(result, matrix(expected, 2, 1), tolerance = 1e-6)
})

test_that("Empty selection defaults to full selection", {
  mu <- c(0.1, 0.2, 0.3)
  second_moment <- diag(3)
  # When selection is empty, full portfolio is used.
  expected <- compute_mve_weights(gamma = 1,
                          mu = mu,
                          second_moment = second_moment,
                          selection = 0:2,
                          do_checks = TRUE)
  result <- compute_mve_weights(gamma = 1,
                        mu = mu,
                        second_moment = second_moment,
                        selection = integer(0),
                        do_checks = TRUE)
  expect_equal(result, expected, tolerance = 1e-6)
})

test_that("Subset selection returns full-length vector with zeros for unselected assets", {
  mu <- c(0.1, 0.2, 0.3)
  second_moment <- diag(3)
  # Select assets 0 and 2 (0-indexed)
  selection <- c(0, 2)
  result <- compute_mve_weights(gamma = 1,
                        mu = mu,
                        second_moment = second_moment,
                        selection = selection,
                        do_checks = TRUE)

  expect_equal(length(result), 3)
  expect_equal(result[2], 0)  # asset 1 (second element) should be zero

  # Compute expected solution for the selected subset:
  first_sel <- mu[selection + 1]  # convert 0-index to 1-index for R
  second_sel <- second_moment[selection + 1, selection + 1]
  expected_sel <- as.vector(solve(second_sel, first_sel))
  expected_full <- rep(0, 3)
  expected_full[selection + 1] <- expected_sel
  expect_equal(result, matrix(expected_full, 3, 1), tolerance = 1e-6)
})

test_that("do_checks detects non-square second_moment", {
  mu <- c(0.1, 0.2)
  second_moment <- matrix(1:6, nrow = 2, ncol = 3)
  expect_error(
    compute_mve_weights(gamma = 1,
                mu = mu,
                second_moment = second_moment,
                selection = 0:1,
                do_checks = TRUE),
    "second_moment must be a square matrix"
  )
})

test_that("do_checks detects mu length not equal to second_moment dimensions", {
  mu <- c(0.1, 0.2, 0.3)
  second_moment <- diag(2)
  expect_error(
    compute_mve_weights(gamma = 1,
                mu = mu,
                second_moment = second_moment,
                selection = 0:1,
                do_checks = TRUE),
    "The length of mu must equal the number of rows of second_moment"
  )
})

test_that("do_checks detects asset selection index out of bounds", {
  mu <- c(0.1, 0.2, 0.3)
  second_moment <- diag(3)
  # Here, index 3 is out-of-bounds (valid indices are 0, 1, 2).
  selection <- c(0, 3)
  expect_error(
    compute_mve_weights(gamma = 1,
                mu = mu,
                second_moment = second_moment,
                selection = selection,
                do_checks = TRUE),
    "Asset selection indices out of bounds"
  )
})

test_that("gamma = 0 produces Inf or NaN values", {
  mu <- rep(1, 2)
  second_moment <- diag(2)
  # With gamma = 0, division by zero is expected.
  result <- compute_mve_weights(gamma = 0,
                        mu = mu,
                        second_moment = second_moment,
                        selection = 0:1,
                        do_checks = TRUE)
  expect_true(all(is.infinite(result)) || any(is.nan(result)))
})
