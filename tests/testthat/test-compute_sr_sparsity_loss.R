library(testthat)

test_that("compute_sr_sparsity_loss returns a list with expected elements", {
  # Create sample population parameters.
  mu <- c(0.1, 0.2, 0.15)
  sigma <- diag(3)

  # Create sample parameters (slightly perturbed).
  mu_hat <- c(0.12, 0.18, 0.16)
  sigma_hat <- diag(3) * 1.1

  mve_sr <- 0.8   # Example population MVE Sharpe ratio.
  max_card <- 2
  greedy_perc <- 1.0

  result <- compute_sr_sparsity_loss(mve_sr, mu, sigma, mu_hat, sigma_hat, max_card, greedy_perc, do_checks = TRUE)

  expect_true(is.list(result))
  expect_true(all(c("sr_loss", "sr_loss_selection", "sr_loss_estimation") %in% names(result)))
  expect_true(is.numeric(result$sr_loss))
  expect_true(is.numeric(result$sr_loss_selection))
  expect_true(is.numeric(result$sr_loss_estimation))
})

test_that("compute_sr_sparsity_loss errors with mismatched sample dimensions", {
  mu <- c(0.1, 0.2, 0.15)
  sigma <- diag(3)
  mu_hat <- c(0.1, 0.2)  # Wrong length
  sigma_hat <- diag(2)
  mve_sr <- 0.8
  max_card <- 2
  greedy_perc <- 1.0

  expect_error(
    compute_sr_sparsity_loss(mve_sr, mu, sigma, mu_hat, sigma_hat, max_card, greedy_perc, do_checks = TRUE)
  )
})

test_that("compute_sr_sparsity_loss errors with invalid max_card", {
  mu <- c(0.1, 0.2, 0.15)
  sigma <- diag(3)
  mu_hat <- c(0.12, 0.18, 0.16)
  sigma_hat <- diag(3) * 1.1
  mve_sr <- 0.8
  max_card <- 4  # Invalid: more than length(mu_hat)
  greedy_perc <- 1.0

  expect_error(
    compute_sr_sparsity_loss(mve_sr, mu, sigma, mu_hat, sigma_hat, max_card, greedy_perc, do_checks = TRUE),
    "max_card cannot exceed the number of assets"
  )
})
