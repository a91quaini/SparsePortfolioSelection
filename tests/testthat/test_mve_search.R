library(testthat)
library(SparsePortfolioSelection)

TOL <- 1e-7

skip_if_no_gurobi <- function() {
  if (Sys.getenv("GUROBI_HOME") == "") skip("Gurobi not available")
  if (!identical(Sys.getenv("RUN_GUROBI_TESTS"), "1")) skip("Gurobi tests disabled (set RUN_GUROBI_TESTS=1 to run)")
}

test_that("exhaustive search returns sane outputs", {
  mu <- c(0.1, 0.2, 0.15, 0.05)
  A  <- matrix(c(0.2,0.01,0.00,0.00,
                 0.01,0.25,0.01,0.00,
                 0.00,0.01,0.18,0.02,
                 0.00,0.00,0.02,0.22), nrow = 4, byrow = TRUE)
  sigma <- A
  res <- mve_exhaustive_search(mu, sigma, k = 2, do_checks = TRUE)
  expect_true(is.numeric(res$sr) && length(res$sr) == 1)
  expect_equal(length(res$selection), 2)
  expect_equal(length(res$weights), length(mu))
  expect_true(res$status %in% c("EXHAUSTIVE", "SAMPLED"))
})

test_that("exhaustive search k=full equals mve", {
  mu <- c(0.1, 0.2, 0.15)
  sigma <- diag(c(0.2, 0.25, 0.18))
  res <- mve_exhaustive_search(mu, sigma, k = length(mu))
  sr_mve <- compute_mve_sr(mu, sigma)
  expect_lt(abs(res$sr - sr_mve), 1e-10)
  expect_equal(sort(res$selection), 1:(length(mu)))
})

test_that("exhaustive search sampling path", {
  mu <- c(0.1, 0.2, 0.15, 0.05, 0.12)
  sigma <- diag(5)
  res <- mve_exhaustive_search(mu, sigma, k = 3,
                               enumerate_all = FALSE, max_samples = 5,
                               dedup_samples = TRUE, compute_weights = TRUE)
  expect_equal(length(res$selection), 3)
  expect_equal(length(res$weights), length(mu))
})

test_that("exhaustive search respects compute_weights=FALSE", {
  mu <- c(0.05, 0.1, 0.02)
  sigma <- diag(3)
  res <- mve_exhaustive_search(mu, sigma, k = 2, compute_weights = FALSE)
  expect_true(all(res$weights == 0))
})

test_that("exhaustive search dedup vs no dedup", {
  mu <- c(0.1, 0.2, 0.15, 0.12)
  sigma <- diag(4)
  res1 <- mve_exhaustive_search(mu, sigma, k = 2,
                                enumerate_all = FALSE, max_samples = 3,
                                dedup_samples = FALSE)
  res2 <- mve_exhaustive_search(mu, sigma, k = 2,
                                enumerate_all = FALSE, max_samples = 3,
                                dedup_samples = TRUE)
  expect_equal(length(res1$selection), 2)
  expect_equal(length(res2$selection), 2)
})

test_that("miqp search basic run (requires Gurobi)", {
  skip_if_no_gurobi()
  mu <- c(0.1, 0.2, 0.15)
  sigma <- diag(3)
  res <- mve_miqp_search(mu, sigma, k = 2, do_checks = TRUE, time_limit = 2)
  expect_true(length(res$selection) <= 2)
  expect_equal(length(res$weights), length(mu))
})

test_that("miqp exactly_k vs banded (requires Gurobi)", {
  skip_if_no_gurobi()
  mu <- c(0.1, 0.2, 0.15, 0.05)
  sigma <- diag(c(0.2, 0.25, 0.18, 0.22))
  res_exact <- mve_miqp_search(mu, sigma, k = 2, exactly_k = TRUE, time_limit = 2)
  res_banded <- mve_miqp_search(mu, sigma, k = 2, exactly_k = FALSE, m = 1, time_limit = 2)
  expect_equal(length(res_exact$selection), 2)
  expect_true(length(res_banded$selection) <= 2)
})

test_that("miqp budget/normalize_weights l1=1 (requires Gurobi)", {
  skip_if_no_gurobi()
  mu <- c(0.05, 0.08, 0.02)
  sigma <- diag(3) * 0.1
  res_budget <- mve_miqp_search(mu, sigma, k = 2, normalize_weights = TRUE, time_limit = 2)
  expect_lt(abs(sum(abs(res_budget$weights)) - 1), 1e-6)
})

test_that("lars search basic path", {
  skip_if_not_installed("lars")
  mu <- c(0.05, 0.08, 0.02, 0.10, 0.03)
  A <- matrix(rnorm(25), 5, 5)
  sigma <- crossprod(A) + 0.1 * diag(5)
  res <- mve_lars_search(mu, sigma, n_obs = 80, k = 2, do_checks = TRUE)
  expect_true(length(res$selection) <= 2)
  expect_equal(length(res$weights), length(mu))
  expect_true(is.finite(res$sr))
  expect_true(res$status %in% c("LARS_OK", "LARS_EXACT_K", "LARS_BELOW_K", "LARS_BELOW_K_FALLBACK", "LARS_DESIGN_FAIL", "LARS_FIT_FAIL"))
})

test_that("lars search respects k edges", {
  skip_if_not_installed("lars")
  mu <- c(0.1, 0.2, 0.15)
  sigma <- diag(c(0.2, 0.25, 0.18))

  res0 <- mve_lars_search(mu, sigma, n_obs = 40, k = 0, do_checks = TRUE)
  expect_equal(res0$selection, integer())
  expect_equal(res0$weights, numeric(length(mu)))
  expect_equal(res0$status, "LARS_EMPTY")

  resfull <- mve_lars_search(mu, sigma, n_obs = 40, k = length(mu), do_checks = TRUE)
  expect_equal(sort(resfull$selection), 1:length(mu))
  expect_equal(length(resfull$weights), length(mu))
  expect_true(resfull$status %in% c("LARS_FULL", "LARS_FULL_FAIL"))
})

test_that("lars search refit and normalization path", {
  skip_if_not_installed("lars")
  set.seed(101)
  mu <- c(0.05, 0.07, 0.02, 0.11)
  A <- matrix(rnorm(16), 4, 4)
  sigma <- crossprod(A) + 0.05 * diag(4)
  res <- mve_lars_search(mu, sigma, n_obs = 60, k = 2,
                         use_refit = TRUE, normalize_weights = TRUE, do_checks = TRUE)
  expect_equal(length(res$weights), length(mu))
  expect_lt(abs(sum(abs(res$weights)) - 1), 1e-6)
  expect_true(is.finite(res$sr))
})

test_that("lars search compute_weights=FALSE yields zero weights", {
  skip_if_not_installed("lars")
  set.seed(202)
  mu <- c(0.04, 0.09, 0.01, 0.08)
  A <- matrix(rnorm(16), 4, 4)
  sigma <- crossprod(A) + 0.05 * diag(4)
  res <- mve_lars_search(mu, sigma, n_obs = 50, k = 2,
                         compute_weights = FALSE, do_checks = TRUE)
  expect_true(all(res$weights == 0))
  expect_true(length(res$selection) <= 2)
  expect_true(is.na(res$sr))
})

test_that("lars search normalize_weights without refit scales weights", {
  skip_if_not_installed("lars")
  set.seed(303)
  mu <- c(0.03, 0.06, 0.02, 0.10, 0.05)
  A <- matrix(rnorm(25), 5, 5)
  sigma <- crossprod(A) + 0.05 * diag(5)
  res <- mve_lars_search(mu, sigma, n_obs = 70, k = 3,
                         normalize_weights = TRUE, use_refit = FALSE, do_checks = TRUE)
  expect_equal(length(res$weights), length(mu))
  expect_true(sum(abs(res$weights)) > 0)
  expect_lt(abs(sum(abs(res$weights)) - 1), 1e-6)
  expect_true(length(res$selection) <= 3)
})
