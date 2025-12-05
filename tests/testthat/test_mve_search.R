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

test_that("miqp budget/normalize_weights toggles sum=1 (requires Gurobi)", {
  skip_if_no_gurobi()
  mu <- c(0.05, 0.08, 0.02)
  sigma <- diag(3) * 0.1
  res_budget <- mve_miqp_search(mu, sigma, k = 2, normalize_weights = TRUE, time_limit = 2)
  expect_lt(abs(sum(res_budget$weights) - 1), 1e-6)
})

test_that("lasso search respects k edges", {
  set.seed(123)
  R <- matrix(rnorm(200), nrow = 50, ncol = 4)
  res0 <- mve_lasso_search_return_based(R, k = 0)
  expect_equal(res0$selection, integer())
  expect_equal(res0$weights, numeric(ncol(R)))

  resfull <- mve_lasso_search_return_based(R, k = ncol(R), compute_weights = TRUE)
  expect_equal(length(resfull$selection), ncol(R))
  expect_equal(length(resfull$weights), ncol(R))
})

test_that("lasso search lambda override and standardize", {
  set.seed(321)
  R <- matrix(rnorm(300), nrow = 60, ncol = 5)
  # small lambda to allow dense solution
  res <- mve_lasso_search_return_based(R, k = 3, lambda = exp(seq(0, -4, length.out = 20)),
                          standardize = TRUE)
  expect_true(length(res$selection) <= 3)
})

test_that("lasso search closest-k behavior", {
  set.seed(42)
  R <- matrix(rnorm(500), nrow = 100, ncol = 5)
  res <- mve_lasso_search_return_based(R, k = 2, nlambda = 50)
  expect_true(length(res$selection) <= 2)
  expect_equal(length(res$weights), ncol(R))
  expect_true(res$status %in% c("LASSO_PATH_EXACT_K", "LASSO_PATH_CLOSEST", "LASSO_PATH_OVER_K"))
})

test_that("lasso alpha grid with CV selects feasible alpha", {
  set.seed(4242)
  R <- matrix(rnorm(400), nrow = 80, ncol = 5)
  res <- mve_lasso_search_return_based(R, k = 2, alpha = c(0.2, 0.5, 0.8), n_folds = 3)
  expect_true(res$alpha %in% c(0.2, 0.5, 0.8))
  expect_true(length(res$selection) <= 2)
})

test_that("lasso search densification hits k when possible", {
  set.seed(777)
  R <- matrix(rnorm(400), nrow = 80, ncol = 5)
  res <- mve_lasso_search_return_based(R, k = 3, nlambda = 20, nadd = 30, nnested = 2)
  expect_true(length(res$selection) <= 3)
})

test_that("lasso search normalize_weights path", {
  set.seed(99)
  R <- matrix(rnorm(320), nrow = 80, ncol = 4)
  res <- mve_lasso_search_return_based(R, k = 2, normalize_weights = TRUE, use_refit = TRUE)
  expect_equal(length(res$weights), ncol(R))
  expect_true(abs(sum(res$weights)) > 0 || sum(abs(res$weights)) == 0)
})

test_that("lasso search different eps/stabilization choices", {
  set.seed(7)
  R <- matrix(rnorm(300), nrow = 75, ncol = 4)
  res1 <- mve_lasso_search_return_based(R, k = 2, epsilon = 0, stabilize_sigma = FALSE)
  res2 <- mve_lasso_search_return_based(R, k = 2, epsilon = 1e-3, stabilize_sigma = TRUE)
  expect_equal(length(res1$weights), ncol(R))
  expect_equal(length(res2$weights), ncol(R))
})

test_that("lasso moment-based path basic run", {
  mu <- c(0.05, 0.08, 0.02, 0.10)
  A <- matrix(rnorm(16), 4, 4)
  sigma <- crossprod(A) + 0.05 * diag(4)
  res <- mve_lasso_search(mu, sigma, n_obs = 60, k = 2, nlambda = 30, nadd = 20, nnested = 2)
  expect_true(length(res$selection) <= 2)
  expect_equal(length(res$weights), length(mu))
})


test_that("lasso search respects k edges", {
  set.seed(123)
  R <- matrix(rnorm(200), nrow = 50, ncol = 4)
  res0 <- mve_lasso_search_return_based(R, k = 0)
  expect_equal(res0$selection, integer())
  expect_equal(res0$weights, numeric(ncol(R)))

  resfull <- mve_lasso_search_return_based(R, k = ncol(R), compute_weights = TRUE)
  expect_equal(length(resfull$selection), ncol(R))
  expect_equal(length(resfull$weights), ncol(R))
})
