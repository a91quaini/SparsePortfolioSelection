library(testthat)
library(SparsePortfolioSelection)

TOL <- 1e-7

test_that("basic sanity", {
  mu <- c(0.10, 0.20, 0.15)
  A  <- matrix(c(0.2, 0.05, 0.0,
                 0.05,0.3, 0.02,
                 0.0, 0.02,0.25), nrow = 3, byrow = TRUE)
  sigma <- A
  w <- rep(1/3, 3)

  expect_true(is.finite(compute_sr(w, mu, sigma)))
  expect_gt(compute_mve_sr(mu, sigma), 0)
  expect_length(compute_mve_weights(mu, sigma), 3)
})

test_that("selection semantics", {
  mu <- c(0.1, 0.2, 0.3, 0.4)
  sigma <- matrix(c(0.08,0.01,0.00,0.00,
                    0.01,0.09,0.01,0.00,
                    0.00,0.01,0.07,0.02,
                    0.00,0.00,0.02,0.10), nrow = 4, byrow = TRUE)
  sel <- c(2,3,4) - 1L  # 0-based for C++ selection
  n <- length(mu)
  w <- runif(n); w <- w / (sum(abs(w)) + .Machine$double.eps)

  sr_sel <- compute_sr(w, mu, sigma, selection = sel, do_checks = TRUE)
  sigma_sub <- sigma[sel + 1L, sel + 1L, drop = FALSE]
  sr_dir <- compute_sr(w[sel + 1L], mu[sel + 1L], sigma_sub, do_checks = TRUE)
  expect_lt(abs(sr_sel - sr_dir), TOL)

  w_sel <- compute_mve_weights(mu, sigma, selection = sel)
  expect_equal(sum(w_sel != 0), length(sel))
  expect_true(all(w_sel[setdiff(seq_len(n), sel + 1L)] == 0))

  sr_all <- compute_mve_sr(mu, sigma)
  sr_fullS <- compute_mve_sr(mu, sigma, selection = 0:(n-1))
  expect_lt(abs(sr_all - sr_fullS), TOL)
})

test_that("scale behavior", {
  set.seed(1)
  n <- 5
  mu <- 0.05 + 0.10 * runif(n)
  A <- matrix(rnorm(n*n), n)
  sigma <- A %*% t(A) + 0.1 * diag(n)

  base <- compute_mve_sr(mu, sigma)
  alpha <- 2; beta <- 3
  scaled <- compute_mve_sr(alpha * mu, beta * sigma)
  expect_lt(abs(scaled - (alpha / sqrt(beta)) * base), 100*TOL)

  w <- runif(n); w <- w / (sum(abs(w)) + .Machine$double.eps)
  sr_base <- compute_sr(w, mu, sigma)
  sr_scale <- compute_sr(w, alpha * mu, beta * sigma)
  expect_lt(abs(sr_scale - (alpha / sqrt(beta)) * sr_base), 100*TOL)

  # scale invariance of compute_sr when weights scaled
  cst <- 3.7
  expect_lt(abs(compute_sr(cst * w, mu, sigma) - compute_sr(w, mu, sigma)), 1e-12)
})

test_that("degenerate covariance & ridge", {
  n <- 6
  v <- rep(1, n)
  sigma <- v %*% t(v)  # rank 1
  mu <- seq(0.1, 0.6, by = 0.1)

  sr_mve <- compute_mve_sr(mu, sigma, epsilon = 0)
  expect_true(is.finite(sr_mve))
  expect_gte(sr_mve, 0)

  w <- compute_mve_weights(mu, sigma, epsilon = 0)
  sr_w <- compute_sr(w, mu, sigma, epsilon = 0)
  expect_lt(abs(sr_w - sr_mve), 1e-5)

  sigma_ns <- sigma + 1e-4 * diag(n)
  sr_eps <- compute_mve_sr(mu, sigma_ns, epsilon = 1e-2)
  expect_true(is.finite(sr_eps))
})

test_that("zero covariance corner case", {
  mu <- c(1, 2)
  sigma <- matrix(0, 2, 2)
  w <- c(0.5, 0.5)

  expect_true(is.nan(compute_sr(w, mu, sigma, epsilon = 0)))
  expect_true(is.nan(compute_sr(w, mu, sigma, epsilon = 1e-2)))

  expect_equal(compute_mve_sr(mu, sigma, epsilon = 0), 0)
  expect_equal(compute_mve_sr(mu, sigma, epsilon = 1e-2), 0)
  expect_equal(as.numeric(compute_mve_weights(mu, sigma, epsilon = 0, stabilize_sigma = FALSE)), c(0, 0))
})

test_that("single-asset edge case", {
  mu <- 0.2
  sigma <- matrix(0.04, 1, 1)
  w <- 1
  expect_lt(abs(compute_sr(w, mu, sigma, epsilon = 0) - (mu / sqrt(0.04))), 1e-12)
  expect_lt(abs(compute_mve_sr(mu, sigma, epsilon = 0) - (abs(mu) / sqrt(0.04))), 1e-12)
  expect_equal(abs(compute_mve_weights(mu, sigma, epsilon = 0, stabilize_sigma = FALSE)), mu / 0.04)
})

test_that("weights consistency", {
  set.seed(2)
  n <- 7
  mu <- 0.02 + 0.05 * runif(n)
  A <- matrix(rnorm(n*n), n)
  sigma <- A %*% t(A) + 0.2 * diag(n)
  wstar <- compute_mve_weights(mu, sigma)
  sr_wstar <- compute_sr(wstar, mu, sigma)
  mve <- compute_mve_sr(mu, sigma)
  expect_lt(abs(sr_wstar - mve), 1e-7)

  # selection consistency with explicit solve
  sel <- sort(sample(seq_len(n), 4)) - 1L
  w_sel <- compute_mve_weights(mu, sigma, selection = sel, stabilize_sigma = FALSE, epsilon = 0)
  sig_s <- sigma[sel + 1L, sel + 1L, drop = FALSE]
  w_direct <- solve(sig_s, mu[sel + 1L])
  expect_true(all.equal(w_sel[sel + 1L], w_direct, tolerance = 1e-10) == TRUE)
})

test_that("stabilize_sigma semantics", {
  set.seed(3)
  n <- 5
  A <- matrix(rnorm(n*n), n)
  sigma_raw <- A %*% t(A) + 1e-1 * diag(n)
  sigma_pert <- sigma_raw + 1e-10 * matrix(rnorm(n*n), n)
  mu <- runif(n)

  s1 <- compute_mve_sr(mu, sigma_pert, stabilize_sigma = TRUE, epsilon = 0)
  s2 <- compute_mve_sr(mu, sigma_pert, stabilize_sigma = FALSE, epsilon = 0)
  expect_true(is.finite(s1) && is.finite(s2))
  expect_lt(abs(s1 - s2), 1e-8)
})

test_that("normalize_weights semantics", {
  set.seed(4)
  n <- 8
  mu <- 0.01 + 0.03 * runif(n)
  A <- matrix(rnorm(n*n), n)
  sigma <- A %*% t(A) + 0.15 * diag(n)

  w0 <- compute_mve_weights(mu, sigma, normalize_weights = FALSE)
  w1 <- compute_mve_weights(mu, sigma, normalize_weights = TRUE)

  denom <- max(abs(sum(w0)), 1e-6 * sum(abs(w0)), 1e-10)
  expect_true(all.equal(w1, w0 / denom, tolerance = 1e-10) == TRUE)
  expect_lt(abs(compute_sr(w0, mu, sigma) - compute_sr(w1, mu, sigma)), 1e-10)

  sel <- sort(sample(seq_len(n), 5)) - 1L
  w_raw <- compute_mve_weights(mu, sigma, selection = sel, normalize_weights = FALSE)
  ws <- compute_mve_weights(mu, sigma, selection = sel, normalize_weights = TRUE)
  expect_true(all(ws[setdiff(seq_len(n), sel + 1L)] == 0))
  denom_sel <- max(abs(sum(w_raw)), 1e-6 * sum(abs(w_raw)), 1e-10)
  expect_true(all.equal(ws, w_raw / denom_sel, tolerance = 1e-10) == TRUE)
  expect_lt(abs(compute_sr(ws, mu, sigma) - compute_sr(w_raw, mu, sigma)), 1e-10)
})

test_that("stabilize_sigma toggling leaves symmetric matrix unchanged", {
  mu <- c(0.1, 0.2, 0.15)
  sigma <- diag(c(0.2, 0.3, 0.25))
  sr1 <- compute_mve_sr(mu, sigma, stabilize_sigma = TRUE, epsilon = 0)
  sr2 <- compute_mve_sr(mu, sigma, stabilize_sigma = FALSE, epsilon = 0)
  expect_lt(abs(sr1 - sr2), 1e-12)
})

test_that("argument checks", {
  mu <- c(0.1, 0.2)
  sigma <- matrix(c(0.04, 0.01, 0.01, 0.09), 2, 2)
  w <- c(0.5, 0.5, 0.5)
  expect_error(compute_sr(w, mu, sigma, do_checks = TRUE))
  bad_sel <- c(0, 3)
  expect_error(compute_mve_sr(mu, sigma, selection = bad_sel, do_checks = TRUE))
  expect_error(compute_mve_weights(mu, sigma, selection = bad_sel, do_checks = TRUE))
  expect_error(compute_mve_sr(mu, matrix(1, 2, 3), do_checks = TRUE))
})
