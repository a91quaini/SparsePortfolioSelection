#' Compute portfolio Sharpe ratio (C++ backend)
#'
#' Computes \eqn{SR = (w^\top \mu) / \sqrt{w^\top \Sigma w}} with optional
#' asset selection and covariance stabilization. Delegates to the C++ routine
#' `compute_sr_cpp`, which symmetrizes and ridge-stabilizes the covariance when
#' requested.
#'
#' @param w Numeric weight vector of length \eqn{N}.
#' @param mu Numeric mean vector of length \eqn{N}.
#' @param sigma Numeric \eqn{N \times N} covariance matrix.
#' @param selection Optional integer vector of asset indices (0-based) on which
#'   to evaluate the Sharpe ratio. Default uses all assets.
#' @param ridge_epsilon Ridge stabilization factor.
#' @param do_checks Logical; if `TRUE`, perform basic input validation.
#'
#' @return Scalar Sharpe ratio (NaN if variance is nonpositive or nonfinite).
#' @export
compute_sr <- function(w, mu, sigma,
                       selection = integer(),
                       ridge_epsilon = 0.0,
                       do_checks = FALSE) {
  compute_sr_cpp(w, mu, sigma, selection, ridge_epsilon, do_checks)
}

#' Compute mean–variance efficient (MVE) Sharpe ratio
#'
#' Evaluates \eqn{\sqrt{\mu^\top \Sigma^{-1} \mu}} on the full universe or a
#' subset. Uses the C++ routine `compute_mve_sr_cpp` with optional covariance
#' stabilization for numerical safety.
#'
#' @inheritParams compute_sr
#'
#' @return Scalar MVE Sharpe ratio (nonnegative).
#' @export
compute_mve_sr <- function(mu, sigma,
                           selection = integer(),
                           ridge_epsilon = 0.0,
                           do_checks = FALSE) {
  compute_mve_sr_cpp(mu, sigma, selection, ridge_epsilon, do_checks)
}

#' Compute MVE weights \eqn{w = \Sigma^{-1} \mu}
#'
#' Solves for the mean–variance efficient weights, optionally normalizing with
#' the stable `normalize_weights_cpp()` scaling. Supports
#' computing on a subset and zero-filling off-support entries.
#'
#' @inheritParams compute_sr
#' @param normalize_weights Logical; if `TRUE`, post-scale weights using
#'   L1-relative normalization.
#'
#' @return Numeric weight vector of length \eqn{N}.
#' @export
compute_mve_weights <- function(mu, sigma,
                                selection = integer(),
                                ridge_epsilon = 0.0,
                                normalize_weights = FALSE,
                                normalization_type = 1L,
                                do_checks = FALSE) {
  as.numeric(compute_mve_weights_cpp(mu, sigma, selection, ridge_epsilon,
                                     normalize_weights, normalization_type,
                                     do_checks))
}
