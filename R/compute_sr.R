#' Compute Sharpe Ratio
#'
#' Given a vector of portfolio weights, a vector of expected returns (mu),
#' and a covariance matrix (sigma), this function computes the Sharpe ratio defined as
#' \deqn{\frac{w^T \mu}{\sqrt{w^T \Sigma w}}.}
#'
#' When \code{do_checks} is TRUE, the wrapper performs input validation in R before calling the C++ function.
#'
#' @param weights Numeric vector of portfolio weights. Must be non-empty.
#' @param mu Numeric vector of expected returns. Must be non-empty and of the same length as \code{weights}.
#' @param sigma Numeric covariance matrix. Must be non-empty, square, and its dimensions must match the length of \code{weights}.
#' @param do_checks Logical flag indicating whether to perform input checks (default = FALSE).
#' @return A numeric scalar representing the Sharpe ratio.
#' @examples
#' # Example: Two-asset portfolio
#' weights <- c(0.5, 0.5)
#' mu <- c(0.1, 0.2)
#' sigma <- diag(2)
#' compute_sr(weights, mu, sigma, do_checks = TRUE)
#' @export
compute_sr <- function(weights, mu, sigma, do_checks = FALSE) {
  if (do_checks) {
    # Validate weights.
    if (missing(weights) || length(weights) == 0) {
      stop("weights must be provided and non-empty")
    }
    if (!is.numeric(weights)) {
      stop("weights must be numeric")
    }

    # Validate mu.
    if (missing(mu) || length(mu) == 0) {
      stop("mu must be provided and non-empty")
    }
    if (!is.numeric(mu)) {
      stop("mu must be numeric")
    }
    if (length(weights) != length(mu)) {
      stop("weights and mu must be of the same length")
    }

    # Validate sigma.
    if (missing(sigma) || !is.matrix(sigma) || nrow(sigma) == 0) {
      stop("sigma must be provided as a non-empty matrix")
    }
    if (!is.numeric(sigma)) {
      stop("sigma must be numeric")
    }
    if (nrow(sigma) != ncol(sigma)) {
      stop("sigma must be a square matrix")
    }
    if (nrow(sigma) != length(weights)) {
      stop("The dimensions of sigma must match the length of weights")
    }
  }

  .Call(`_SparsePortfolioSelection_compute_sr_cpp`, weights, mu, sigma, FALSE)
}
