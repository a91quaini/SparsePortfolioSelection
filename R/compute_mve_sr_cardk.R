#' Compute Mean-Variance Efficient Sharpe Ratio with Cardinality K
#'
#' This function takes as inputs the mean vector \code{mu}, the covariance matrix \code{sigma},
#' the maximum active portfolio cardinality \code{max_card},
#' and the maximum number of combinations per cardinality to evaluate \code{max_comb}.
#' With these inputs, it searches over all combinations of assets with cardinality from 1 up to \code{max_card}
#' and computes the Sharpe ratio defined as
#' \eqn{\mu^T \, \sigma^{-1}\, \mu}.
#' It returns the highest Sharpe ratio found along with the associated asset selection.
#'
#' @param mu Mean vector.
#' @param sigma Coveriance matrix.
#' @param max_card Maximum investment cardinality (from 1 up to the number of assets).
#' @param max_comb Maximum number of combinations to consider. If 0 (default),
#'                 all combinations are computed.
#' @param do_checks Logical flag indicating whether to perform input checks (default = \code{FALSE}).
#' @return A list with \code{sr} (the optimal Sharpe ratio) and \code{selection} (the asset indices of the optimal selection).
#' @examples
#' # Consider a portfolio with 4 assets
#' mu <- c(0.1, 0.2, 0.15, 0.12)
#' sigma <- matrix(c(0.01, 0.002, 0.001, 0.0005,
#'                   0.002, 0.02, 0.0015, 0.001,
#'                   0.001, 0.0015, 0.015, 0.0007,
#'                   0.0005, 0.001, 0.0007, 0.012), nrow = 4, byrow = TRUE)
#' result <- compute_mve_sr_cardk(mu = mu,
#'                                sigma = sigma,
#'                                max_card = 2,
#'                                do_checks = TRUE)
#' @export
compute_mve_sr_cardk <- function(mu, sigma, max_card, max_comb = 0, do_checks = FALSE) {

  # Check inputs
  if (do_checks) {
    # Validate mu
    if (missing(mu) || length(mu) == 0) {
      stop("mu must be provided and non-empty")
    }
    if (!is.numeric(mu)) {
      stop("mu must be numeric")
    }

    # Validate sigma
    if (missing(sigma) || !is.matrix(sigma) || nrow(sigma) == 0) {
      stop("sigma must be provided as a non-empty matrix")
    }
    if (!is.numeric(sigma)) {
      stop("sigma must be numeric")
    }
    if (nrow(sigma) != ncol(sigma)) {
      stop("sigma must be a square matrix")
    }
    if (length(mu) != nrow(sigma)) {
      stop("The length of mu must equal the number of rows of sigma")
    }

    # Validate max_card
    if (missing(max_card) || length(max_card) == 0) {
      stop("max_card must be provided")
    }
    if (max_card < 1 || max_card > length(mu)) {
      stop("max_card must be between 1 and the number of assets (length of mu)")
    }

    # Validate max_comb"
    if (length(max_comb) == 0) {
      stop("max_comb must be provided")
    }
    if (!is.numeric(max_comb)) {
      stop("max_comb must be numeric")
    }
    if (max_comb < 0) {
      stop("max_comb must be non-negative")
    }
  }

  .Call(`_SparsePortfolioSelection_compute_mve_sr_cardk_cpp`, mu, sigma, as.integer(max_card), as.integer(max_comb), FALSE)
}
