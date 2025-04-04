#' Compute Optimal Portfolio Sharpe Ratio
#'
#' Given a mean vector \eqn{\mu}, a covariance matrix \eqn{\Sigma},
#' and an asset selection vector, this function computes the optimal square
#' Sharpe ratio defined as
#' \deqn{\mu^T \Sigma^{-1}\mu}
#' over the selected assets.
#' If the provided asset selection vector has length less than N,
#' the function uses the subset of assets and computes the ratio using those assets.
#'
#' @param mu Mean vector.
#' @param sigma Covariance matrix.
#' @param selection Asset selection vector (default = full set, 0-indexed).
#' @param do_checks Logical flag to perform input checks (default = FALSE).
#' @return A scalar value corresponding to \eqn{\sqrt{\mu^T \Sigma^{-1}\mu}}.
#' @examples
#' # Full selection example:
#' mve_sqsr(mu = c(0.1, 0.2, 0.15),
#'                        sigma = matrix(c(1, 0.2, 0.1,
#'                                           0.2, 1, 0.3,
#'                                           0.1, 0.3, 1), nrow = 3),
#'                        selection = c(0, 1, 2),
#'                        do_checks = TRUE)
#'
#' # Subset selection example:
#' mve_sqsr(mu = c(0.1, 0.2, 0.15, 0.12),
#'                        sigma = diag(4),
#'                        selection = c(0, 2),
#'                        do_checks = TRUE)
#' @export
mve_sqsr <- function(mu, sigma, selection, do_checks = FALSE) {

  # Check inputs
  if (do_checks) {
    # Check that mu is provided and numeric.
    if (missing(mu) || length(mu) == 0) {
      stop("mu must be provided and non-empty")
    }
    if (!is.numeric(mu)) {
      stop("mu must be numeric")
    }

    # Check that sigma is provided, is a numeric matrix, and is non-empty.
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

    # Check selection.
    if (missing(selection) || length(selection) == 0) {
      # If no selection is provided, default to full selection.
      selection <- 0:(length(mu) - 1)
    } else {
      if (!is.numeric(selection)) {
        stop("selection must be numeric or integer")
      }
      selection <- as.integer(selection)
      if (min(selection) < 0 || max(selection) > (length(mu) - 1)) {
        stop("Asset selection indices out of bounds")
      }
    }
  }

  .Call(`_SparsePortfolioSelection_mve_sqsr`, mu, sigma, as.integer(selection), FALSE)
}
