#' Compute MVE Weights
#'
#' Given the risk-aversion parameter \eqn{\gamma}, the first moment vector, the covariance matrix,
#' and the selection index vector, this function computes the mean variance efficient portfolio weights
#' \deqn{w = \frac{1}{\gamma}\, \text{sigma}^{-1}\, \mu},
#' over the selected assets.
#' If the provided asset selection vector has length less than N,
#' the function returns an N-length weight vector with zeros for the unselected assets.
#'
#' @param mu First moment vector
#' @param sigma Covariance matrix
#' @param selection Asset selection vector (default = full selection)
#' @param gamma Risk aversion parameter (default = 1)
#' @param do_checks Logical flag to perform input checks (default = FALSE)
#' @return A vector of portfolio weights.
#' @examples
#' # Full selection example:
#' compute_mve_weights(mu = c(0.1, 0.2, 0.15),
#'             sigma = matrix(c(1, 0.2, 0.1,
#'                                      0.2, 1, 0.3,
#'                                      0.1, 0.3, 1), nrow = 3),
#'             selection = c(1, 2, 3),
#'             gamma = 1,
#'             do_checks = TRUE)
#'
#' # Subset selection example:
#' compute_mve_weights(mu = c(0.1, 0.2, 0.15, 0.12),
#'             sigma = diag(4),
#'             selection = c(1, 3),
#'             gamma = 1,
#'             do_checks = TRUE)
#' @export
compute_mve_weights <- function(mu,
                                sigma,
                                selection = c(),
                                gamma = 1,
                                do_checks = FALSE) {

  # Check if selection is empty
  if (missing(selection) || length(selection) == 0) {
    # Default to full selection
    selection <- 1:length(mu)
  }

  # Check inputs
  if (do_checks) {
    # Validate mu
    if (missing(mu) || length(mu) == 0) {
      stop("mu must be provided and non-empty")
    }
    if (!is.numeric(mu)) {
      stop("mu must be numeric")
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
    if (length(mu) != nrow(sigma)) {
      stop("The length of mu must equal the number of rows of sigma")
    }

    # Validate selection.
    if (!is.numeric(selection)) {
      stop("selection must be numeric or integer")
    }
    selection <- as.integer(selection)
    if (min(selection) < 1 || max(selection) > length(mu)) {
      stop("Asset selection indices out of bounds")
    }
  }

  .Call(`_SparsePortfolioSelection_compute_mve_weights_cpp`, mu, sigma, as.integer(selection - 1), gamma, FALSE)
}
