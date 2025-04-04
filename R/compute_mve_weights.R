#' Compute MVE Weights
#'
#' Given a risk-aversion parameter \eqn{\gamma}, a first moment vector, a second moment matrix,
#' and a selection index vector, this function computes the mean variance efficient portfolio weights
#' \deqn{w = \frac{1}{\gamma}\, \text{second\_moment}^{-1}\, \text{first\_moment}},
#' over the selected assets.
#' If the provided asset selection vector has length less than N,
#' the function returns an N-length weight vector with zeros for the unselected assets.
#'
#' @param first_moment First moment vector
#' @param second_moment Second moment matrix
#' @param selection Asset selection vector (0-indexed)
#' @param gamma Risk aversion parameter (default = 1)
#' @param do_checks Logical flag to perform input checks (default = FALSE)
#' @return A vector of portfolio weights.
#' @examples
#' # Full selection example:
#' compute_mve_weights(first_moment = c(0.1, 0.2, 0.15),
#'             second_moment = matrix(c(1, 0.2, 0.1,
#'                                      0.2, 1, 0.3,
#'                                      0.1, 0.3, 1), nrow = 3),
#'             selection = c(0, 1, 2),
#'             gamma = 1,
#'             do_checks = TRUE)
#'
#' # Subset selection example:
#' compute_mve_weights(first_moment = c(0.1, 0.2, 0.15, 0.12),
#'             second_moment = diag(4),
#'             selection = c(0, 2),
#'             gamma = 1,
#'             do_checks = TRUE)
#' @export
compute_mve_weights <- function(first_moment,
                        second_moment,
                        selection,
                        gamma = 1,
                        do_checks = FALSE) {

  # Check inputs
  if (do_checks) {
    # Validate first_moment.
    if (missing(first_moment) || length(first_moment) == 0) {
      stop("first_moment must be provided and non-empty")
    }
    if (!is.numeric(first_moment)) {
      stop("first_moment must be numeric")
    }

    # Validate second_moment.
    if (missing(second_moment) || !is.matrix(second_moment) || nrow(second_moment) == 0) {
      stop("second_moment must be provided as a non-empty matrix")
    }
    if (!is.numeric(second_moment)) {
      stop("second_moment must be numeric")
    }
    if (nrow(second_moment) != ncol(second_moment)) {
      stop("second_moment must be a square matrix")
    }
    if (length(first_moment) != nrow(second_moment)) {
      stop("The length of first_moment must equal the number of rows of second_moment")
    }

    # Validate selection.
    if (missing(selection) || length(selection) == 0) {
      # Default to full selection (0-indexed)
      selection <- 0:(length(first_moment) - 1)
    } else {
      if (!is.numeric(selection)) {
        stop("selection must be numeric or integer")
      }
      selection <- as.integer(selection)
      if (min(selection) < 0 || max(selection) > (length(first_moment) - 1)) {
        stop("Asset selection indices out of bounds")
      }
    }
  }

  .Call(`_SparsePortfolioSelection_compute_mve_weights`, first_moment, second_moment, as.integer(selection), gamma, FALSE)
}
