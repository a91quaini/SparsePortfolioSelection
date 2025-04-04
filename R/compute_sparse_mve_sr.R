#' Compute Optimal Sharpe Ratio with Cardinality Constraint
#'
#' This function takes as inputs the mean vector \code{mu}, the covariance matrix \code{sigma},
#' the maximum active portfolio cardinality \code{max_card},
#' and the fraction of combinations to evaluate \code{greedy_perc}.
#' With these inputs, it searches over all combinations of assets with cardinality from 1 up to \code{max_card}
#' and computes the square Sharpe ratio defined as
#' \eqn{\mu^T \, \sigma^{-1}\, \mu}.
#' It returns the highest square Sharpe ratio found along with the associated asset selection
#' If \code{greedy_perc} is less than 1, then for each cardinality the search is performed over a random
#' subset (a fraction equal to \code{greedy_perc}) of all possible combinations.
#'
#' @param mu Mean vector
#' @param sigma Coveriance matrix
#' @param max_card Maximum cardinality to consider (from 1 up to the number of assets)
#' @param greedy_perc If less than 1, the fraction of combinations to evaluate for each cardinality;
#'                    if 1 or greater, all combinations are evaluated
#' @param do_checks Logical flag indicating whether to perform input checks (default = FALSE)
#' @return A list with \code{sqsr} (the optimal square Sharpe ratio) and \code{selection} (the asset indices of the optimal selection)
#' @examples
#' # Consider a portfolio with 4 assets
#' mu <- c(0.1, 0.2, 0.15, 0.12)
#' sigma <- matrix(c(0.01, 0.002, 0.001, 0.0005,
#'                   0.002, 0.02, 0.0015, 0.001,
#'                   0.001, 0.0015, 0.015, 0.0007,
#'                   0.0005, 0.001, 0.0007, 0.012), nrow = 4, byrow = TRUE)
#' result <- compute_sparse_mve_sr(mu = mu,
#'                                 sigma = sigma,
#'                                 max_card = 2,
#'                                 greedy_perc = 1.0,
#'                                 do_checks = TRUE)
#' @export
compute_sparse_mve_sr <- function(mu, sigma, max_card, greedy_perc, do_checks = FALSE) {

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
    if (!is.numeric(max_card) || max_card < 1) {
      stop("max_card must be a positive number")
    }
    max_card <- as.integer(max_card)
    if (max_card > length(mu)) {
      stop("max_card cannot exceed the number of assets (length of mu)")
    }

    # Validate greedy_perc"
    if (missing(greedy_perc) || length(greedy_perc) == 0) {
      stop("greedy_perc must be provided")
    }
    if (!is.numeric(greedy_perc)) {
      stop("greedy_perc must be numeric")
    }
    if (greedy_perc < 0) {
      stop("greedy_perc must be non-negative")
    }
  }

  .Call(`_SparsePortfolioSelection_compute_sparse_mve_sr`, mu, sigma, as.integer(max_card), greedy_perc, FALSE)
}
