#' Load Returns or Factors Data
#'
#' This function loads and assembles data based on the specified type. The input
#' \code{type} can be one of the following:
#'
#' "p" loads all double-sorted and industry portfolio objects from the SparsePortfolioSelection package.
#' The function removes the Date column from each portfolio object and aggregates their returns (columns represent portfolios)
#' into a single numeric matrix.
#'
#' "i" loads the CRSP individual asset returns data from the SparsePortfolioSelection package.
#' It removes the Date column and returns a numeric matrix of individual asset returns.
#'
#' "f" loads the Fama–French 5 factors (tradable version) data from the SparsePortfolioSelection package.
#' It removes the Date column and returns a numeric matrix of factors.
#'
#' All data are assumed to have been filtered to include only observations with dates from May 2008 to December 2022.
#'
#' @param type A character string indicating the type of data to load. Valid options are "p" for portfolios, "i" for individual assets, and "f" for factors.
#' @param do_checks Logical. If \code{TRUE} (the default), the function checks that \code{type} is one of "p", "i", or "f".
#'
#' @return A numeric matrix named \code{data} containing the data with the Date column removed. For \code{type = "p"}, each column corresponds to a portfolio's excess returns; for \code{type = "i"}, the matrix contains individual asset excess returns; for \code{type = "f"}, the matrix contains the Fama–French 5 factors (tradable).
#'
#' @examples
#' \dontrun{
#'   # Load portfolio returns
#'   data <- load_data(type = "p")
#'
#'   # Load individual asset returns (CRSP)
#'   data <- load_data(type = "i")
#'
#'   # Load factor data (Fama–French 5)
#'   data <- load_data(type = "f")
#' }
#'
#' @export
load_data <- function(type = "p", do_checks = TRUE) {
  # Check that type is valid
  if (do_checks) {
    if (!is.character(type) || !(type %in% c("p", "i", "f"))) {
      stop("Invalid 'type'. It must be 'p' (for portfolios), 'i' (for individual assets), or 'f' (for factors).")
    }
  }

  if (type == "p") {
    # Load portfolio objects from the SparsePortfolioSelection package and combine them.
    data <- cbind(
      SparsePortfolioSelection::returns_mebeme25[,-1, drop = FALSE],
      SparsePortfolioSelection::returns_meop25[,-1, drop = FALSE],
      SparsePortfolioSelection::returns_opinv25[,-1, drop = FALSE],
      SparsePortfolioSelection::returns_bemeinv25[,-1, drop = FALSE],
      SparsePortfolioSelection::returns_bemeop25[,-1, drop = FALSE],
      SparsePortfolioSelection::returns_meac25[,-1, drop = FALSE],
      SparsePortfolioSelection::returns_mebeta25[,-1, drop = FALSE],
      SparsePortfolioSelection::returns_meinv25[,-1, drop = FALSE],
      SparsePortfolioSelection::returns_meprior10[,-1, drop = FALSE],
      SparsePortfolioSelection::returns_meprior122[,-1, drop = FALSE],
      SparsePortfolioSelection::returns_meprior6013[,-1, drop = FALSE],
      SparsePortfolioSelection::returns_mevar25[,-1, drop = FALSE],
      SparsePortfolioSelection::returns_ind17[,-1, drop = FALSE]
    )
    return(data)

  } else if (type == "i") {
    data <- SparsePortfolioSelection::returns_crsp[,-1, drop = FALSE]
    return(data)

  } else if (type == "f") {
    data <- SparsePortfolioSelection::factors_ff5[,-1, drop = FALSE]
    return(data)
  }
}

#' Simulate Sharpe Ratio Loss
#'
#' This function simulates a sample of size \code{n_obs} from a multivariate normal
#' distribution with mean vector \code{mu} and covariance matrix \code{sigma}.
#' It computes the sample mean vector (\code{mu_sample}) and the sample covariance matrix (\code{sigma_sample}),
#' then with the population and sample parameters,
#' the maximum cardinality (\code{max_card}), and the maximum number of combinations (\code{max_comb}),
#' it computes: the sample MVE Sharpe ratio, the sample MVE Sharpe ratio with
#' maximum cardinality k, and the population MVE Sharpe ratio with maximum
#' cardinality k decomposed in estimation and selection component.
#'
#' @param mu A numeric vector; the population mean vector.
#' @param sigma A numeric matrix; the population covariance matrix.
#' @param n_obs An integer specifying the sample size to simulate.
#' @param max_card Maximum cardinality to consider (from 1 up to the number of assets).
#' @param max_comb Maximum number of combinations to consider. If 0 (default),
#'                 all combinations are computed.
#' @param do_checks Logical; if TRUE, input checks are performed.
#'
#' @return A list with:
#'    - \code{sample_mve_sr} computed as the optimal sample mve Sharpe ratio,
#'    - \code{sample_mve_sr_cardk} computed as the optimal sample mve Sharpe ratio
#'      with cardinality \code{max_card},
#'    - \code{mve_sr_cardk_est_term} computed as \eqn{w^T \mu^T / \sqrt{w^T\sigma w}}
#'      where \code{w} are the optimal sample mve weights,
#'    - \code{mve_sr_cardk_sel_term} computed as \eqn{\mu_S^T  \sigma_S^{-1} \mu_S}
#'      where \code{S} is the set of assets yielding the optimal sample mve Sharpe ratio.
#' @examples
#' \dontrun{
#' # Example with three assets:
#' mu <- c(0.1, 0.2, 0.15)
#'   sigma <- diag(3)
#'   mve_sr <- 0.8
#'   n_obs <- 100
#'   max_card <- 2
#'   result <- simulate_mve_sr(mu,
#'                             sigma,
#'                             n_obs,
#'                             max_card,
#'                             do_checks = TRUE)
#'   print(result)
#' }
#' @export
simulate_mve_sr <- function(mu, sigma, n_obs, max_card, max_comb = 0, do_checks = FALSE) {
  # Perform input validation if do_checks is TRUE
  if (do_checks) {
    if (missing(mu) || length(mu) == 0 || !is.numeric(mu)) {
      stop("mu must be provided and be a non-empty numeric vector")
    }
    if (missing(sigma) || !is.matrix(sigma) || nrow(sigma) == 0 || !is.numeric(sigma)) {
      stop("sigma must be provided as a non-empty numeric matrix")
    }
    if (nrow(sigma) != ncol(sigma)) {
      stop("sigma must be a square matrix")
    }
    if (length(mu) != nrow(sigma)) {
      stop("The length of mu must equal the number of rows of sigma")
    }
    if (missing(n_obs) || !is.numeric(n_obs) || n_obs < 1) {
      stop("n_obs must be a positive integer")
    }
    if (missing(max_card) || !is.numeric(max_card) || max_card < 1) {
      stop("max_card must be a positive integer")
    }
    max_card <- as.integer(max_card)
    if (max_card > length(mu)) {
      stop("max_card cannot exceed the number of assets (length of mu)")
    }
    if (missing(max_comb) || !is.numeric(max_comb) || max_comb < 0) {
      stop("max_comb must be a nonnegative numeric scalar")
    }
  }

  .Call(`_SparsePortfolioSelection_simulate_mve_sr_cpp`, mu, sigma, n_obs, max_card, max_comb, FALSE)

}
