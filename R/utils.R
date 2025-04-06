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
#' This function serves as a wrapper that simulates a sample of returns from a multivariate normal
#' distribution with population mean vector \code{mu} and covariance matrix \code{sigma} (sample size \code{n_obs}),
#' then computes the sample mean vector and sample covariance matrix internally using a C++ implementation.
#' It then calls the Rcpp function \code{simulate_sr_loss_rcpp} (which uses the population parameters, the sample
#' estimates, a maximum cardinality (\code{max_card}), and a greedy percentage (\code{greedy_perc}))
#' to compute the Sharpe ratio loss.
#'
#' When \code{do_checks} is TRUE, input validation is performed to ensure that:
#' \code{mu} is a non-empty numeric vector; \code{sigma} is a non-empty numeric square matrix with dimensions matching \code{mu};
#' \code{n_obs} is a positive integer; \code{max_card} is a positive integer no larger than the length of \code{mu};
#' \code{greedy_perc} is nonnegative; and \code{mve_sr} is a finite numeric scalar.
#'
#' @param mve_sr A numeric scalar representing the population MVE Sharpe ratio.
#' @param mu Numeric vector; the population mean vector.
#' @param sigma Numeric matrix; the population covariance matrix.
#' @param n_obs Integer; the sample size to simulate.
#' @param max_card Positive integer; the maximum cardinality to consider.
#' @param greedy_perc Numeric scalar (default 1.0) indicating the fraction of combinations to evaluate.
#' @param do_checks Logical; if TRUE, performs input validation (default = FALSE).
#'
#' @return A list containing the Sharpe ratio loss measures as computed by the Rcpp function.
#'
#' @examples
#' \dontrun{
#'   # Example with three assets:
#'   mu <- c(0.1, 0.2, 0.15)
#'   sigma <- diag(3)
#'   mve_sr <- 0.8
#'   n_obs <- 100
#'   max_card <- 2
#'   result <- simulate_sr_loss(mve_sr,
#'                              mu,
#'                              sigma,
#'                              n_obs,
#'                              max_card,
#'                              greedy_perc = 1.0,
#'                              do_checks = TRUE)
#'   print(result)
#' }
#'
#' @export
simulate_sr_loss <- function(mve_sr, mu, sigma, n_obs, max_card, greedy_perc = 1.0, do_checks = FALSE) {
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
    if (missing(greedy_perc) || !is.numeric(greedy_perc) || greedy_perc < 0) {
      stop("greedy_perc must be a nonnegative numeric scalar")
    }
    if (!is.numeric(mve_sr) || length(mve_sr) != 1 || !is.finite(mve_sr)) {
      stop("mve_sr must be a finite numeric scalar")
    }
  }

  .Call(`_SparsePortfolioSelection_simulate_sr_loss`, mve_sr, mu, sigma, n_obs, max_card, greedy_perc, FALSE)

}

#' Calibrate a Factor Model
#'
#' This function calibrates a factor model given a TxN matrix of asset returns and a TxK matrix of factor returns.
#' The function computes the population mean vector
#' and covariance matrix of the factors, and uses these to compute regression coefficients via
#' \eqn{\beta = (\Sigma_f)^{-1} \, \mathrm{Cov}(factors, returns)^T,}
#' where \eqn{\Sigma_f} is the covariance matrix of the factors.
#' The model-implied mean return vector is computed as
#' \eqn{\mu = \beta \, \mu_f,}
#' where \eqn{\mu_f} is the mean vector of the factors.
#' Residuals are then computed (from demeaned returns) and their variances are
#' used to form a diagonal matrix \eqn{\Sigma_0}.
#' Finally, the covariance matrix of asset returns is estimated as
#' \eqn{\Sigma = \beta \, \Sigma_f \, \beta^T + \Sigma_0.}
#'
#' The additional parameter \code{weak_coeff} adjusts the factor betas weakness by dividing them by
#' \eqn{N^{(\text{weak\_coeff}/2)}}, where \eqn{N} is the number of assets.
#' A value of \code{0} (default) indicates no adjustment,
#' while a value of \code{1} indicates full weakness.
#'
#' The parameter \code{idiosy_vol_type} determines the structure of the idiosyncratic volatility:
#' \itemize{
#'   \item \code{0} (default): homoskedastic volatility, so \eqn{\Sigma_0} is the average residual variance times the identity.
#'   \item \code{1}: heteroskedastic volatility (no correlation), so \eqn{\Sigma_0} is a diagonal matrix of asset-specific residual variances.
#' }
#'
#' @param returns A numeric matrix (T x N) of asset returns.
#' @param factors A numeric matrix (T x K) of factor returns.
#' @param weak_coeff A numeric scalar between 0 and 1 indicating the weakness of the factors.
#'                   A value of 0 (default) implies no weakness adjustment,
#'                   while a value of 1 implies full adjustment.
#' @param idiosy_vol_type A numeric scalar representing the type of idiosyncratic volatility:
#'                   0 (default) for homoskedastic; 1 for heteroskedastic.
#' @param do_checks Logical flag indicating whether to perform input validation (default is FALSE).
#'
#' @return A list with two components:
#' \describe{
#'   \item{\eqn{\mu}}{The model-implied mean return vector (an N x 1 matrix or a vector).}
#'   \item{\eqn{\Sigma}}{The model-implied covariance matrix of asset returns (an N x N matrix).}
#' }
#'
#' @examples
#' \dontrun{
#'   set.seed(123)
#'   T <- 100  # number of time periods
#'   N <- 5    # number of assets
#'   K <- 3    # number of factors
#'
#'   # Create a TxK matrix of factor returns
#'   factors <- matrix(rnorm(T * K), nrow = T, ncol = K)
#'
#'   # Create a TxN matrix of asset returns from a linear factor model plus noise
#'   beta_true <- matrix(runif(N * K), nrow = N, ncol = K)
#'   returns <- factors %*% t(beta_true) + matrix(rnorm(T * N, sd = 0.5), nrow = T, ncol = N)
#'
#'   # Calibrate the model with no weak factor adjustment (weak_coeff = 0)
#'   # under homoskedastic volatility
#'   model1 <- calibrate_factor_model(returns,
#'                                    factors,
#'                                    weak_coeff = 0,
#'                                    idiosy_vol_type = 0,
#'                                    do_checks = TRUE)
#'
#'   # Calibrate the model with moderate weak factor adjustment (weak_coeff = 0.5)
#'   # under heteroskedastic volatility
#'   model2 <- calibrate_factor_model(returns,
#'                                    factors,
#'                                    weak_coeff = 0.5,
#'                                    idiosy_vol_type = 1,
#'                                    do_checks = TRUE)
#'
#'   # Display the calibrated mean return vectors and covariance matrices
#'   print(model1$mu)
#'   print(model1$sigma)
#'   print(model2$mu)
#'   print(model2$sigma)
#' }
#'
#' @export
calibrate_factor_model <- function(returns, factors, weak_coeff = 0.0, idiosy_vol_type = 0, do_checks = FALSE) {

  # Perform input validation if do_checks is TRUE
  if (do_checks) {
    if (missing(returns) || !is.matrix(returns) || nrow(returns) == 0 || !is.numeric(returns)) {
      stop("returns must be provided and be a non-empty numeric matrix")
    }
    if (missing(factors) || !is.matrix(factors) || nrow(factors) == 0 || !is.numeric(factors)) {
      stop("factors must be provided as a non-empty numeric matrix")
    }
    if (nrow(returns) != nrow(factors)) {
      stop("The number of rows in returns must equal the number of rows in factors")
    }
    if (!is.numeric(weak_coeff) || length(weak_coeff) != 1) {
      stop("weak_coeff must be a single numeric value")
    }
    if (weak_coeff < 0 || weak_coeff > 1) {
      stop("weak_coeff must be between 0 and 1")
    }
    if (!is.numeric(idiosy_vol_type) || length(idiosy_vol_type) != 1 || !(idiosy_vol_type %in% c(0, 1))) {
      stop("idiosy_vol_type must be a numeric scalar equal to 0 or 1")
    }
  }

  # Compute the mean vector and covariance matrix of the factors
  mu_f <- colMeans(factors)
  sigma_f <- stats::cov(factors)

  # Compute the regression coefficients beta (each column corresponds to an asset)
  beta <- t(solve(sigma_f, stats::cov(factors, returns))) / ncol(returns)^(weak_coeff / 2.0)

  # Compute the model-implied mean returns
  mu <- beta %*% mu_f

  # Compute the residuals (subtract the fitted values and the asset means)
  residuals <- returns - matrix(1, nrow(returns), 1) %*% colMeans(returns) - factors %*% t(beta)

  # Compute the residual variance for each asset and create a diagonal matrix
  if (idiosy_vol_type == 0) {
    sigma0 <- mean(apply(residuals, 2, stats::var)) * diag(ncol(returns))
  } else if (idiosy_vol_type == 1) {
    sigma0 <- diag(apply(residuals, 2, stats::var))
  }

  # Construct the covariance matrix of asset returns
  sigma <- beta %*% sigma_f %*% t(beta) + sigma0

  return(list(mu = mu, sigma = sigma))
}
