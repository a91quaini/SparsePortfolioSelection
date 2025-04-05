#' Compute Sharpe Ratio Loss due to Estimation and Selection Errors
#'
#' This function takes as input a population MVE Sharpe ratio (\eqn{mve\_sr}),
#' the population parameters (\code{mu} and \code{sigma}) and the corresponding sample parameters
#' (\code{mu_hat} and \code{sigma_hat}), together with a maximum cardinality (\code{max_card})
#' and a greedy percentage (\code{greedy_perc}). It computes the sample sparse MVE portfolio,
#' then uses it to compute two versions of the population Sharpe ratio: one computed over the full
#' population using the sample selection and one computed on the sample.
#' Finally, it returns a list containing:
#'   sr_loss: The loss between the population MVE Sharpe ratio and the sample portfolio SR.
#'   sr_loss_selection: The loss between the population MVE SR and the population SR computed using the sample selection.
#'   sr_loss_estimation: The difference between the two computed population SRs.
#'
#' When \code{do_checks} is TRUE, the wrapper performs input validation to ensure:
#' \itemize{
#'   \item Population parameters \code{mu} and \code{sigma} are non-empty, numeric, and that \code{sigma} is square with dimensions matching \code{mu}.
#'   \item Sample parameters \code{mu_hat} and \code{sigma_hat} are non-empty, numeric, and that \code{sigma_hat} is square with dimensions matching \code{mu_hat}.
#'   \item \code{max_card} is a positive integer no larger than \code{length(mu_hat)}.
#'   \item \code{greedy_perc} is nonnegative.
#'   \item \code{mve_sr} is a finite number.
#' }
#'
#' @param mve_sr A numeric scalar representing the population MVE Sharpe ratio.
#' @param mu Numeric vector; the population mean vector.
#' @param sigma Numeric matrix; the population covariance matrix.
#' @param mu_hat Numeric vector; the sample mean vector.
#' @param sigma_hat Numeric matrix; the sample covariance matrix.
#' @param max_card Positive integer; the maximum cardinality to consider (from 1 up to the number of assets in the sample).
#' @param greedy_perc Numeric; if less than 1, the fraction of combinations to evaluate for each cardinality;
#'                    if 1 or greater, all combinations are evaluated.
#' @param do_checks Logical flag indicating whether to perform input checks (default = FALSE).
#' @return A list with three elements:
#'   sr_loss: Loss between the population MVE Sharpe ratio and the sample portfolio Sharpe ratio.
#'   sr_loss_selection: Loss between the population MVE Sharpe ratio and the population Sharpe ratio computed using the sample selection.
#'   sr_loss_estimation: Difference between the two computed population Sharpe ratios.
#' @examples
#' # Example with three assets:
#' mu <- c(0.1, 0.2, 0.15)
#' sigma <- diag(3)
#' mu_hat <- c(0.12, 0.18, 0.16)
#' sigma_hat <- diag(3) * 1.1
#' mve_sr <- 0.8  # Example population MVE Sharpe ratio
#' max_card <- 2
#' result <- compute_sr_sparsity_loss(mve_sr,
#'                                    mu,
#'                                    sigma,
#'                                    mu_hat,
#'                                    sigma_hat,
#'                                    max_card,
#'                                    greedy_perc = 1.0,
#'                                    do_checks = TRUE)
#' print(result)
#' @export
compute_sr_sparsity_loss <- function(mve_sr, mu, sigma, mu_hat, sigma_hat, max_card, greedy_perc = 1.0, do_checks = FALSE) {

  # Optinal check of the inputs.
  if (do_checks) {
    # Validate population parameters.
    if (missing(mu) || length(mu) == 0) {
      stop("mu must be provided and non-empty")
    }
    if (!is.numeric(mu)) {
      stop("mu must be numeric")
    }
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

    # Validate sample parameters.
    if (missing(mu_hat) || length(mu_hat) == 0) {
      stop("mu_hat must be provided and non-empty")
    }
    if (!is.numeric(mu_hat)) {
      stop("mu_hat must be numeric")
    }
    if (missing(sigma_hat) || !is.matrix(sigma_hat) || nrow(sigma_hat) == 0) {
      stop("sigma_hat must be provided as a non-empty matrix")
    }
    if (!is.numeric(sigma_hat)) {
      stop("sigma_hat must be numeric")
    }
    if (nrow(sigma_hat) != ncol(sigma_hat)) {
      stop("sigma_hat must be a square matrix")
    }
    if (length(mu_hat) != nrow(sigma_hat)) {
      stop("The length of mu_hat must equal the number of rows of sigma_hat")
    }

    # Validate max_card.
    if (missing(max_card) || length(max_card) == 0) {
      stop("max_card must be provided")
    }
    if (!is.numeric(max_card) || max_card < 1) {
      stop("max_card must be a positive integer")
    }
    max_card <- as.integer(max_card)
    if (max_card > length(mu_hat)) {
      stop("max_card cannot exceed the number of assets (length of mu_hat)")
    }

    # Validate greedy_perc.
    if (missing(greedy_perc) || length(greedy_perc) == 0) {
      stop("greedy_perc must be provided")
    }
    if (!is.numeric(greedy_perc)) {
      stop("greedy_perc must be numeric")
    }
    if (greedy_perc < 0) {
      stop("greedy_perc must be non-negative")
    }

    # Validate mve_sr.
    if (!is.numeric(mve_sr) || length(mve_sr) != 1 || !is.finite(mve_sr)) {
      stop("mve_sr must be a finite numeric scalar")
    }
  }

  .Call(`_SparsePortfolioSelection_compute_sr_sparsity_loss`,
        mve_sr,
        mu,
        sigma,
        mu_hat,
        sigma_hat,
        as.integer(max_card),
        greedy_perc,
        FALSE)
}
