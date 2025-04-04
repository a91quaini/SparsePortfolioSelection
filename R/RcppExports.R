# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Compute Mean-Variance Efficient Portfolio Sharpe Ratio
#'
#' Given a mean vector \eqn{\mu}, a covariance matrix \eqn{\Sigma},
#' and an asset selection vector, this function computes the optimal Sharpe ratio defined as
#' \deqn{\sqrt{\mu^T \Sigma^{-1}\mu}}
#' over the selected assets.
#' If the provided asset selection vector has length less than N,
#' the function uses the subset of assets and computes the ratio using those assets.
#'
#' @param mu Mean vector.
#' @param sigma Covariance matrix.
#' @param selection Unsigned integer vector with asset indices.
#' @param do_checks Logical flag indicating whether to perform input checks (default = FALSE).
#' @return A scalar value corresponding to \eqn{\sqrt{\mu^T \Sigma^{-1}\mu}}.
#' @export
compute_mve_sr <- function(mu, sigma, selection, do_checks = FALSE) {
    .Call(`_SparsePortfolioSelection_compute_mve_sr`, mu, sigma, selection, do_checks)
}

#' Compute Mean-Variance Efficient (MVE) Portfolio Weights
#'
#' Given a risk-aversion parameter \eqn{\gamma}, a first moment vector, a second moment matrix,
#' and a selection index vector, this function computes the mean variance efficient portfolio weights
#' \deqn{w = \frac{1}{\gamma}\, \text{second\_moment}^{-1}\, \text{first\_moment}},
#' over the selected assets.
#' If the provided asset selection vector has length less than N,
#' the function returns an N-length weight vector with zeros for the unselected assets.
#'
#' @param mu First moment vector.
#' @param second_moment Second moment matrix.
#' @param selection Unsigned integer vector with asset indices.
#' @param gamma Risk aversion parameter. Default is 1.
#' @param do_checks Logical flag indicating whether to perform input checks (default = FALSE).
#' @return An N-length vector of mean variance efficient weights.
#' @export
compute_mve_weights <- function(mu, second_moment, selection, gamma = 1.0, do_checks = FALSE) {
    .Call(`_SparsePortfolioSelection_compute_mve_weights`, mu, second_moment, selection, gamma, do_checks)
}

#' Compute Mean-Variance Efficient Sharpe Ratio with Cardinality Constraint
#'
#' This function takes as inputs the mean vector \code{mu}, the covariance matrix \code{sigma},
#' the maximum active portfolio cardinality \code{max_card},
#' and the fraction of combinations to evaluate \code{greedy_perc}.
#' With these inputs, it searches over all combinations of assets with cardinality from 1 up to \code{max_card}
#' and computes the square Sharpe ratio defined as
#' \eqn{\mu^T \, \sigma^{-1}\, \mu}.
#' It returns the highest square Sharpe ratio found along with the associated asset selection.
#' If \code{greedy_perc} is less than 1, then for each cardinality the search is performed over a random
#' subset (a fraction equal to \code{greedy_perc}) of all possible combinations.
#'
#' @param mu Mean vector.
#' @param sigma Coveriance matrix.
#' @param max_card Maximum cardinality to consider (from 1 up to the number of assets).
#' @param greedy_perc If less than 1, the fraction of combinations to evaluate for each cardinality;
#'                    if 1 or greater, all combinations are evaluated;
#'                    if less than 0, no combinations are evaluated.
#' @param do_checks Logical flag indicating whether to perform input checks (default = FALSE).
#' @return A list with \code{sqsr} (the optimal square Sharpe ratio) and \code{selection} (the asset indices of the optimal selection).
#' @export
compute_sparse_mve_sr <- function(mu, sigma, max_card = 1L, greedy_perc = 1.0, do_checks = FALSE) {
    .Call(`_SparsePortfolioSelection_compute_sparse_mve_sr`, mu, sigma, max_card, greedy_perc, do_checks)
}

#' Compute Sharpe Ratio Loss due to Estimation and Selection Errors
#'
#' This function takes as input a population MVE Sharpe ratio (mve_sr),
#' the population parameters (mu and sigma) and corresponding sample parameters
#' (mu_sample and sigma_sample), together with a maximum cardinality (max_card)
#' and a greedy percentage (greedy_perc). It computes the sample sparse MVE portfolio,
#' then uses it to compute two versions of the population Sharpe ratio (one computed on the
#' full population using the sample selection and one computed on the sample).
#' Finally, it returns a list containing:
#'   sr_loss: The loss between the population MVE Sharpe ratio and the sample portfolio SR.
#'   sr_loss_selection: The loss between the population MVE SR and the population SR computed using the sample selection.
#'   sr_loss_estimation: The difference between the two computed population SRs.
#'
#' @param mve_sr Optimal population Sharpe ratio.
#' @param mu Mean vector.
#' @param sigma Covariance matrix.
#' @param mu_hat Sample mean vector.
#' @param sigma Sample covariance matrix.
#' @param max_card Maximum cardinality to consider (from 1 up to the number of assets).
#' @param greedy_perc If less than 1, the fraction of combinations to evaluate for each cardinality;
#'                    if 1 or greater, all combinations are evaluated;
#'                    if less than 0, no combinations are evaluated.
#' @param do_checks Logical flag indicating whether to perform input checks (default = FALSE).
#' @return A scalar value corresponding to \eqn{\sqrt{\mu^T \Sigma^{-1}\mu}}.
#' @export
compute_sr_sparsity_loss <- function(mve_sr, mu, sigma, mu_sample, sigma_sample, max_card, greedy_perc = 1.0, do_checks = FALSE) {
    .Call(`_SparsePortfolioSelection_compute_sr_sparsity_loss`, mve_sr, mu, sigma, mu_sample, sigma_sample, max_card, greedy_perc, do_checks)
}

#' Compute Portfolio Sharpe Ratio
#'
#' Given a vector of portfolio weights, a vector of mean returns, and a covariance matrix,
#' this function computes the Sharpe ratio defined as
#' \deqn{\frac{w^T \mu}{\sqrt{w^T \Sigma w}}.}
#'
#' @param weights A numeric vector of portfolio weights.
#' @param mu A numeric vector of expected returns.
#' @param sigma A numeric covariance matrix.
#' @param do_checks Logical flag indicating whether to perform input checks (default = false).
#' @return A double representing the Sharpe ratio.
#' @export
compute_sr <- function(weights, mu, sigma, do_checks = FALSE) {
    .Call(`_SparsePortfolioSelection_compute_sr`, weights, mu, sigma, do_checks)
}

