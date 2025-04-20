#ifndef MVE_PORTFOLIO_H
#define MVE_PORTFOLIO_H

#include <RcppArmadillo.h>

//' Compute Portfolio Sharpe Ratio
//'
//' Given a vector of portfolio weights, a vector of mean returns, and a covariance matrix,
//' this function computes the Sharpe ratio defined as
//' \deqn{\frac{w^T \mu}{\sqrt{w^T \Sigma w}}.}
//'
//' @param weights A numeric vector of portfolio weights.
//' @param mu A numeric vector of expected returns.
//' @param sigma A numeric covariance matrix.
//' @param do_checks Logical flag indicating whether to perform input checks (default = false).
//' @return A double representing the Sharpe ratio.
// [[Rcpp::export]]
double compute_sr_cpp(const arma::vec& weights,
                      const arma::vec& mu,
                      const arma::mat& sigma,
                      const bool do_checks = false);

//' Compute Mean-Variance Efficient Portfolio Sharpe Ratio
//'
//' Given a mean vector \eqn{\mu}, a covariance matrix \eqn{\Sigma},
//' and an asset selection vector, this function computes the optimal Sharpe ratio defined as
//' \deqn{\sqrt{\mu^T \Sigma^{-1}\mu}}
//' over the selected assets.
//' If the provided asset selection vector has length less than N,
//' the function uses the subset of assets and computes the ratio using those assets.
//'
//' @param mu Mean vector.
//' @param sigma Covariance matrix.
//' @param selection Index vector with asset indices. Default is an empty
//'        vector, which means all assets are selected.
//' @param do_checks Logical flag indicating whether to perform input checks (default = FALSE).
//' @return A scalar value corresponding to \eqn{\sqrt{\mu^T \Sigma^{-1}\mu}}.
// [[Rcpp::export]]
double compute_mve_sr_cpp(const arma::vec& mu,
                          const arma::mat& sigma,
                          const arma::uvec& selection = arma::uvec(),
                          const bool do_checks = false);

//' Compute Mean-Variance Efficient (MVE) Portfolio Weights
//'
//' Given the risk-aversion parameter \eqn{\gamma}, the first moment vector, the covariance matrix,
//' and the selection index vector, this function computes the mean variance efficient portfolio weights
//' \deqn{w = \frac{1}{\gamma}\, \text{sigma}^{-1}\, \text{mu}},
//' over the selected assets.
//' If the provided asset selection vector has length less than N,
//' the function returns an N-length weight vector with zeros for the unselected assets.
//'
//' @param mu First moment vector.
//' @param sigma Covariance matrix.
//' @param selection Unsigned integer vector with asset indices. Default is an empty
//'        vector, which means all assets are selected.
//' @param gamma Risk aversion parameter. Default is 1.
//' @param do_checks Logical flag indicating whether to perform input checks (default = FALSE).
//' @return An N-length vector of mean variance efficient weights.
// [[Rcpp::export]]
arma::vec compute_mve_weights_cpp(const arma::vec& mu,
                                  const arma::mat& sigma,
                                  const arma::uvec& selection = arma::uvec(),
                                  const double gamma = 1.0,
                                  const bool do_checks = false);

//' Compute Mean-Variance Efficient Sharpe Ratio with Cardinality K
//'
//' This function takes as inputs the mean vector \code{mu}, the covariance matrix \code{sigma},
//' the maximum active portfolio cardinality \code{max_card},
//' the maximum number of combinations per cardinality to evaluate \code{max_comb},
//' and the risk-aversion parameter \eqn{\gamma}.
//' With these inputs, it searches over all combinations of assets with cardinality from 1 up to \code{max_card}
//' and computes the Sharpe ratio defined as
//' \eqn{\mu^T \, \sigma^{-1}\, \mu},
//' and computes the mean-variance efficient weights
//' \deqn{w = \frac{1}{\gamma}\, \text{second\_moment}^{-1}\, \text{first\_moment}},
//' over the selected assets.
//' It returns the highest Sharpe ratio found along with the associated weights and asset selection.
//'
//' @param mu Mean vector.
//' @param sigma Coveriance matrix.
//' @param max_card Maximum investment cardinality (from 1 up to the number of assets).
//' @param max_comb Maximum number of combinations to consider. If 0 (default),
//'                 all combinations are computed.
//' @param gamma Risk aversion parameter. Default is 1.
//' @param do_checks Logical flag indicating whether to perform input checks (default = \code{FALSE}).
//' @return A list with \code{sr} (the optimal Sharpe ratio), \code{mve_weights}
//'         (the optimal weights) and \code{selection} (the optimal asset selection).
// [[Rcpp::export]]
Rcpp::List compute_mve_sr_cardk_cpp(const arma::vec& mu,
                                     const arma::mat& sigma,
                                     const unsigned int max_card,
                                     const unsigned int max_comb = 0,
                                     const double gamma = 1.0,
                                     const bool do_checks = false);

// Compute the Mean-Variance Efficient Sharpe Ratio decomposition in Estimation and Selection term
//
// This function computes the sample MVE Sharpe ratio, the sample MVE Sharpe ratio
// with maximum cardinality k, and the population MVE Sharpe ratio with maximum
// cardinality k decomposed in estimation and selection component, taking as input
// the population parameters (\code{mu} and \code{sigma})
// and corresponding sample parameters (\code{mu_sample} and \code{sigma_sample}),
// together with a maximum cardinality (\code{max_card}) and a maximum number of
// combinations (\code{max_comb}).
//
// @param mu Mean vector.
// @param sigma Covariance matrix.
// @param mu_hat Sample mean vector.
// @param sigma Sample covariance matrix.
// @param max_card Maximum cardinality to consider (from \code{1} up to the number of assets).
// @param max_comb Maximum number of combinations to consider. If \code{0{} (default),
//                 all combinations are computed.
// @param do_checks Logical flag indicating whether to perform input checks (default = \code{FALSE}).
// @return A list with:
//    - \code{sample_mve_sr} computed as the optimal sample mve Sharpe ratio,
//    - \code{sample_mve_sr_cardk} computed as the optimal sample mve Sharpe ratio
//      with cardinality \code{max_card},
//    - \code{mve_sr_cardk_est_term} computed as \eqn{w^T \mu^T / \sqrt{w^T\sigma w}}
//      where \code{w} are the optimal sample mve weights,
//    - \code{mve_sr_cardk_sel_term} computed as \eqn{\mu_S^T  \sigma_S^{-1} \mu_S}
//      where \code{S} is the set of assets yielding the optimal sample mve Sharpe ratio.
Rcpp::List compute_mve_sr_decomposition_cpp(const arma::vec& mu,
                                            const arma::mat& sigma,
                                            const arma::vec& mu_sample,
                                            const arma::mat& sigma_sample,
                                            const unsigned int max_card,
                                            const unsigned int max_comb = 0,
                                            const bool do_checks = false);

//' Simulate Sharpe Ratio Loss
//'
//' This function simulates a sample of size \code{n_obs} from a multivariate normal
//' distribution with mean vector \code{mu} and covariance matrix \code{sigma}.
//' It computes the sample mean vector (\code{mu_sample}) and the sample covariance matrix (\code{sigma_sample}),
//' then calls \code{compute_sr_sparsity_loss} with the population and sample parameters,
//' the maximum cardinality (\code{max_card}), and the maximum number of combinations (\code{max_comb}).
//'
//' @param mu A numeric vector; the population mean vector.
//' @param sigma A numeric matrix; the population covariance matrix.
//' @param n_obs An integer specifying the sample size to simulate.
//' @param max_card Maximum cardinality to consider (from 1 up to the number of assets).
//' @param max_comb Maximum number of combinations to consider. If 0 (default),
//'                 all combinations are computed.
//' @param do_checks Logical; if TRUE, input checks are performed.
//'
//' @return // @return A list with:
//'    - \code{sample_mve_sr} computed as the optimal sample mve Sharpe ratio,
//'    - \code{sample_mve_sr_cardk} computed as the optimal sample mve Sharpe ratio
//'      with cardinality \code{max_card},
//'    - \code{mve_sr_cardk_est_term} computed as \eqn{w^T \mu^T / \sqrt{w^T\sigma w}}
//'      where \code{w} are the optimal sample mve weights,
//'    - \code{mve_sr_cardk_sel_term} computed as \eqn{\mu_S^T  \sigma_S^{-1} \mu_S}
//'      where \code{S} is the set of assets yielding the optimal sample mve Sharpe ratio.
// [[Rcpp::export]]
Rcpp::List simulate_mve_sr_cpp(const arma::vec& mu,
                               const arma::mat& sigma,
                               const unsigned int n_obs,
                               const unsigned int max_card,
                               const unsigned int max_comb = 0,
                               const bool do_checks = false);

#endif // MVE_PORTFOLIO_H
