#include "mve_portfolio.h"
#include "utils.h"
#include <RcppArmadillo.h>
#include <cmath>
#include <limits>

// [[Rcpp::depends(RcppArmadillo)]]

////////////////////////////////////////////////////////////////////////////////
// Compute Portfolio Sharpe Ratio: (w^T mu) / sqrt(w^T sigma w)
double compute_sr_cpp(const arma::vec& weights,
                  const arma::vec& mu,
                  const arma::mat& sigma,
                  const bool do_checks) {

  // Optional input checks.
  if (do_checks) {
    if (weights.n_elem == 0) {
      Rcpp::stop("weights must be non-empty");
    }
    if (mu.n_elem == 0) {
      Rcpp::stop("mu must be non-empty");
    }
    if (sigma.n_elem == 0) {
      Rcpp::stop("sigma must be provided");
    }
    if (weights.n_elem != mu.n_elem) {
      Rcpp::stop("weights and mu must be of the same length");
    }
    if (sigma.n_rows != sigma.n_cols) {
      Rcpp::stop("sigma must be a square matrix");
    }
    if (sigma.n_rows != weights.n_elem) {
      Rcpp::stop("The dimensions of sigma must match the length of weights");
    }
  }

  // Compute the Sharpe ratio: (w^T mu) / sqrt(w^T sigma w)
  return arma::dot(weights, mu) / std::sqrt( arma::dot(weights, sigma * weights) );
}

////////////////////////////////////////////////////////////////////////////////
// Compute Optimal Portfolio Sharpe Ratio (MVE)
double compute_mve_sr_cpp(const arma::vec& mu,
                          const arma::mat& sigma,
                          const arma::uvec& selection,
                          const bool do_checks) {

  // Optional input checks.
  if (do_checks) {
    if (mu.n_elem == 0 || sigma.n_elem == 0) {
      Rcpp::stop("Mean vector mu and covariance matrix sigma must be supplied");
    }
    if (sigma.n_rows != sigma.n_cols) {
      Rcpp::stop("Covariance matrix sigma must be square");
    }
    if (mu.n_elem != sigma.n_rows) {
      Rcpp::stop("Mean vector mu must have a length equal to the dimensions of sigma");
    }
    if (selection.n_elem > 0 && arma::max(selection) >= mu.n_elem) {
      Rcpp::stop("Asset selection index out of bounds");
    }
  }

  // If selection is not supplied or is full, use full mu and sigma.
  if (selection.n_elem == 0 || selection.n_elem == mu.n_elem) {
    //return std::sqrt( arma::dot(mu, arma::solve(sigma, mu)) );
    return std::sqrt( arma::dot(mu, solve_sympd(sigma, mu)) );
  }

  // Otherwise, subset the inputs according to the asset selection.
  const arma::vec mu_sel = mu.elem(selection);
  // return std::sqrt( arma::dot(mu_sel, arma::solve(
  //     sigma.submat(selection, selection), mu_sel)) );
  return std::sqrt( arma::dot(mu_sel, solve_sympd(
      sigma.submat(selection, selection), mu_sel)) );
}

///////////////////////////////////////////////////////////////////////////////

// Compute Mean-Variance Efficient Portfolio Weights
arma::vec compute_mve_weights_cpp(const arma::vec& mu,
                                  const arma::mat& second_moment,
                                  const arma::uvec& selection,
                                  const double gamma,
                                  const bool do_checks) {

  // Optional input checks.
  if (do_checks) {
    if (mu.n_elem == 0 || second_moment.n_elem == 0) {
      Rcpp::stop("First moment vector and second moment matrix must be supplied");
    }
    if (second_moment.n_rows != second_moment.n_cols) {
      Rcpp::stop("Second moment matrix must be square");
    }
    if (mu.n_elem != second_moment.n_rows) {
      Rcpp::stop("First moment vector and second moment matrix must have conforming dimensions");
    }
    if (selection.n_elem > 0 && arma::max(selection) > mu.n_elem + 1) {
      Rcpp::stop("Asset selection index out of bounds");
    }
  }

  // If the selection vector has the same length as mu, or if it is not supplied,
  // return the full-sample solution.
  if (selection.n_elem == 0 || selection.n_elem == mu.n_elem) {
    // return (1.0 / gamma) * arma::solve(second_moment, mu);
    return (1.0 / gamma) * solve_sympd(second_moment, mu);
  }

  // Initialize full weight vector (length = N) with zeros.
  arma::vec full_weights(mu.n_elem, arma::fill::zeros);

  // Place the computed weights in the positions corresponding to the selected assets.
  // full_weights.elem(selection) = (1.0 / gamma) * arma::solve(
  //   second_moment.submat(selection, selection), mu.elem(selection));
  full_weights.elem(selection) = (1.0 / gamma) * solve_sympd(
    second_moment.submat(selection, selection), mu.elem(selection));

  return full_weights;
}

///////////////////////////////////////////////////////////////////////////////

// Compute Highest Sharpe Ratio with Cardinality K
Rcpp::List compute_mve_sr_cardk_cpp(const arma::vec& mu,
                                    const arma::mat& sigma,
                                    const unsigned int max_card,
                                    const unsigned int max_comb,
                                    const double gamma,
                                    const bool do_checks) {
  // Input checks.
  if (do_checks) {
    if (mu.n_elem == 0 || sigma.n_elem == 0) {
      Rcpp::stop("Mean vector and covariance matrix must be supplied");
    }
    if (sigma.n_rows != sigma.n_cols) {
      Rcpp::stop("Covariance matrix must be square");
    }
    if (mu.n_elem != sigma.n_rows) {
      Rcpp::stop("Mean vector length must equal the dimensions of the covariance matrix");
    }
    if (max_card < 1 || max_card > mu.n_elem) {
      Rcpp::stop("max_card must be between 1 and the number of assets");
    }
    if (max_comb < 0) {
      Rcpp::stop("max_comb must be non-negative");
    }
  }

  // Set number of assets
  const unsigned int n = mu.n_elem;

  // Initialize variables to store the best Sharpe ratio and selection.
  double mve_sr = 0.0;
  arma::uvec mve_selection;

  // If max_comb is equal to 0, evaluate all combinations for each cardinality.
  if (max_comb == 0) {
    for (unsigned int k = 1; k <= max_card; k++) {
      // Generate all combinations of k indices from {0,1,...,n-1}.
      std::vector<arma::uvec> all_combs;
      arma::uvec current;
      // Generate combinations.
      generate_combinations(n, k, 0, current, all_combs);

      // For each combination, compute the Sharpe ratio.
      for (size_t i = 0; i < all_combs.size(); i++) {
        // Select the current combination.
        const arma::uvec selection = all_combs[i];
        // Compute the Sharpe ratio for the selected combination.
        const double current_sr = compute_mve_sr_cpp(mu, sigma, selection, false);
        // Update the maximum Sharpe ratio and selection if current is better.
        if (current_sr > mve_sr) {
          mve_sr = current_sr;
          mve_selection = selection;
        }
      }
    }

    // Otherwise, evaluate a random sample of combinations for each cardinality.
  } else {
    // Set indices from 0 to n-1.
    const arma::uvec indices = arma::regspace<arma::uvec>(0, n - 1);

    for (unsigned int k = 1; k <= max_card; k++) {
      for (unsigned int i = 0; i < max_comb; i++) {
        // Generate a random combination of k distinct indices from 0 to n-1.
        const arma::uvec shuffled_indices = arma::shuffle(indices);
        const arma::uvec selection_k = shuffled_indices.head(k);
        const double current_sr = compute_mve_sr_cpp(mu, sigma, selection_k, false);
        // Update the maximum Sharpe ratio and selection if current is better.
        if (current_sr > mve_sr) {
          mve_sr = current_sr;
          mve_selection = selection_k;
        }
      }
    }
  }

  // Compute the MVE weights using the selected assets.
  const arma::vec mve_weights = compute_mve_weights_cpp(mu, sigma + mu * mu.t(), mve_selection, gamma, false);

  // Return the results as a list.
  return Rcpp::List::create(Rcpp::Named("sr") = mve_sr,
                            Rcpp::Named("weights") = mve_weights,
                            Rcpp::Named("selection") = mve_selection);
}

///////////////////////////////////////////////////////////////////////////////

// Compute MVE Sharpe Ratio Decomposition
Rcpp::List compute_mve_sr_decomposition_cpp(const arma::vec& mu,
                                            const arma::mat& sigma,
                                            const arma::vec& mu_sample,
                                            const arma::mat& sigma_sample,
                                            const unsigned int max_card,
                                            const unsigned int max_comb,
                                            const bool do_checks) {

  // Check if the input parameters are valid
  if (do_checks) {
    if (mu.n_elem == 0) {
      Rcpp::stop("mu must be non-empty");
    }
    if (sigma.n_elem == 0) {
      Rcpp::stop("sigma must be provided");
    }
    if (sigma.n_rows != sigma.n_cols) {
      Rcpp::stop("sigma must be a square matrix");
    }
    if (mu.n_elem != sigma.n_rows) {
      Rcpp::stop("Length of mu must equal number of rows of sigma");
    }
    if (mu_sample.n_elem == 0) {
      Rcpp::stop("mu_sample must be non-empty");
    }
    if (sigma_sample.n_elem == 0) {
      Rcpp::stop("sigma_sample must be provided");
    }
    if (sigma_sample.n_rows != sigma_sample.n_cols) {
      Rcpp::stop("sigma_sample must be a square matrix");
    }
    if (mu_sample.n_elem != sigma_sample.n_rows) {
      Rcpp::stop("Length of mu_sample must equal number of rows of sigma_sample");
    }
    if (max_card < 1 || max_card > mu.n_elem) {
      Rcpp::stop("max_card must be between 1 and the number of assets");
    }
    if (max_comb < 0) {
      Rcpp::stop("max_comb must be non-negative");
    }
  }

  // Compute the sample MVE Sharpe ratio
  const double sample_mve = compute_mve_sr_cpp(mu_sample,
                                               sigma_sample,
                                               arma::uvec(),
                                               false);

  // Compute the sample MVE Sharpe ratio with max cardinality k
  const Rcpp::List sample_mve_cardk = compute_mve_sr_cardk_cpp(mu_sample,
                                                               sigma_sample,
                                                               max_card,
                                                               max_comb,
                                                               1.0,
                                                               false);

  // 1. Compute the population MVE Sharpe ratio estimation term
  const double mve_sr_cardk_est_term = compute_sr_cpp(sample_mve_cardk["weights"],
                                                      mu,
                                                      sigma,
                                                      false);

  // 2. Compute the population MVE Sharpe ratio selection term
  const double mve_sr_cardk_sel_term = compute_mve_sr_cpp(mu,
                                                          sigma,
                                                          sample_mve_cardk["selection"],
                                                          false);

  // Compute the Sharpe ratio loss
  return Rcpp::List ::create(Rcpp::Named("sample_mve_sr") = sample_mve,
                             Rcpp::Named("sample_mve_sr_cardk") = sample_mve_cardk["sr"],
                             Rcpp::Named("mve_sr_cardk_est_term") = mve_sr_cardk_est_term,
                             Rcpp::Named("mve_sr_cardk_sel_term") = mve_sr_cardk_sel_term);

}

////////////////////////////////////////////////////////////////////////////////

Rcpp::List simulate_mve_sr_cpp(const arma::vec& mu,
                               const arma::mat& sigma,
                               const unsigned int n_obs,
                               const unsigned int max_card,
                               const unsigned int max_comb,
                               const bool do_checks) {

  // Check if the input parameters are valid
  if (do_checks) {
    if (mu.n_elem == 0) {
      Rcpp::stop("mu must be non-empty");
    }
    if (sigma.n_elem == 0) {
      Rcpp::stop("sigma must be provided");
    }
    if (sigma.n_rows != sigma.n_cols) {
      Rcpp::stop("sigma must be a square matrix");
    }
    if (mu.n_elem != sigma.n_rows) {
      Rcpp::stop("Length of mu must equal number of rows of sigma");
    }
    if (max_card < 1 || max_card > mu.n_elem) {
      Rcpp::stop("max_card must be between 1 and the number of assets");
    }
    if (max_comb < 0) {
      Rcpp::stop("max_comb must be non-negative");
    }
  }

 // Simulate sample data from a multivariate normal with mean mu and covariance sigma.
 const arma::mat sample = arma::mvnrnd(mu, sigma, n_obs);

 // Compute the sample mean (as a column vector)
 const arma::vec mu_sample = arma::mean(sample, 1);

 // Compute the sample covariance matrix.
 const arma::mat sigma_sample = arma::cov(sample.t());

 // Compute the MVE Sharpe ratio decomposition
 return compute_mve_sr_decomposition_cpp(mu, sigma, mu_sample, sigma_sample, max_card, max_comb, false);
}
