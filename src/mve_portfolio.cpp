#include "mve_portfolio.h"
#include "utils.h"
#include <RcppArmadillo.h>
#include <cmath>
#include <limits>

// [[Rcpp::depends(RcppArmadillo)]]

// Compute Optimal Portfolio Sharpe Ratio (MVE)
double compute_mve_sr(const arma::vec& mu,
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
    return std::sqrt( arma::dot(mu, arma::solve(sigma, mu)) );
  }

  // Otherwise, subset the inputs according to the asset selection.
  const arma::vec mu_sel = mu.elem(selection);
  const arma::mat sigma_sel = sigma.submat(selection, selection);

  return std::sqrt( arma::dot(mu_sel, arma::solve(sigma_sel, mu_sel)) );
}

// Compute Mean-Variance Efficient Portfolio Weights
arma::vec compute_mve_weights(const arma::vec& first_moment,
                              const arma::mat& second_moment,
                              const arma::uvec& selection,
                              const double gamma,
                              const bool do_checks) {

  // Optional input checks.
  if (do_checks) {
    if (first_moment.n_elem == 0 || second_moment.n_elem == 0) {
      Rcpp::stop("First moment vector and second moment matrix must be supplied");
    }
    if (second_moment.n_rows != second_moment.n_cols) {
      Rcpp::stop("Second moment matrix must be square");
    }
    if (first_moment.n_elem != second_moment.n_rows) {
      Rcpp::stop("First moment vector and second moment matrix must have conforming dimensions");
    }
    if (selection.n_elem > 0 && arma::max(selection) > first_moment.n_elem + 1) {
      Rcpp::stop("Asset selection index out of bounds");
    }
  }

  // If the selection vector has the same length as first_moment, or if it is not supplied,
  // return the full-sample solution.
  if (selection.n_elem == first_moment.n_elem || selection.n_elem == 0) {
    return (1.0 / gamma) * arma::solve(second_moment, first_moment);
  }

  // Otherwise, subset first_moment and second_moment according to the asset selection.
  const arma::vec first_sel = first_moment.elem(selection);
  const arma::mat second_sel = second_moment.submat(selection, selection);

  // Initialize full weight vector (length = N) with zeros.
  arma::vec full_weights(first_moment.n_elem, arma::fill::zeros);

  // Place the computed weights in the positions corresponding to the selected assets.
  full_weights.elem(selection) = (1.0 / gamma) * arma::solve(second_sel, first_sel);

  return full_weights;
}

// Compute Highest Sharpe Ratio with Cardinality Constraint (sparse MVE)
Rcpp::List compute_sparse_mve_sr(const arma::vec& mu,
                                 const arma::mat& sigma,
                                 unsigned int max_card,
                                 const double greedy_perc,
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
  }

  // If max_card is 0, return empty selection.
  if (greedy_perc <= 0.0 || max_card == 0) {
    return Rcpp::List::create(Rcpp::Named("sqsr") = 0.0,
                              Rcpp::Named("selection") = arma::uvec());
  }

  // If max_card is greater than the number of assets, set it to the number of assets.
  const unsigned int n = mu.n_elem;
  if (max_card > n) {
    max_card = n;
  }

  // Initialize variables to store the best square Sharpe ratio and selection.
  double mve_sr = -std::numeric_limits<double>::infinity();
  arma::uvec mve_selection;

  // If greedy_perc is greater or equal to 1.0, evaluate all combinations for each cardinality.
  if (greedy_perc >= 1.0) {
    for (unsigned int k = 1; k <= max_card; k++) {
      // Generate all combinations of k indices from {0,1,...,n-1}.
      std::vector<arma::uvec> all_combs;
      arma::uvec current;
      // Generate combinations.
      generate_combinations(n, k, 0, current, all_combs);

      // For each combination, compute the square Sharpe ratio.
      for (size_t i = 0; i < all_combs.size(); i++) {
        // Select the current combination.
        const arma::uvec sel = all_combs[i];
        // Compute the square Sharpe ratio for the selected combination.
        const arma::vec mu_sel = mu.elem(sel);
        const arma::mat sigma_sel = sigma.submat(sel, sel);
        const double current_sr = std::sqrt( arma::dot(mu_sel, arma::solve(sigma_sel, mu_sel)) );
        // Update the maximum square Sharpe ratio and selection if current is better.
        if (current_sr > mve_sr) {
          mve_sr = current_sr;
          mve_selection = sel;
        }
      }
    }

    // Otherwise, evaluate a random sample of combinations for each cardinality.
  } else {
    for (unsigned int k = 1; k <= max_card; k++) {
      // Generate a random sample of combinations.
      const unsigned long long total_comb = nCk(n, k);
      // Calculate the number of combinations to evaluate.
      const unsigned long long num_to_eval = std::max(
        static_cast<unsigned long long>(std::floor(greedy_perc * total_comb)),
        static_cast<unsigned long long>(1));

      for (unsigned long long i = 0; i < num_to_eval; i++) {
        // Generate a random combination of k distinct indices from 0 to n-1.
        const arma::uvec sel = random_combination(n, k);
        const arma::vec mu_sel = mu.elem(sel);
        const arma::mat sigma_sel = sigma.submat(sel, sel);
        // Compute the square Sharpe ratio for the selected combination.
        const double current_sr = std::sqrt( arma::dot(mu_sel, arma::solve(sigma_sel, mu_sel)) );
        // Update the maximum square Sharpe ratio and selection if current is better.
        if (current_sr > mve_sr) {
          mve_sr = current_sr;
          mve_selection = sel;
        }
      }
    }
  }

  // Return the results as a list.
  return Rcpp::List::create(Rcpp::Named("sr") = mve_sr,
                            Rcpp::Named("selection") = mve_selection);
}
