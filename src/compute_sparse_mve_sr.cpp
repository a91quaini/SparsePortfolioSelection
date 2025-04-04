#include <RcppArmadillo.h>
#include <vector>
#include <algorithm>
#include <random>
#include <cmath>
#include <limits>

// [[Rcpp::depends(RcppArmadillo)]]

// Helper: Compute nCk (combinations count) as an unsigned long long.
unsigned long long nCk(const unsigned int n, const unsigned int k) {
  if (k > n) return 0;
  if (k == 0 || k == n) return 1;
  unsigned long long res = 1;
  for (unsigned int i = 1; i <= k; i++) {
    res = res * (n - i + 1) / i;
  }
  return res;
}

// Helper: recursively generate all combinations of k indices from {0,1,...,n-1}.
void generate_combinations(const unsigned int n, const unsigned int k, const unsigned int offset,
                           arma::uvec &current, std::vector<arma::uvec> &all) {
  if (current.n_elem == k) {
    all.push_back(current);
    return;
  }
  for (unsigned int i = offset; i < n; i++) {
    current.insert_rows(current.n_elem, 1);
    current(current.n_elem - 1) = i;
    generate_combinations(n, k, i + 1, current, all);
    current.shed_row(current.n_elem - 1);
  }
}

// Helper: generate a random combination of k distinct indices from 0 to n-1.
arma::uvec random_combination(const unsigned int n, const unsigned int k) {
  std::vector<unsigned int> indices(n);
  for (unsigned int i = 0; i < n; i++) {
    indices[i] = i;
  }
  std::random_device rd;
  std::mt19937 g(rd());
  std::shuffle(indices.begin(), indices.end(), g);
  std::sort(indices.begin(), indices.begin() + k);
  arma::uvec comb(k);
  for (unsigned int i = 0; i < k; i++) {
    comb(i) = indices[i];
  }
  return comb;
}

//' Compute Highest Sharpe Ratio with Cardinality Constraint
//'
//' This function takes as inputs the mean vector \code{mu}, the covariance matrix \code{sigma},
//' the maximum active portfolio cardinality \code{max_card},
//' and the fraction of combinations to evaluate \code{greedy_perc}.
//' With these inputs, it searches over all combinations of assets with cardinality from 1 up to \code{max_card}
//' and computes the square Sharpe ratio defined as
//' \eqn{\mu^T \, \sigma^{-1}\, \mu}.
//' It returns the highest square Sharpe ratio found along with the associated asset selection.
//' If \code{greedy_perc} is less than 1, then for each cardinality the search is performed over a random
//' subset (a fraction equal to \code{greedy_perc}) of all possible combinations.
//'
//' @param mu Mean vector.
//' @param sigma Coveriance matrix.
//' @param max_card Maximum cardinality to consider (from 1 up to the number of assets).
//' @param greedy_perc If less than 1, the fraction of combinations to evaluate for each cardinality;
//'                    if 1 or greater, all combinations are evaluated;
//'                    if less than 0, no combinations are evaluated.
//' @param do_checks Logical flag indicating whether to perform input checks (default = FALSE).
//' @return A list with \code{sqsr} (the optimal square Sharpe ratio) and \code{selection} (the asset indices of the optimal selection).
//' @export
// [[Rcpp::export]]
Rcpp::List compute_sparse_mve_sr(const arma::vec& mu,
                           const arma::mat& sigma,
                           unsigned int max_card = 1,
                           const double greedy_perc = 1.0,
                           const bool do_checks = false) {
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
       const arma::uvec sel = random_combination(n, k);
       const arma::vec mu_sel = mu.elem(sel);
       const arma::mat sigma_sel = sigma.submat(sel, sel);
       const double current_sr = std::sqrt( arma::dot(mu_sel, arma::solve(sigma_sel, mu_sel)) );
       if (current_sr > mve_sr) {
         mve_sr = current_sr;
         mve_selection = sel;
       }
     }
   }
 }

 return Rcpp::List::create(Rcpp::Named("sr") = mve_sr,
                           Rcpp::Named("selection") = mve_selection);
}
