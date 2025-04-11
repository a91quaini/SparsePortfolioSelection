// Author: Alberto Quaini

#ifndef UTILS_H
#define UTILS_H

#include <RcppArmadillo.h>
#include <vector>

// Set template T, representing both arma::mat or arma::vec in subsequent functions.
template<typename T>

// Function for internal use.
// Solve a system of linear equations
// A * solution = b
// where A is likely symmetric and positive definite,
// and b is a vector or a matrix.
arma::mat solve_sympd(const arma::mat& A, const T& b) {

  try {
    // Try arma::solve
    return arma::solve(
      A, b, arma::solve_opts::likely_sympd + arma::solve_opts::no_approx
    );

  } catch (const std::runtime_error&) {
    // Fallback to generalized inverse
    return arma::solve(
      A, b, arma::solve_opts::force_approx
    );

  }

}

// Helper: Recursively generate all combinations of k indices from {0,1,...,n-1}
void generate_combinations(const unsigned int n, const unsigned int k, const unsigned int offset,
                           arma::uvec &current, std::vector<arma::uvec> &all);

// // Helper: Compute nCk (number of combinations of n elements in k places)
// unsigned long long nCk(const unsigned int n, const unsigned int k);

// // Helper: Generate a random combination of k distinct indices from 0 to n-1
// arma::uvec random_combination(const unsigned int n, const unsigned int k);

#endif // UTILS_H
