#ifndef UTILS_H
#define UTILS_H

#include <RcppArmadillo.h>
#include <vector>

// Helper: Recursively generate all combinations of k indices from {0,1,...,n-1}
void generate_combinations(const unsigned int n, const unsigned int k, const unsigned int offset,
                           arma::uvec &current, std::vector<arma::uvec> &all);

// // Helper: Compute nCk (number of combinations of n elements in k places)
// unsigned long long nCk(const unsigned int n, const unsigned int k);

// // Helper: Generate a random combination of k distinct indices from 0 to n-1
// arma::uvec random_combination(const unsigned int n, const unsigned int k);

#endif // UTILS_H
