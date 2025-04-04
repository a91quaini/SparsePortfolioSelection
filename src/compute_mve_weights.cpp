#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Compute Mean-Variance Efficient (MVE) Portfolio Weights
//'
//' Given a risk-aversion parameter \eqn{\gamma}, a first moment vector, a second moment matrix,
//' and a selection index vector, this function computes the mean variance efficient portfolio weights
//' \deqn{w = \frac{1}{\gamma}\, \text{second\_moment}^{-1}\, \text{first\_moment}},
//' over the selected assets.
//' If the provided asset selection vector has length less than N,
//' the function returns an N-length weight vector with zeros for the unselected assets.
//'
//' @param first_moment First moment vector.
//' @param second_moment Second moment matrix.
//' @param selection Unsigned integer vector with asset indices.
//' @param gamma Risk aversion parameter. Default is 1.
//' @param do_checks Logical flag indicating whether to perform input checks (default = FALSE).
//' @return An N-length vector of mean variance efficient weights.
//' @export
// [[Rcpp::export]]
arma::vec compute_mve_weights(const arma::vec& first_moment,
                      const arma::mat& second_moment,
                      const arma::uvec& selection,
                      const double gamma = 1.0,
                      const bool do_checks = false) {

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
