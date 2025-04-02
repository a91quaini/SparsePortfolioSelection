#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int add_ints(int x, int y) {
  return x + y;
}
