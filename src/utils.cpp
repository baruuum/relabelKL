#include "utils.h"
using namespace Rcpp;

//' Log of the sum of exponential
//'
//' @param x a vector
//' @return log of the sum of exponentiated elements in \code{x}
//[[Rcpp::export]]
double lse(const NumericVector & x) {

    return log_sum_exp(x);

}
