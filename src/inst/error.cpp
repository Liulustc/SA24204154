#include <Rcpp.h>
using namespace Rcpp;
//' @title Compute Community Assignment Error
//' @description Calculate the error between the true community vector and the estimated community vector
//' @param n the length of the input community vector
//' @param c the true community assignment vector
//' @param es_c the estimated community assignment vector
//' @return the computational error between two input vectors
//' @export
// [[Rcpp::export]]
double error(int n, IntegerVector c, IntegerVector es_c) {
  int error1 = 0;
  int error2 = 0;
  for (int j = 0; j < n; j++) {
    for (int i = 0; i < j; i++) {
      if ((c[i] == c[j]) && (es_c[i] != es_c[j])) {
        error1++;
      }
      if ((c[i] != c[j]) && (es_c[i] == es_c[j])) {
        error2++;
      }
    }
  }
  return 2.0 * (error1 + error2) / (n * (n - 1));
}
