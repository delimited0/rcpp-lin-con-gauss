#include "../src/nesting.cpp"

// [[Rcpp::export]]
arma::uvec test_update_fix_shift(arma::mat A, arma::vec b, double fraction,
                            double shift, arma::vec shiftvals) {
  LinearConstraints lincon = LinearConstraints(A, b, true);
  SubsetNesting sn(lincon, fraction);
  return sn.update_fix_shift(shift, shiftvals);
}