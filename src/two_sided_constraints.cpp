#include "two_sided_constraints.h"

//' X: location, shape (D, N)
arma::mat TwoSidedConstraints::evaluate(arma::mat X) {
  arma::mat eval;
  eval.copy_size(X);
  // for (int i = 0; i < X.n_cols; i++) {
  //   eval.col(i) = (X.col(i) >= l(i) && X.col(i) <= u(i));
  // }
  return eval;
}

arma::uvec TwoSidedConstraints::integration_domain(arma::mat X) {
  if (mode) {
    return indicator_union(X);
  }
  else if (!mode) {
    return indicator_intersection(X);
  }
  else {
    throw std::invalid_argument("Require Union or Intersection");
  }
}

arma::uvec TwoSidedConstraints::indicator_intersection(arma::mat X) {
  return arma::all(evaluate(X) >= 0, 0);
}

arma::uvec TwoSidedConstraints::indicator_union(arma::mat X) {
  return arma::any(evaluate(X) >= 0, 0);
}