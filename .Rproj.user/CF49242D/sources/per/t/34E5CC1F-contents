#include ""two_sided_constraints.h"

//' X: location, shape (D, N)
TwoSidedConstraints::evaluate(arma::mat X) {
  arma::vec eval(X.ncol);
  for (int i = 0; i < X.ncol; i++) {
    eval(i) = X.col(i) >= l(i) && X.col(i) <= u(i);
  }
  return eval;
}

TwoSidedConstraints::integration_domain(arma::mat X) {
  if (mode == "Union") {
    return indicator_union(X);
  }
  else if (mode == "Intersection") {
    return indicator_intersection(X);
  }
  else {
    throw std::invalid_argument("Require Union or Intersection");
  }
}

TwoSidedConstraints::indicator_intersection(arma::mat X) {
  return arma::all(evaluate(X) >= 0);
}

TwoSidedConstraints::indicator_union(arma::mat X) {
  return arma::any(evaluate(X) >= 0);
}