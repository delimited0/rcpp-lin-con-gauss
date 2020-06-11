#include "linear_constraints.h"

//'
//' X: location, shape (D, N)
LinearConstraints::evaluate(arma::mat X) {
  arma::vec eval(X.ncol);
  for (int i = 0; i < X.ncol; i++) {
    eval(i) = A * X.col(i) + b;
  }
  return eval;
}

//'
//' is 1 if x is in the integration domain, else 0
//' X: location, shape (D, N)
LinearConstraints::integration_domain(arma::mat X) {
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

LinearConstraints::indicator_intersection(arma::mat X) {
  return arma::all(evaluate(X) >= 0);
}

LinearConstraints::indicator_union(arma::mat X) {
  return arma::any(evaluate(X) >= 0);
}

