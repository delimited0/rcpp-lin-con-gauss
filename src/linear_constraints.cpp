#include "linear_constraints.h"

//'
//' X: location, shape (D, N)
arma::mat LinearConstraints::evaluate(arma::mat X) {
  arma::mat eval(this->n_dim, X.n_rows);
  // eval.copy_size(X);
  for (int i = 0; i < X.n_cols; i++) {
    eval.col(i) = A * X.col(i) + b;
  }
  return eval;
}

//'
//' is 1 if x is in the integration domain, else 0
//' X: location, shape (D, N)
arma::uvec LinearConstraints::integration_domain(arma::mat X) {
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

arma::uvec LinearConstraints::indicator_intersection(arma::mat X) {
  return arma::all(evaluate(X) >= 0, 0);
}

arma::uvec LinearConstraints::indicator_union(arma::mat X) {
  return arma::any(evaluate(X) >= 0, 0);
}

