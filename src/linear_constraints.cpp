#include "linear_constraints.h"

//'
//' X: location, shape (D, N)
arma::mat LinearConstraints::evaluate(arma::mat X) {
  arma::mat eval(this->n_constraints, X.n_cols);
  for (int i = 0; i < eval.n_cols; i++) {
    eval.col(i) = this->A * X.col(i) + this->b;
  }
  return eval;
}

//' is 1 if x is in the integration domain, else 0
//' X: location, shape (D, N)
arma::vec LinearConstraints::integration_domain(arma::mat X) {
  if (mode) {
    return indicator_intersection(X);
  }
  else if (!mode) {
    return indicator_union(X);
  }
  else {
    throw std::invalid_argument("Require Union or Intersection");
  }
}

arma::vec LinearConstraints::indicator_intersection(arma::mat X) {
  arma::umat eval = evaluate(X) >= 0;
  arma::vec result(X.n_cols);
  for (int i = 0; i < result.n_elem; i++) {
    result(i) = arma::prod(eval.col(i));
  }
  return result;
}

arma::vec LinearConstraints::indicator_union(arma::mat X) {
  arma::umat eval = evaluate(X) < 0;
  arma::vec result(X.n_cols);
  for (int i = 0; i < result.n_elem; i++) {
    result(i) = arma::prod(eval.col(i));
  }
  return 1 - result;
}

