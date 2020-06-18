#include "linear_constraints.h"

//'
//' X: location, shape (D, N)
arma::mat LinearConstraints::evaluate(arma::mat X) {
  arma::mat eval(this->n_dim, X.n_cols);
  // eval.copy_size(X);
  Rcpp::Rcout << "eval dim: " << arma::size(eval) << std::endl;
  Rcpp::Rcout << "X dim: " << arma::size(X) << std::endl;
  Rcpp::Rcout << "A dim: " << arma::size(A) << std::endl;
  Rcpp::Rcout << "b dim: " << arma::size(b) << std::endl;
  for (int i = 0; i < eval.n_cols; i++) {
    Rcpp::Rcout << "ax+b dim: " << arma::size(this->A * X.col(i) + this->b) << std::endl;
    eval.col(i) = this->A * X.col(i) + this->b;
  }
  Rcpp::Rcout << "finished looping" << std::endl;
  return eval;
}

//' is 1 if x is in the integration domain, else 0
//' X: location, shape (D, N)
arma::vec LinearConstraints::integration_domain(arma::mat X) {
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

arma::vec LinearConstraints::indicator_intersection(arma::mat X) {
  return arma::prod(arma::conv_to<arma::mat>::from(evaluate(X) >= 0), 0);
}

arma::vec LinearConstraints::indicator_union(arma::mat X) {
  return 1 - arma::prod(arma::conv_to<arma::mat>::from(evaluate(X) < 0), 0);
}

