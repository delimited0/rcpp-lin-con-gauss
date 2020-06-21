#include "ellipse.h"

arma::mat Ellipse::x(const arma::vec & theta) {
  arma::mat result(this->a1.n_elem, theta.n_elem);
  for (int i = 0; i < theta.n_elem; i++) {
    result.col(i) = this->a1 * std::cos(theta(i)) + this->a2 * std::sin(theta(i));
  }
  return result;
}

arma::vec Ellipse::x(double theta) {
  return this->a1 * std::cos(theta) + this->a2 * std::sin(theta);
}