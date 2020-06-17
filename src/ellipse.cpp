#include "ellipse.h"

arma::mat Ellipse::x(arma::mat theta) {
  return this->a1 * arma::cos(theta) + this->a2 * arma::sin(theta);
}

arma::vec Ellipse::x(double theta) {
  return this->a1 * arma::cos(theta) + this->a2 * arma::sin(theta);
}