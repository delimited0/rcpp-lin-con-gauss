#include "ellipse.h"

arma::vec x(arma::vec theta) {
  return a1 * std::cos(theta) + a2 * std::sin(theta);
}