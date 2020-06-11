#ifndef ELLIPSE_H
#define ELLIPSE_H

#include "util.h"

class Ellipse {
  arma::vec a1;
  arma::vec a2;
  
  Ellipse(arma::vec a1, arma::vec a2) : a1(a1), a2(a2) {
  }
  
  arma::vec x(double theta);
};

#endif ELLIPSE_H