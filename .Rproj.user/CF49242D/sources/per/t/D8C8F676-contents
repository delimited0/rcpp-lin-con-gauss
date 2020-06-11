#ifndef LINEARCONSTRAINTS_H
#define LINEARCONSTRAINTS_H

#include "util.h"

class LinearConstraints {
public
  arma::mat A;
  arma::vec b;
  int n_constraints;
  int n_dim;
  String mode;
  
  LinearConstraints(arma::mat A, arma::vec b, String mode) : 
    A(A), b(b), mode(mode), n_constraints(b.size()), n_dim(A.ncols) {
  }
  
  arma::vec evaluate(arma::mat X);
  bool integration_domain(arma::mat X);
  bool indicator_intersection(arma::mat X);
  bool indicator_union(arma::mat X);
};

#endif