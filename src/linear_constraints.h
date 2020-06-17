#ifndef LINEARCONSTRAINTS_H
#define LINEARCONSTRAINTS_H

#include "util.h"

class LinearConstraints {
public:
  arma::mat A;
  arma::vec b;
  int n_constraints;
  int n_dim;
  bool mode;
  
  LinearConstraints(arma::mat A, arma::vec b, bool mode) : 
    A(A), b(b), mode(mode), n_constraints(b.size()), n_dim(A.n_cols) {
  }
  
  arma::mat evaluate(arma::mat X);
  arma::uvec integration_domain(arma::mat X);
  arma::uvec indicator_intersection(arma::mat X);
  arma::uvec indicator_union(arma::mat X);
};

#endif