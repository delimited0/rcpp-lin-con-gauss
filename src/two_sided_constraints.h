#ifndef TWOSIDEDCONSTRAINTS_H
#define TWOSIDEDCONSTRAINTS_H

#include "util.h"

class TwoSidedConstraints {
public:
  arma::vec l;
  arma::vec u;
  int n_constraints;
  int n_dim;
  bool mode;
  
  TwoSidedConstraints(arma::vec l, arma::vec u) : 
    l(l), u(u), mode(mode), n_constraints(l.n_elem), n_dim(l.n_elem) {
  }
  
  arma::mat evaluate(arma::mat X);
  arma::uvec integration_domain(arma::mat X);
  arma::uvec indicator_intersection(arma::mat X);
  arma::uvec indicator_union(arma::mat X);
};

#endif TWOSIDEDCONSTRAINTS_H