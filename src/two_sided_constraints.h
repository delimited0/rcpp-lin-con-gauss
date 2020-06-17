#ifndef TWOSIDEDCONSTRAINTS_H
#define TWOSIDEDCONSTRAINTS_H

#include "util.h"

class TwoSidedConstraints {
public:
  arma::vec l;
  arma::vec u;
  int n_constraints;
  int n_dim;
  
  TwoSidedConstraints(arma::vec l, arma::vec u) : 
    l(l), u(u), n_constraints(l.n_elem, n_dim(l.n_elem){
  }
  
  arma::uvec evaluate(arma::mat X);
  bool integration_domain(arma::mat X);
  bool indicator_intersection(arma::mat X);
  bool indicator_union(arma::mat X);
};