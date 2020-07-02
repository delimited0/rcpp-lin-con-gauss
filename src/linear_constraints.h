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
  
  arma::mat evaluate(const arma::mat & X);
  arma::vec integration_domain(const arma::mat & X);
  arma::vec indicator_intersection(const arma::mat & X);
  arma::vec indicator_union(const arma::mat & X);
};

class ShiftedLinearConstraints : public LinearConstraints {
public:
  double shift;
  
  ShiftedLinearConstraints(arma::mat A, arma::vec b, double shift, bool mode) :
    LinearConstraints(A, b + shift, mode), shift(shift) {
  }
};

#endif