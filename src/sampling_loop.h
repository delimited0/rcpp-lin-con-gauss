#ifndef SAMPLINGLOOP_H
#define SAMPLINGLOOP_H

#include "loop.h"

class SamplingLoop: public Loop {
public
  int n_interations;
  
  SamplingLoop(int n_iterations, LinearConstraints linear_constraints, int n_skip) :
    Loop(linear_constraints, n_skip), n_iterations(n_iterations) {
  }
  
  virtual arma::vec compute_next_point(arma::vec x0);
  virtual bool is_converged();
};

#endif SAMPLINGLOOP_H