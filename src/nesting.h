#ifndef NESTING_H
#define NESTING_H

#include "linear_constraints.h"
#include "elliptical_slice_sampler.h"

class HDRNesting {
public:
  ShiftedLinearConstraints shifted_lincon;
  double shift;
  int dim;
  double log_conditional_probability;
  
  HDRNesting(LinearConstraints lincon, double shift) :
    shift(shift), shifted_lincon(lincon.A, lincon.b, shift, true) {
    // shifted_lincon = ShiftedLinearConstraints(lincon.A, lincon.b, shift, true);
    dim = shifted_lincon.n_dim;
    double log_conditional_probability = 0.0;
  }
  
  arma::mat sample_from_nesting(int n_samples, arma::vec x_init, int n_skip);
  void compute_log_nesting_factor(arma::mat x);
};

class SubsetNesting {
public: 
  double fraction;
  LinearConstraints lincon;
  int n_save;
  int n_inside;
  double log_conditional_probability;
  double shift;
  arma::vec x_in;
  ShiftedLinearConstraints shifted_lincon;
  
  SubsetNesting(LinearConstraints lincon, double fraction, int n_save = 1) :
    lincon(lincon), fraction(fraction), n_save(n_save), 
    shifted_lincon(lincon.A, lincon.b, 0.0, true) {
  }
  
  void update_properties_from_samples(arma::mat X);
  arma::mat sample_from_nesting(int n, arma::vec x_init, int n_skip);
  void compute_log_nesting_factor(arma::mat X);
  std::pair<double, arma::uvec> update_find_shift(arma::mat shiftvals);
  arma::uvec update_fix_shift(double shift, arma::mat shiftvals);
};
  

#endif