#include "nesting.h"

arma::mat HDRNesting::sample_from_nesting(int n_samples, arma::vec x_init, int n_skip) {
  return ess(n_samples, this->shifted_lincon, x_init, false);
}

void HDRNesting::compute_log_nesting_factor(arma::mat X) {
  this->log_conditional_probability = 
    std::log(arma::sum(this->shifted_lincon.integration_domain(X)));
}


arma::uvec SubsetNesting::update_fix_shift(double shift, arma::mat shiftvals) {
  arma::umat less_than = shiftvals < shift;
  arma::uvec nonzeros = arma::find(less_than > 0);
  this->n_inside = nonzeros.n_elem;
  
  arma::ucolvec row_mask = arma::any(less_than, 1);
  return arma::find(row_mask > 0);
}

std::pair<double, arma::uvec> SubsetNesting::update_find_shift(arma::mat shiftvals) {
  // find the elements of shiftvals with order at least n_inside
}

