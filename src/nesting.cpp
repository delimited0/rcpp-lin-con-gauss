#include "nesting.h"

arma::mat HDRNesting::sample_from_nesting(int n_samples, arma::vec x_init, int n_skip) {
  return ess(n_samples, this->shifted_lincon, x_init, false);
}

void HDRNesting::compute_log_nesting_factor(arma::mat X) {
  this->log_conditional_probability = 
    std::log(arma::sum(this->shifted_lincon.integration_domain(X)));
}

void SubsetNesting::update_properties_from_samples(arma::mat X) {
  this->n_inside = X.n_cols * this->fraction;
  
  // Update log conditional probability
  this->compute_log_nesting_factor(X);
  arma::vec shiftvals = arma::min(this->lincon.evaluate(X), 0);
  
  arma::uvec idx_inside;
  if (arma::sum(shiftvals < 0) > this->n_inside) {
    // consider failure domain directly
    this->shift = 0.0;
    idx_inside = this->update_fix_shift(this->shift, shiftvals);
  }
  else {
    std::pair<double, arma::uvec> results = this->update_find_shift(shiftvals);
    this->shift = results.first;
    idx_inside = results.second;
  }
}

void SubsetNesting::compute_log_nesting_factor(arma::mat X) {
  this->log_conditional_probability = 
    std::log(X.n_cols * this->fraction) - std::log(X.n_cols);
}

arma::mat SubsetNesting::sample_from_nesting(int n, arma::vec x_init, int n_skip) {
  return ess(n, this->shifted_lincon, x_init, false);
}

arma::uvec SubsetNesting::update_fix_shift(double shift, arma::vec shiftvals) {
  arma::uvec less_than = shiftvals < shift;
  arma::uvec nonzeros = arma::find(less_than > 0);
  this->n_inside = nonzeros.n_elem;
  
  return arma::find(less_than > 0);
}

std::pair<double, arma::uvec> SubsetNesting::update_find_shift(arma::vec shiftvals) {
  arma::uvec idx = arma::sort_index(shiftvals);
  arma::vec sorted_shiftvals = shiftvals(idx);
  arma::vec sorted_inside_shiftvals = sorted_shiftvals.head(this->n_inside);
  double shift = sorted_inside_shiftvals(sorted_inside_shiftvals.n_elem - 1);
  std::pair <double, arma::uvec> result;
  result.first = shift;
  result.second = idx.head(idx.n_elem - 1);
  
  return result;
}

// [[Rcpp::export]]
arma::uvec test_update_fix_shift(arma::mat A, arma::vec b, double fraction,
                                 double shift, arma::vec shiftvals) {
  LinearConstraints lincon = LinearConstraints(A, b, true);
  SubsetNesting sn(lincon, fraction);
  return sn.update_fix_shift(shift, shiftvals);
}
