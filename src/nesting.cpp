#include "nesting.h"

arma::mat HDRNesting::sample_from_nesting(int n_samples, arma::vec x_init, int n_skip) {
  arma::mat samples = ess(n_samples, this->shifted_lincon, x_init, false);
  return samples.cols(arma::regspace<arma::uvec>(0, n_skip, samples.n_cols));
}

void HDRNesting::compute_log_nesting_factor(arma::mat X) {
  this->log_conditional_probability = 
    std::log(arma::sum(this->shifted_lincon.integration_domain(X))) - 
    std::log(X.n_cols);
}

void SubsetNesting::update_properties_from_samples(arma::mat X) {
  
  this->n_inside = X.n_cols * this->fraction;
  
  // Update log conditional probability
  this->compute_log_nesting_factor(X);
  // Rcpp::Rcout << "Computed log nesting factor" << std::endl;
  arma::mat eval = this->lincon.evaluate(X);
  // Rcpp::Rcout << "Evaluated X" << std::endl;
  arma::vec shiftvals = - arma::trans(arma::min(eval, 0));
  
  Rcpp::Rcout << "Shiftvals: " << shiftvals << std::endl;
  Rcpp::Rcout << "n_inside: " << this->n_inside << std::endl;
  
  // Rcpp::Rcout << "Got shiftvals" << std::endl;
  
  arma::uvec idx_inside;
  arma::uvec neg_shiftvals = arma::find(shiftvals < 0);
  if (neg_shiftvals.n_elem > this->n_inside) {
    // consider failure domain directly
    this->shift = 0.0;
    idx_inside = this->update_fix_shift(this->shift, shiftvals);
  }
  else {
    std::pair<double, arma::uvec> results = this->update_find_shift(shiftvals);
    this->shift = results.first;
    idx_inside = results.second;
  }
  
  // Rcpp::Rcout << "Did update" << std::endl;
  
  // always assume that we save one sample from inside domain
  // arma::vec probs = (1 / idx_inside.n_elem) * arma::ones(idx_inside.n_elem);
  // arma::uvec sequence = arma::linspace<arma::uvec>(0, idx_inside.n_elem-1, 
  //                                                  idx_inside.n_elem);
  // int idx = Rcpp::sample(idx_inside.n_elem, 1, false, probs, false)[0];
  // int idx = Rcpp::RcppArmadillo::sample(sequence, 1, false, probs)(0);
  
  Rcpp::IntegerVector choices = Rcpp::seq(0, idx_inside.n_elem-1);
  int idx = Rcpp::sample(choices, 1, false, R_NilValue)[0];
  
  // Rcpp::Rcout << "sampled index" << std::endl;
  
  this->x_in = X.col(idx);
  this->shifted_lincon = ShiftedLinearConstraints(this->lincon.A,
                                                  this->lincon.b,
                                                  this->shift,
                                                  true);
  return;
}

void SubsetNesting::compute_log_nesting_factor(arma::mat X) {
  this->log_conditional_probability = 
    std::log(X.n_cols * this->fraction) - std::log(X.n_cols);
}

arma::mat SubsetNesting::sample_from_nesting(int n, arma::vec x_init, int n_skip) {
  arma::mat samples =  ess(n, this->shifted_lincon, x_init, false);
  return samples.cols(arma::regspace<arma::uvec>(0, n_skip, samples.n_cols));
}

arma::uvec SubsetNesting::update_fix_shift(double shift, arma::vec shiftvals) {
  // arma::uvec less_than = shiftvals < shift;
  // arma::uvec nonzeros = arma::find(less_than > 0);
  arma::uvec nonzeros = arma::find(shiftvals < shift);
  this->n_inside = nonzeros.n_elem;
  
  return nonzeros;
}

std::pair<double, arma::uvec> SubsetNesting::update_find_shift(arma::vec shiftvals) {
  arma::uvec idx = arma::sort_index(shiftvals);
  arma::vec sorted_inside_shiftvals = shiftvals(idx.head(this->n_inside));
  // arma::vec sorted_shiftvals = shiftvals(idx);
  // arma::vec sorted_inside_shiftvals = sorted_shiftvals.head(this->n_inside);
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

void test_update_find_shift(arma::mat A, arma::vec b, double fraction,
                            double shift, arma::vec shiftvals) {
  LinearConstraints lincon = LinearConstraints(A, b, true);
  SubsetNesting sn(lincon, fraction);
  std::pair<double, arma::uvec> result = sn.update_find_shift(shiftvals);
}

