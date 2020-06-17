#include "angle_sampler.h"
#include <algorithm>

//' Compute the cumulatve lengths of the slices, with a zero prepended 
//' :return: cumulative lengths of the slices
arma::vec AngleSampler::get_slices_cumulative_length() {
  // Rcpp::Rcout << this->rotated_slices << std::endl;
  arma::vec lengths = this->rotated_slices.col(1) - this->rotated_slices.col(0);
  arma::vec cum_len = arma::cumsum(lengths);
  cum_len.insert_rows(0, 1);
  return cum_len;
}

//' Draw one sample angle from the given slice
//' :return: random angle from slice(s)
double AngleSampler::draw_angle() {
  arma::vec cum_len = this->get_slices_cumulative_length();
  double l = cum_len(cum_len.n_elem - 1);
  
  // std::vector<double> std_cum_len = 
  //   arma::conv_to<std::vector<double>>::from(cum_len);
  
  arma::vec::iterator it = cum_len.begin();
  
  double sample = l * arma::randu();
  auto lower = std::lower_bound(cum_len.begin(), cum_len.end(), sample);
  int idx = std::distance(cum_len.begin(), lower) - 1;
  
  return this->rotated_slices(idx, 0) + sample - 
    cum_len(idx) + this->rotation_angle;
}