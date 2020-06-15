#include "angle_sampler.h"

//' Compute the cumulatve lengths of the slices, with a zero prepended 
//' :return: cumulative lengths of the slices
arma::vec AngleSampler::get_slices_cumulative_length() {
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
  
  double sample = l * arma::randu();
  auto lower = std::lower_bound(cum_len)
  int idx = std::distance(cum_len.begin(), lower);
  
  return 2.;
}