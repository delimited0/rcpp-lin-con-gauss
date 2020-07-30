#include "nesting.h"
#include "elliptical_slice_sampler.h"

//' Subset simulation to find a linearly constrained probability of failure in
//' a Gaussian space
//' @param A 
//' @param n_samples 
//' @param domain_fraction fraction of samples that should lie in the new domain
//' @param n_skip
// [[Rcpp::export]]
arma::vec hdr_prob(arma::mat A, arma::vec b, bool mode,
                 int n_sub_samples, double domain_fraction, int n_sub_skip,
                 int n_hdr_samples, int n_hdr_skip) {
  
  LinearConstraints lincon = LinearConstraints(A, b, mode);
  
  // subset nesting
  
  arma::mat X = arma::randn(lincon.n_dim, n_sub_samples);
  SubsetNesting sn = SubsetNesting(lincon, domain_fraction);
  sn.update_properties_from_samples(X);
  
  // outputs of subset nesting that we'll populate for HDR
  arma::mat X_inits(lincon.n_dim, 1);
  arma::vec shift_sequence(1);
  shift_sequence.fill(sn.shift);
  X_inits.col(0) = sn.x_in;
  
  while (sn.shift > 1e-6) {
    Rcpp::Rcout << "Subset sim looping" << std::endl;
    
    // sample from domain using elliptical slice sampler
    X = sn.sample_from_nesting(n_sub_samples, sn.x_in, n_sub_skip);
    
    // create new nesting and save the shift
    sn = SubsetNesting(lincon, domain_fraction);
    sn.update_properties_from_samples(X);
    arma::vec shift = {sn.shift};
    shift_sequence.insert_rows(shift_sequence.n_elem, shift);
    X_inits.insert_cols(X_inits.n_cols, sn.x_in);
  }
  
  Rcpp::Rcout << "Shift sequence: " << shift_sequence << std::endl;
  
  // HDR
  
  X = arma::randn(lincon.n_dim, n_sub_samples);
  
  HDRNesting current_nesting = HDRNesting(lincon, shift_sequence[0]);
  current_nesting.compute_log_nesting_factor(X);
  arma::vec logprobs(shift_sequence.n_elem);
  logprobs[0] = current_nesting.log_conditional_probability;
  
  for (int i = 1; i < shift_sequence.n_elem; i++) {
    
    X = current_nesting.sample_from_nesting(n_hdr_samples, X_inits.col(i), n_hdr_skip);
    
    current_nesting = HDRNesting(lincon, shift_sequence[i]);
    current_nesting.compute_log_nesting_factor(X);
    
    logprobs[i] = current_nesting.log_conditional_probability;
  }
  
  return logprobs;
}



