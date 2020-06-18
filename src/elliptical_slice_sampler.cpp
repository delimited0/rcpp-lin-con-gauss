#include "angle_sampler.h"

// [[Rcpp::export]]
arma::mat ess(int n_iterations, arma::mat A, arma::vec b, arma::vec x_init, 
              bool mode = true, bool verbose = false) {
  Rcpp::Rcout << "Getting started" << std::endl;
  LinearConstraints lincon = LinearConstraints(A, b, mode);
  arma::mat samples(lincon.n_dim, n_iterations);
  
  Rcpp::Rcout << "Set up constraints" << std::endl;
  
  arma::vec x0;
  arma::vec x1 = x_init;
  for (int i = 0; i < n_iterations; i++) {
    if (verbose) Rcpp::Rcout << "Iteration: " << i << std::endl;
    
    x0 = x1;
    x1 = arma::randn(lincon.n_dim, 1);
    
    Ellipse ellipse = Ellipse(x0, x1);
    Rcpp::Rcout << "Made ellipse" << std::endl;
    Rcpp::Rcout << "a1: " << ellipse.a1 << std::endl;
    Rcpp::Rcout << "a2: " << ellipse.a2 << std::endl;
    ActiveIntersections active_intersections = ActiveIntersections(ellipse, lincon);
    Rcpp::Rcout << "Found active intersections" << std::endl;
    AngleSampler slice_sampler = AngleSampler(active_intersections);
    Rcpp::Rcout << "Made slice sampler" << std::endl;
    
    if (!active_intersections.ellipse_in_domain) {
      throw std::domain_error("At least one point should be in the domain");
    }
    
    double t_new = slice_sampler.draw_angle();
    Rcpp::Rcout << "Drew angle" << std::endl;
    x1 = ellipse.x(t_new);
    Rcpp::Rcout << "Storing sample" << std::endl;
    samples.col(i) = x1;
  }
  
  return samples;
}