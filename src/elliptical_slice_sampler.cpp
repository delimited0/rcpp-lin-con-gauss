#include "angle_sampler.h"

// [[Rcpp::export]]
arma::mat ess(int n_iterations, arma::mat A, arma::vec b, arma::vec x_init, 
              bool mode = true) {
  LinearConstraints lincon = LinearConstraints(A, b, mode);
  arma::mat samples(lincon.n_dim, n_iterations);
  
  arma::vec x0;
  arma::vec x1 = x_init;
  for (int i = 0; i < n_iterations; i++) {
    x0 = x1;
    x1 = arma::randn(lincon.n_dim, 1);
    Ellipse ellipse = Ellipse(x0, x1);
    ActiveIntersections active_intersections = ActiveIntersections(ellipse, lincon);
    AngleSampler slice_sampler = AngleSampler(active_intersections);
    
    if (!active_intersections.ellipse_in_domain) {
      throw std::domain_error("At least one point should be in the domain");
    }
    
    double t_new = slice_sampler.draw_angle();
    x1 = ellipse.x(t_new);
    samples.col(i) = x1;
  }
  
  return samples;
}