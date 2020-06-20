#include "active_intersections.h"

arma::vec ActiveIntersections::intersection_angles() {
  arma::vec g1 = lincon.A * ellipse.a1;
  arma::vec g2 = lincon.A * ellipse.a2;
  
  // Rcpp::Rcout << "g1, g2: (" << g1 << ", " << g2 << ")" << std::endl;
  
  arma::vec r = arma::sqrt(arma::pow(g1, 2) + arma::pow(g2, 2));
  arma::vec phi = 2 * arma::atan(g2 / (r + g1));
  
  // Rcpp::Rcout << "r: " << r << std::endl;
  // Rcpp::Rcout << "phi: " << phi << std::endl;
  
  arma::vec arg = - lincon.b / r;
  
  // Rcpp::Rcout << "arg: " << arg << std::endl;

  arma::mat theta = arma::mat(n_constraints, 2, arma::fill::zeros);
  
  arma::uvec no_intersect_idx = arma::find(arma::abs(arg) > 1);
  arg(no_intersect_idx).fill(arma::datum::nan);
  theta.col(0) = arma::acos(arg) + phi;
  theta.col(1) = - arma::acos(arg) + phi;
  
  // Rcpp::Rcout << "Theta: " << theta << std::endl;
  
  arma::vec thetaf = theta.elem(arma::find_finite(theta));
  
  // Rcpp::Rcout << "thetaf: " << thetaf << std::endl;
  
  arma::uvec neg_idx = arma::find(thetaf < 0);
  thetaf(neg_idx) =  thetaf(neg_idx) + 2.*arma::datum::pi;
  
  return arma::sort(thetaf);
}

//' Compute indices of angles on the ellipse that are on the boundary of the integration domain
//' :param t: angle theta, shape (2M,)
//' :param dt: infinitesimal angle delta_theta (integer)
//' :return: indices where result is non-zero
arma::vec ActiveIntersections::index_active(arma::vec t, double dt) {
  arma::mat ellplus = ellipse.x(t + dt);
  arma::mat ellminus = ellipse.x(t - dt);
  arma::vec ellplus_int_domain = lincon.integration_domain(ellplus);
  arma::vec ellminus_int_domain = lincon.integration_domain(ellminus);
  return ellplus_int_domain - ellminus_int_domain;
}

arma::vec ActiveIntersections::find_active_intersections() {
  double delta_theta = 1e-10 * 2. * arma::datum::pi;
  arma::mat theta = this->intersection_angles();
  
  // Rcpp::Rcout << "theta: " << theta << std::endl;
  
  arma::vec active_directions = this->index_active(theta, delta_theta);
  
  // Rcpp::Rcout << "active directions: " << active_directions << std::endl;
  
  arma::uvec active_nonzero = arma::find(active_directions != 0.);
  arma::vec theta_active = theta(active_nonzero);
  
  // Rcpp::Rcout << "theta active: " << theta_active << std::endl;
  
  while (theta_active.n_elem % 2 == 1) {
    // Almost tangential ellipses, reduce delta_theta
    delta_theta = 1.e-1 * delta_theta;
    active_directions = this->index_active(theta, delta_theta);
    active_nonzero = arma::find(active_directions != 0);
    theta_active = theta(active_nonzero);
  }
  
  // Rcpp::Rcout << "reduced delta theta active: " << theta_active << std::endl;
  
  if (theta_active.n_elem == 0) {
    theta_active = {0.0, 2 * arma::datum::pi};
    double u = arma::randu();
    arma::vec x = this->ellipse.x(2 * arma::datum::pi * u);
    if (!this->lincon.integration_domain(x)(0))
      // entire ellipse is outside of domain
      this->ellipse_in_domain = false;
  }
  else {
    active_nonzero = arma::find(active_directions != 0);
    arma::vec ad_nonzero = active_directions(active_nonzero);
    if ((ad_nonzero(0) + 1) < .0001) {
      double theta_active_first = theta_active(0);
      theta_active.head(theta_active.n_elem-1) = theta_active.tail(theta_active.n_elem-1);
      theta_active.tail(1) = theta_active_first;
    }
  }
  
  return theta_active;
}

std::pair<double, arma::vec> ActiveIntersections::rotated_intersections() {
  arma::vec slices = this->find_active_intersections();
  double rotation_angle = slices[0];
  slices = slices - rotation_angle;
  
  std::pair<double, arma::vec> result;
  result.first = rotation_angle;
  
  arma::uvec neg_idx = arma::find(slices < 0);
  slices(neg_idx) = slices(neg_idx) + 2 * arma::datum::pi;
  result.second = slices;
  
  return result;
}