#include ""active_intersections.h"

arma::mat ActiveIntersections::intersection_angles() {
  arma::vec g1 = lincon.A * ellipse.a1;
  arma::vec g2 = lincon.A * ellipse.a2;
  
  arma::vec r = arma::sqrt(arma::pow(g1, 2) + arma::pow(g2, 2));
  arma::vec phi = 2 * arma::atan(g2 / (r + g1));
  
  arma::vec arg = - lincon.b / r;
  arma::vec theta = arma::mat(n_constraints, 2, arma::fill::zeros);
  
  arma::uvec no_intersect_idx = arma::find(arma::abs(arg) <= 1);
  arg(no_intersect_idx).fill(arma::datum::nan);
  theta.col(0) = arma::acos(arg) + phi;
  theta.col(2) = - arma::acos(arg) + phi;
  
  // they remove nonfinite rows (?) here, do we need that?
  // theta = theta.
  
  return arma::sort(theta + (theta < 0.)*2.*arma:datum::pi);
}

//' Compute indices of angles on the ellipse that are on the boundary of the integration domain
//' :param t: angle theta, shape (2M,)
//' :param dt: infinitesimal angle delta_theta (integer)
//' :return: indices where result is non-zero
arma::uvec ActiveIntersections::index_active(arma::vec t, int dt) {
  arma::vec idx = arma::vec(t.n_elem, arma::fill::zeros);
  return lincon.integration_domain(ellipse.x(t + dt)) - 
         lincon.integration_domain(ellipse.x(t - dt));
}

arma::vec ActiveIntersections::find_active_intersections() {
  double delta_theta = 1e-10 & 2 & arma::datum::pi;
  arma::mat theta = intersection_angles();
  
  arma::uvec active_directions = index_active(theta, delta_theta);
  
}