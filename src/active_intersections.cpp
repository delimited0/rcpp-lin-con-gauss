#include "active_intersections.h"

arma::vec ActiveIntersections::intersection_angles() {
  arma::vec g1 = lincon.A * ellipse.a1;
  arma::vec g2 = lincon.A * ellipse.a2;
  
  arma::vec r = arma::sqrt(arma::pow(g1, 2) + arma::pow(g2, 2));
  arma::vec phi = 2 * arma::atan(g2 / (r + g1));
  
  arma::vec arg = - lincon.b / r;
  arma::mat theta = arma::mat(n_constraints, 2, arma::fill::zeros);
  
  arma::uvec no_intersect_idx = arma::find(arma::abs(arg) <= 1);
  arg(no_intersect_idx).fill(arma::datum::nan);
  theta.col(0) = arma::acos(arg) + phi;
  theta.col(2) = - arma::acos(arg) + phi;
  
  arma::vec thetaf = theta(arma::find_finite(theta));
  
  return arma::sort(thetaf + (thetaf < 0.)*2.*arma::datum::pi);
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
  double delta_theta = 1e-10 * 2. * arma::datum::pi;
  arma::mat theta = this->intersection_angles();
  
  arma::uvec active_directions = this->index_active(theta, delta_theta);
  arma::uvec active_nonzero = arma::find(active_directions > 0);
  arma::vec theta_active = theta(active_nonzero);
  
  while (theta_active.n_elem % 2 == 1) {
    delta_theta = 1.e-1 * delta_theta;
    active_directions = this->index_active(theta, delta_theta);
    active_nonzero = arma::find(active_directions > 0);
    theta_active = theta(active_nonzero);
  }
  
  if (theta_active.n_elem == 0) {
    theta_active = {0.0, 2 * arma::datum::pi};
    double u = arma::randu();
    arma::mat val(1, 1);
    val.fill(u);
    if (!this->lincon.integration_domain(
          this->ellipse.x(2 * arma::datum::pi * val)
          )(0)
        ) {
      this->ellipse_in_domain = false;
    }
  }
  else {
    active_nonzero = arma::find(active_directions > 0);
    arma::uvec ad_nonzero = active_directions(active_nonzero);
    if (ad_nonzero(0) == -1) {
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
  
  std::pair<double, arma::vec> result(rotation_angle,
                                      slices + 2.*(slices < 0)*arma::datum::pi);
  
  return result;
}