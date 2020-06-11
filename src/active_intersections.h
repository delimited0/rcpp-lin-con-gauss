#ifndef ACTIVEINTERSECTIONS_H
#define ACTIVEINTERSECTIONS_H

#include "linear_constraints.h"
#include "ellipse.h"

class ActiveIntersections {
public
  Ellipse ellipse;
  LinearConstraints lincon;
  int n_constraints;
  bool ellipse_in_domain;
  
  ActiveIntersections(Ellipse ellipse, LinearConstraints lincon) :
    ellipse(ellipse), lincon(lincon), 
    n_constraints(lincon.b.n_elem), ellipse_in_domain(true) {
  }
  
  arma::mat intersection_angles();
  arma::vec find_active_intersections();
  std::pair<double, arma::vec> rotated_intersections();
  arma::uvec index_active(arma::vec t, int dt);
};

#endif ACTIVEINTERSECTIONS_H