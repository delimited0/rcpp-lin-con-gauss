#include "angle_sampler.h"

// [[Rcpp::export]]
arma::mat ellipse_x(arma::vec x0, arma::vec x1, arma::vec theta) {
  Ellipse ell = Ellipse(x0, x1);
  return ell.x(theta);
}

// [[Rcpp::export]]
arma::vec intersection_angles(arma::mat A, arma::vec b,
                              arma::vec x0, arma::vec x1) {
  Ellipse ell = Ellipse(x0, x1);
  LinearConstraints lincon = LinearConstraints(A, b, true);
  ActiveIntersections actint = ActiveIntersections(ell, lincon);
  return actint.intersection_angles();
}

// [[Rcpp::export]]
arma::vec index_active(arma::mat A, arma::vec b,
                               arma::vec x0, arma::vec x1, 
                               arma::vec t, int dt) {
  Ellipse ell = Ellipse(x0, x1);
  LinearConstraints lincon = LinearConstraints(A, b, true);
  ActiveIntersections actint = ActiveIntersections(ell, lincon);
  return actint.index_active(t, dt);
}

// [[Rcpp::export]]
arma::vec find_active_intersections(arma::mat A, arma::vec b,
                                    arma::vec x0, arma::vec x1) {
  Ellipse ell = Ellipse(x0, x1);
  LinearConstraints lincon = LinearConstraints(A, b, true);
  ActiveIntersections actint = ActiveIntersections(ell, lincon);
  return actint.find_active_intersections();
}

// [[Rcpp::export]]
arma::vec integration_domain(arma::mat A, arma::vec b, arma::mat X,
                              bool mode) {
  LinearConstraints lincon = LinearConstraints(A, b, mode);
  return lincon.integration_domain(X);
}