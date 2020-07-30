#include "../src/angle_sampler.h"

// [[Rcpp::export]]
arma::mat ellipse_x(arma::vec x0, arma::vec x1, arma::vec theta) {
  Ellipse ell = Ellipse(x0, x1);
  return ell.x(theta);
}

arma::mat ellpse_dx(arma::vec x0, arma::vec x1, double theta) {
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
                               arma::vec t, double dt) {
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

void rotated_intersections(arma::mat A, arma::vec b,
                           arma::vec x0, arma::vec x1) {
  Ellipse ell = Ellipse(x0, x1);
  LinearConstraints lincon = LinearConstraints(A, b, true);
  ActiveIntersections actint = ActiveIntersections(ell, lincon);
  std::pair<double, arma::vec> result = actint.rotated_intersections();
  Rcpp::Rcout << "Rotation angle: " << result.first;
}

// [[Rcpp::export]]
arma::vec integration_domain(arma::mat A, arma::vec b, arma::mat X,
                              bool mode) {
  LinearConstraints lincon = LinearConstraints(A, b, mode);
  return lincon.integration_domain(X);
}

// [[Rcpp::export]]
arma::mat evaluate(arma::mat A, arma::vec b, arma::mat X) {
  LinearConstraints lincon = LinearConstraints(A, b, true);
  return lincon.evaluate(X);
}

// [[Rcpp::export]]
void angle_sampler(arma::mat A, arma::vec b,
                  arma::vec x0, arma::vec x1) {
  Ellipse ell = Ellipse(x0, x1);
  LinearConstraints lincon = LinearConstraints(A, b, true);
  ActiveIntersections actint = ActiveIntersections(ell, lincon);
  AngleSampler as = AngleSampler(actint);
  Rcpp::Rcout << "Rotation angle: " << as.rotation_angle << std::endl;
  Rcpp::Rcout << "Rotated slices: " << as.rotated_slices << std::endl;
}

// [[Rcpp::export]]
double draw_angle(arma::mat A, arma::vec b,
                     arma::vec x0, arma::vec x1) {
  Ellipse ell = Ellipse(x0, x1);
  LinearConstraints lincon = LinearConstraints(A, b, true);
  ActiveIntersections actint = ActiveIntersections(ell, lincon);
  AngleSampler as = AngleSampler(actint);
  return as.draw_angle();
}

// [[Rcpp::export]]
arma::vec get_slices_cumulative_length(arma::mat A, arma::vec b,
                                       arma::vec x0, arma::vec x1) {
  Ellipse ell = Ellipse(x0, x1);
  LinearConstraints lincon = LinearConstraints(A, b, true);
  ActiveIntersections actint = ActiveIntersections(ell, lincon);
  AngleSampler as = AngleSampler(actint);
  return as.get_slices_cumulative_length();
}