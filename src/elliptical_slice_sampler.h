#ifndef ELLIPTICAL_SLICE_SAMPLER_H
#define ELLIPTICAL_SLICE_SAMPLER_H

#include "angle_sampler.h"

arma::mat ess(int n_iterations, arma::mat A, arma::vec b, arma::vec x_init, 
              bool mode, bool verbose);

arma::mat ess(int n_iterations, LinearConstraints lincon, arma::vec x_init, 
              bool verbose);

#endif