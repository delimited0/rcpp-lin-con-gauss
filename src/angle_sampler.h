#ifndef ANGLESAMPLER_H
#define ANGLESAMPLER_H

#include "active_intersections.h"

class AngleSampler {
public
  ActiveIntersections active_intersections;
  double rotation_angle;
  arma::vec rotated_slices;
  
  AngleSampler(ActiveIntersections active_intersections) :
    active_intersections(active_intersections) {
    std::pair<double, arma::vec> 
  }
};



#endif ANGLESAMPLER_H

