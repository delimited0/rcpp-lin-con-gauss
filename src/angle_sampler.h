#ifndef ANGLESAMPLER_H
#define ANGLESAMPLER_H

#include "active_intersections.h"

class AngleSampler {
public:
  ActiveIntersections active_intersections;
  double rotation_angle;
  arma::mat rotated_slices;
  
  AngleSampler(ActiveIntersections active_intersections) :
    active_intersections(active_intersections) {
    std::pair<double, arma::vec> result = active_intersections.rotated_intersections();
    rotation_angle = result.first;
    rotated_slices = arma::reshape(result.second, result.second.size() / 2, 2);
  }
  
  arma::vec get_slices_cumulative_length();
  double draw_angle();
};

#endif ANGLESAMPLER_H

