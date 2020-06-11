#ifndef LOOP_H
#define LOOP_H

#include "util.h"

class Loop {
public
  LinearConstraints lincon;
  int n_skip;
  
  Loop(LinearConstraints lincon, int n_skip) : lincon(lincon), n_skip(n_skip) {
  }
  void run();
};

#endif LOOP_H