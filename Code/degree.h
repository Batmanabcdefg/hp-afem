#pragma once

#include "math.h"

// TODO: convert this to a namespace rather than a class.
class Degree {
public:
  static const int Constant   = 1;
  static const int Linear     = 3;
  static const int Quadratic  = 6;
  static const int Cubic      = 10;
  static const int Quartic    = 15;
  static const int Quintic    = 21;
  static int dofToDegree(int dof) { return Math::degree(dof); }
  static int degreeToDim(int degree) { return Math::triNum(degree); }
  static int dofToDim(int dof) { return degreeToDim(dofToDegree(dof)); }
};
