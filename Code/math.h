#pragma once
#include <limits>
#include <vector>
#include <cmath>

#include "config.h"

// TODO: convert this to a namespace instead of a class.
class Math {
public:
  static bool almost_equal(scalar x, scalar y, int ulp);
  static scalar recip(scalar a, scalar b);

  static int triNum( int n);
  static int degree( int dof);

  static int pow2roundup(int n);
};
