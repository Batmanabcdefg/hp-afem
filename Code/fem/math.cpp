#include <iostream>

#include "math.h"

using namespace std;

bool Math::almost_equal(scalar x, scalar y, int ulp) {
  return std::abs(x-y) < std::numeric_limits<scalar>::epsilon() * std::abs(x+y) * ulp
      || std::abs(x-y) < std::numeric_limits<scalar>::min();
}

scalar Math::recip(scalar a, scalar b) {
  if( a == 0.0 || b == 0.0) return 0.0;
  else return 1/(1/a + 1/b);
}

int Math::degree( int dof) {
  int p = 0;
  while( p*(p+1)/2 <= dof) {
    //cout << p*(p+1)/2 << " " << dof << endl;
    p++;
  }
  return p-2;
}

int Math::triNum( int n) { return (n+2)*(n+1)/2; }


int Math::pow2roundup( int n) {
  if (n < 0) {
    return 0;
  }
  --n;
  n |= n >> 1;
  n |= n >> 2;
  n |= n >> 4;
  n |= n >> 8;
  n |= n >> 16;
  return n+1;
}
