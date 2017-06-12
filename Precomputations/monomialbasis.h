#pragma once

#include <stdio.h>
#include "basis.h"
#include "basisfunction.h"

namespace Precomputations {
class MonomialBasis : public Basis {
protected:
public:
  MonomialBasis(size_t deg) : Basis(deg) {
    for( size_t xpow = 0; xpow <= deg; xpow++) {
      for( size_t ypow = 0; xpow + ypow <= deg; ypow++) {
        _bfs.push_back(MonomialBasisFunction(xpow, ypow));
      }
    }
  }
};
}
