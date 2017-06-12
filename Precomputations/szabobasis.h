#pragma once

#include <stdio.h>
#include <iostream>
#include "basis.h"
#include "tritype.h"
#include "basisfunction.h"

using namespace std;
using namespace GiNaC;

namespace Precomputations {
class SzaboBasis : public Basis {
protected:
  Tritype _tt;

  size_t pairCantor( size_t x, size_t y) { return (x+y)*(x+y+1)/2 + y; }

  vector<FaceBasisFunction> generateFaceBasisFunctions() {
    vector<FaceBasisFunction> ff;
    if( _deg <= 2) return ff;
    size_t numFuncs = (_deg - 2)*(_deg - 1)/2;

    ff.resize(numFuncs);
    for( size_t r1 = 0; r1 < _deg - 2; r1++) {
      for( size_t r2 = 0; r1 + r2 < _deg - 2; r2++) {
        ff[pairCantor(r1, r2)] = FaceBasisFunction(r1, r2);
        cout << ff[pairCantor(r1, r2)].bf() << endl;
      }
    }

    return ff;
  }

public:
  SzaboBasis(size_t deg, Tritype tt) : Basis(deg), _tt(tt) {
    for( size_t i = 0; i < 3; i++) {
      _bfs.push_back(VertexBasisFunction(i+1));
    }

    auto facefuncs = generateFaceBasisFunctions();

    size_t j = 0;
    if( deg >= 2) {
      for( size_t d = 1; d <= deg; d++) {
        if( d >= 3) for( size_t r1 = 0; r1 < d-2; r1++) 
            _bfs.push_back(facefuncs[j++]);

        if( d < deg) for( size_t e = 0; e < 3; e++) 
            _bfs.push_back(EdgeBasisFunction(e+1, d, tt.i(e)));
      }
    }

    for( auto &bf : _bfs) cout << bf.bf() << endl;
  }
  SzaboBasis & operator=(SzaboBasis &other) = default;
};
}
