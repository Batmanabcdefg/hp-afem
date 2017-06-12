#pragma once

#include <stdio.h>
#include <iostream>
#include "basis.h"
#include "tritype.h"
#include "basisfunction.h"

using namespace std;
using namespace GiNaC;

namespace Precomputations {
class BezierBasisFunction : public BasisFunction {
 public:
  BezierBasisFunction(size_t deg, const std::vector<numeric>& pt)
    : BasisFunction(deg), _pt(pt) {
      _ex = multinomial() * GiNaC::pow(L1(), pt[0]) * GiNaC::pow(L2(), pt[1]) * GiNaC::pow(L3(), pt[2]);
    }
 private:
  std::vector<numeric> _pt;

  ex multinomial() {
    ex result = factorial(_deg)/(factorial(_pt[0]) * factorial(_pt[1]) * factorial(_pt[2]));
    cout << "nu punt " << result << endl;
    return result;
  }
};

class BezierBasis : public Basis {
protected:
  std::vector<std::vector<numeric>> points(size_t deg) const {
    std::vector<std::vector<numeric> > p;

    for( size_t i = 0; i <= deg; i++) {
      for( size_t j = 0; j+i <= deg; j++) {
        size_t k = deg - i - j;
        p.push_back({i,j,k});
      }
    }

    return p;
  }

public:
  BezierBasis(size_t deg) : Basis(deg) {
    auto pts = points(deg);
    for (const auto &pt : pts) {
      _bfs.push_back(BezierBasisFunction(deg, pt));
    }

    for (auto &bf : _bfs) cout << bf.bf() << endl;
  }
  BezierBasis & operator=(BezierBasis &other) = default;
};
}
