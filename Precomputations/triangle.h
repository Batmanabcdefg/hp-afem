#pragma once

#include <ginac/ginac.h>
#include <vector>

using namespace GiNaC;

namespace Precomputations {
class Triangle {
protected:
  std::vector<std::vector<numeric> > _pts;

public:
  numeric vol() const {
    auto v = _pts;
    auto twice = -v[1][0]*v[0][1] - v[2][0]*v[0][1] +
                  v[0][0]*v[1][1] - v[2][0]*v[1][1] -
                  v[0][0]*v[2][1] + v[1][0]*v[2][1];
    return abs(twice/2);
  }

  std::vector<std::vector<numeric>> points(size_t deg) const {
    auto v = _pts;
    std::vector<std::vector<numeric> > p;

    for( size_t i = 0; i <= deg; i++) {
      for( size_t j = 0; j+i <= deg; j++) {
        size_t k = deg - i - j;
        numeric p0 = (k*v[0][0] + j*v[1][0] + i*v[2][0])/deg;
        numeric p1 = (k*v[0][1] + j*v[1][1] + i*v[2][1])/deg;
        p.push_back({p0, p1});
      }
    }

    return p;
  }

  lst inverse( ex x, ex y) const {
    auto v = _pts;
    matrix m(2, 2);
    m = v[1][0] - v[0][0], v[2][0] - v[0][0],
        v[1][1] - v[0][1], v[2][1] - v[0][1];
    matrix b(2, 1, lst(x - v[0][0], y - v[0][1]));
    matrix minv = m.inverse();
    matrix res = minv.mul(b);
    return lst(res(0,0), res(1,0));
  }

  std::vector<std::vector<numeric> > pts() { return _pts;}

  Triangle(const std::vector<std::vector<numeric> > &pts) : _pts(pts) {}
  Triangle() {
    std::vector<numeric> p1 = {0, 0};
    std::vector<numeric> p2 = {1, 0};
    std::vector<numeric> p3 = {0, 1};
    std::vector<std::vector<numeric> > pts = {p1, p2, p3};
    _pts = pts;
  }
};
}
