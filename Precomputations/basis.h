#pragma once

#include <stdio.h>
#include <vector>
#include <assert.h>
#include "basisfunction.h"
#include "ginacsupportforeigen.h"
#include "triangle.h"

using namespace std;
using namespace GiNaC;

namespace Precomputations {
class Basis {
protected:
  size_t _deg;
  vector<BasisFunction> _bfs;
  Triangle _tri;

  matrix _A0, _A1, _A2, _M, _b;

public:
  size_t dim() const {
    return (_deg + 1)*(_deg + 2)/2;
  }

  void transformTo(const Triangle &tri) {
    _tri = tri;
    for(size_t i = 0; i < dim(); i++) {
      symbol x = _bfs[i].x;
      symbol y = _bfs[i].y;
      _bfs[i]._ex = _bfs[i].sub(tri.inverse(x, y));
    }
    cout << "hier nu" << endl;
    for( auto &bf : _bfs) cout << bf.bf() << endl;
  }

  void genAs() {
    size_t d = dim();
    _A0 = matrix(d, d);
    _A1 = matrix(d, d);
    _A2 = matrix(d, d);
    _M = matrix(d, d);
    _b = matrix(d, 1);
    for( size_t i = 0; i < d; i++) {
      for( size_t j = 0; j < d; j++) {
        if( i <= j) {
          _A0(i, j) = _bfs[i].a0(_bfs[j]);
          _A1(i, j) = _bfs[i].a1(_bfs[j]);
          _A2(i, j) = _bfs[i].a2(_bfs[j]);
          _M(i, j) = _bfs[i].m(_bfs[j]);
        } else {
          _A0(i, j) = _A0(j, i);
          _A1(i, j) = _A1(j, i);
          _A2(i, j) = _A2(j, i);
          _M(i, j) = _M(j, i);
        }
      }
      ex b = 0;
      for( size_t j = 0; j < 3; j++) {
        b = b + _M.op(i*d + j);
      }
      _b.set(i, 0, b.normal());
    }
  }

  matrix A0() { if( _A0.rows() != dim()) genAs(); return _A0; }
  matrix A1() { if( _A1.rows() != dim()) genAs(); return _A1; }
  matrix A2() { if( _A2.rows() != dim()) genAs(); return _A2; }
  matrix M() { if( _M.rows() != dim()) genAs(); return _M; }
  matrix b() { if( _M.rows() != dim()) genAs(); return _b; }

  matrix genFPhi_ks( Basis &b) {
    size_t d = dim();
    matrix FPhi_ks(d,d);
    auto pts = b._tri.points(_deg);
    for( size_t i = 0; i < d; i++) {
      for( size_t j = 0; j < d; j++) {
        auto phi_k_j = b.bf(j).sub(pts[i][0], pts[i][1]);
        FPhi_ks(i,j) = phi_k_j;
      }
    }
    return FPhi_ks;
  }

  matrix genFPhis( Basis &b) {
    size_t d = dim();
    matrix FPhis(d,d);
    auto pts = b._tri.points(_deg);
    for( size_t i = 0; i < d; i++) {
      for( size_t j = 0; j < d; j++) {
        auto phi_j = bf(j).sub(pts[i][0], pts[i][1]);
        FPhis(i,j) = phi_j;
      }
    }
    return FPhis;
  }
  
  matrix transferMatrixs( Basis &b) {
    size_t d = dim();
    matrix FPhi_ks = genFPhi_ks(b);
    matrix FPhis = genFPhis(b);
    matrix vars(d,d);
    for( size_t i = 0; i < d; i++) for( size_t j = 0; j < d; j++) vars(i,j) = symbol();

    cout << "Gonna solve system of size " << d << endl;
    matrix sol = FPhi_ks.solve( vars, FPhis);
    return sol;
  }

  Eigen::MatrixXd transferMatrix( Basis &b) {
    size_t d = dim();
    Eigen::MatrixXd FPhi_k(d, d);
    Eigen::MatrixXd FPhi(d, d);
    auto pts = b._tri.points(_deg);
    for( size_t i = 0; i < d; i++) {
      for( size_t j = 0; j < d; j++) {
        auto phi_k_j = b.bf(j).sub(pts[i][0], pts[i][1]);
        auto phi_j = bf(j).sub(pts[i][0], pts[i][1]);
        FPhi_k(i, j) = Eigen::ex2double(phi_k_j);
        FPhi(i, j) = Eigen::ex2double(phi_j);
      }
    }

    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> dec(FPhi_k);
    return dec.solve(FPhi);
  }

  void printBasis() {
    size_t d = dim();
    for( size_t j = 0; j < d; j++) {
      cout << bf(j).bf() << endl;
    }
    cout << endl;
  }

  BasisFunction &bf(size_t i) { assert(i < dim()); return _bfs[i]; }

  Basis(size_t deg) : _deg(deg) {}
  Basis &operator=(Basis &other) = default;
};
}
