#pragma once

#include <assert.h>
#include <cmath>
#include <stdexcept>
#include <stdio.h>
#include <ginac/ginac.h>

using namespace std;
using namespace GiNaC;

namespace Precomputations {
ex LegendrePoly(const symbol &x, int n) {
  ex PKer = pow(x * x - 1, n);
  return normal(1 / (GiNaC::pow(2, n) * factorial(n)) * diff(PKer, x, n));
}

class BasisFunction {
public:
  symbol x;
  symbol y;
  size_t _deg;
  ex _ex;

  ex L1() { return 1 - x - y; }
  ex L2() { return x; }
  ex L3() { return y; }

  ex bf() { return _ex; }

  ex sub(ex a, ex b) {
    return _ex.subs(lst(x == a, y == b));
  }

  ex sub(lst l) {
    return _ex.subs(lst(x, y), l);
  }

  ex integrate(ex integrand) {
    return integral(x, 0, 1, integral(y, 0, 1-x, integrand).eval_integ()).eval_integ().normal();
  }

  ex a(BasisFunction &that) {
    ex you = that.sub( lst(x,y));
    ex integrand = diff(_ex, x)*diff(you, x) + diff(_ex, y)*diff(you, y);
    return integrate(integrand);
  }

  ex a0(BasisFunction &that) {
    ex you = that.sub( lst(x,y));
    ex integrand = diff(_ex, x)*diff(you, x);
    return integrate(integrand);
  }

  ex a1(BasisFunction &that) {
    ex you = that.sub( lst(x,y));
    ex integrand = diff(_ex, x)*diff(you, y) + diff(_ex, y)*diff(you, x);
    return integrate(integrand);
  }

  ex a2(BasisFunction &that) {
    ex you = that.sub( lst(x,y));
    ex integrand = diff(_ex, y)*diff(you, y);
    return integrate(integrand);
  }

  ex m(BasisFunction &that) {
    ex you = that.sub( lst(x,y));
    ex integrand = _ex * you;
    return integrate(integrand);
  }

  ex energyNorm() {
    return sqrt(a(*this));
  }

  BasisFunction() {}
  BasisFunction(size_t deg) : x("x"), y("y"), _deg(deg) {}
};

class VertexBasisFunction : public BasisFunction {
protected:
  size_t _v;
public:
  VertexBasisFunction(size_t v) : BasisFunction(1), _v(v) {
    switch(v) {
      case 1:
        _ex = L1();
        break;
      case 2:
        _ex = L2();
        break;
      case 3:
        _ex = L3();
        break;
      default:
        throw runtime_error("vertex must be between 1 and 3");
    }
  }
};

class EdgeBasisFunction : public BasisFunction {
protected:
  size_t _e;
  bool _reverse;
public:
  EdgeBasisFunction(size_t e, size_t deg, bool reverse) : 
    BasisFunction(deg), _e(e), _reverse(reverse) {
    size_t j1 = e;
    size_t j2 = 1 + (e % 3);
    ex L[] = {L1(), L2(), L3()};
    ex Lj1 = L[j1-1], Lj2 = L[j2-1];
    ex arg = (GiNaC::pow(-1,reverse))*(Lj2 - Lj1);
    symbol z("z");

    ex diff_leg = diff(LegendrePoly(z, deg), z).subs(z == arg);
    ex E = -8 * sqrt(ex(4*deg + 2))/(deg * (deg + 1)) * diff_leg;
    _ex = (Lj1 * Lj2 * E).normal();
  }
};

class FaceBasisFunction : public BasisFunction {
protected:
  size_t _r1, _r2;
public:
  FaceBasisFunction() {}
  FaceBasisFunction(size_t r1, size_t r2) : 
    BasisFunction(r1+r2+3), _r1(r1), _r2(r2) {
    ex F_r1r2 = pow((L2() - L1()), r1) * pow((2 * L3() - 1), r2);
    _ex = (L1() * L2() * L3() * F_r1r2).normal();
  }
};

class MonomialBasisFunction : public BasisFunction {
protected:
  size_t _xpow, _ypow;
public:
  MonomialBasisFunction(size_t xpow, size_t ypow) : 
    BasisFunction(xpow+ypow), _xpow(xpow), _ypow(ypow) {
    _ex = pow(x, xpow) * pow(y, ypow);
  }
};
}
