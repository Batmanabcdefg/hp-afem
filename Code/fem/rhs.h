#pragma once

#include "../hasdof.h"
#include "../piecewisepolynomial.h"
#include "../triangleset.h"

namespace FEM {

class Rhs : public PiecewisePolynomial, public HasDOF {
 public:
  Rhs() = default;
  Rhs(int dof) : _dof(dof) {}

  virtual int dof() const override { return _dof; }
  void retainOnly(ElementSet &elts);

 protected:
  int _dof;
};

}
