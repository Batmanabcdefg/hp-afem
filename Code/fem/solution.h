#pragma once
#include <iostream>
#include <vector>

#include "rhs.h"
#include "../dofs.h"
#include "../element.h"
#include "../matrix.h"
#include "../piecewisepolynomial.h"
#include "../triangleset.h"

class Solver;

namespace FEM {

class Solution : public PiecewisePolynomial {
 public:
  Solution() = default;
  Solution(Vector &sol, const ElementDofsSet &elts, const Rhs &rhs);
  Solution(Vector &sol, const ElementDofsMap &elts, const Rhs &rhs);

  void insert_dofs(Element *elt, Dofs local, bool definedOn);
  void erase(Element *elt);

  const Vector& global() const { return _global; }
  Vector global() { return _global; }
  scalar operator[](int i) { return _global[i]; }

  int numDOFs() { return _global.size(); }

  static Solution fromSolFile(std::string filename, ElementDofsSet elts, Rhs &rhs);

  Vector residualVector();
  scalar residual() { return residualVector().norm(); }

  StiffnessMatrix _systemMat;
  Vector _systemRhs;

 private:
  Vector _global;
  ElementDofsMap _g;

  Rhs _rhs;

  Vector indexSol(const Dofs &g);
};

}
