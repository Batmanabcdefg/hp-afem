#pragma once

#include "solution.h"
#include "rhs.h"
#include "../triangleset.h"
#include "../dofhandler.h"

namespace FEM {
class Solver {
public:
  Solver(const ElementDofsMap &elts, Rhs &rhs);
  Solver(const DOFHandler &handler, Rhs &rhs) : Solver(handler.map(), rhs) {}

  int numDOFs();

  const Solution        &sol()          const { return _sol; }
  const StiffnessMatrix &systemMatrix() const { return _sysmat; }
  const LoadVector      &systemRhs()    const { return _sysrhs; }

protected:
  Solution _sol;
  const ElementDofsMap &_elts;
  const Rhs &_femrhs;
  int _numDOFs = -1;

  LoadVector _sysrhs;
  StiffnessMatrix _sysmat;

  int estimateNonZeros();
  void computeSystemMatrix();
  void computeSystemRhs();
  void solveSystem();
};
}
