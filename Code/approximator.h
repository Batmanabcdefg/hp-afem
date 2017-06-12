#pragma once
#include <algorithm>
#include <map>

#include "piecewisepolynomial.h"
#include "element.h"
#include "matrix.h"

// Approximator class.  This implements the Binev hp-error functional.  Given
// a partition, an element K, a solution w, and a degree p, we find the 
// best-approximating polynomial Q_p(w) to w (in H^1_0-seminorm) on K.
//
// TODO: make this more efficient, in the sense that it can do multiple error
// queries (different K and p) for a constant partition and solution.
class Approximator {
public:
  Approximator(PiecewisePolynomial &sol);
  scalar error(Element *elt, int dof);

protected:
  // The solution we are approximating.
  PiecewisePolynomial &_sol;

  // Integral of the solution locally at these elements; if _sol is polynomial
  // on the element, this is computed using Poly::area().  Else, it is the 
  // summation of such values.
  std::map<Element *, scalar> _areas;

  // Getter and setter for the above data structure.
  scalar getArea(Element *elt);
  scalar setArea(Element *elt, scalar val);
  void setAreaRecursive(Element *node);

  // Computes the matrix necessary for finding the best approximant.
  Matrix approximationMatrix(Element *elt, int dim);
  
  // Computes the RHS necessary for finding the best approximant.
  Vector approximationRhsRecursive(Element *parent, int dim);

  // Solves a system given matrix and right-hand side using Householder QR.
  Vector solveSystem(Matrix &mat, Vector &rhs);
};
