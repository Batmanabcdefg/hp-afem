#pragma once

#include "element.h"
#include "matrix.h"

class Poly {
public:
  /**
   *  If we view `v` as the coefficients of a polynomial expressed in the basis
   *  locally at `elt`, then we see:
   *    u(x,y) = v^\top \Phi_elt(x,y)
   *    => | u |^2_{H^1_0(elt)} = v^\top A^K v,
   *  where A^K is the local stiffness matrix.
   */
  static scalar squaredH1Norm(Element *elt, const Vector &v);

  /**
   *  If we view `v` as the coefficients of a polynomial expressed in the basis
   *  locally at `elt`, then we see:
   *    u(x,y) = v^\top \Phi_elt(x,y)
   *    => | u |^2_{L^2(elt)} = v^\top M^K v,
   *  where A^K is the local mass matrix.
   */
  static scalar squaredL2Norm(Element *elt, const Vector &v);

  /**
   *  If we view `v` as the coefficients of a polynomial expressed in the basis
   *  locally at `elt`, then we see that
   *    u(x,y) = v^\top \Phi_elt(x,y)
   *    => \int_elt u = v^\top \int_elt \Phi_elt = v^\top E_elt,
   *  where E_elt is the local element vector.
   */
  static scalar area(Element *elt, const Vector &v);

  /**
   *  If we view `v` as the coefficients of a polynomial expressed in the basis
   *  locally at `elt`, then we see that on the children of `elt`, `v` is still
   *  a polynomial.  This function finds the coefficients of this polynomial
   *  expressed in the childrens' basis.
   */
  static Vector onChild(Element *elt, const Vector &v, bool rightChild);

  /**
   * Simple method to find the difference of two Vectors that are not necessarily
   * of the same length.
   */
  static Vector minus(const Vector &v1, const Vector &v2);

  /**
   * Simple method to find the sum of two Vectors that are not necessarily
   * of the same length.
   */
  static Vector plus(const Vector &v1, const Vector &v2);
};

