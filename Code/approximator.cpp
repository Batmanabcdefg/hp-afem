#include <Eigen/Dense>
#include <iostream>
#include <limits>

#include "approximator.h"
#include "poly.h"
#include "degree.h"

using namespace std;

scalar Approximator::getArea(Element *elt) {
  auto it = _areas.find(elt);
  assert(it != _areas.end());
  return it->second;
}

scalar Approximator::setArea(Element *elt, scalar val) {
  _areas.insert(make_pair(elt, val));
  return val;
}

Matrix Approximator::approximationMatrix(Element *elt, int dim) {
  Matrix A(dim+1, dim+1);
  Vector eltvec = elt->elementVector(dim);
  A.block(0, 0, dim, dim) = elt->elementMatrix(dim);
  A.block(0, dim, dim, 1) = eltvec;
  A.block(dim, 0, 1, dim) = eltvec.transpose();
  A(dim, dim) = 0;
  return A;
}

void Approximator::setAreaRecursive(Element *node) {
  if( _sol.has(node)) {
    setArea(node, Poly::area(node, _sol.locallyAt(node)));
  } else {
    setAreaRecursive(node->left());
    setAreaRecursive(node->right());
    setArea(node, getArea(node->left()) + getArea(node->right()));
  }
}

Vector Approximator::approximationRhsRecursive(Element *node, int dim) {
  Vector rhs = Vector::Zero(dim+1);
  if (_sol.has(node)) {
    int soldim = _sol.local_dim(node);
    Vector localSol = Vector::Zero(max(dim, soldim));
    localSol.head(soldim) = _sol.locallyAt(node);
    rhs.head(dim) = (node->elementMatrix(max(dim,soldim))*localSol).head(dim);
    // TODO: is this correct?
    setArea(node, Poly::area(node, localSol));
  } else {
    assert(!node->isLeaf());

    // Fetch child vectors.
    Vector leftrhs = approximationRhsRecursive(node->left(), dim);
    Vector rightrhs = approximationRhsRecursive(node->right(), dim);

    setArea(node, getArea(node->left()) + getArea(node->right()));

    // Combine them to get the RHS for node.
    int leftdim = leftrhs.rows()-1;
    auto ltmt = node->transferMatrix(0, leftdim);
    ltmt.transposeInPlace();
    auto leftrhsTransfered = ltmt * leftrhs.head(leftdim);
    rhs.head(leftdim) += leftrhsTransfered;

    int rightdim = rightrhs.rows()-1;
    auto rtmt = node->transferMatrix(1, rightdim);
    rtmt.transposeInPlace();
    auto rightrhsTransfered = rtmt * rightrhs.head(rightdim);
    rhs.head(rightdim) += rightrhsTransfered;
  }

  // We need to fill in the last element separately.
  rhs[dim] = getArea(node);
  return rhs;
}

Vector Approximator::solveSystem(Matrix &mat, Vector &rhs) {
  Eigen::ColPivHouseholderQR<Matrix> solver(mat);
  Vector res = solver.solve(rhs);
  return res.head(res.rows()-1);
}

scalar Approximator::error(Element *elt, int dof) {
  Vector x;
  int dim = Degree::dofToDim(dof);
  // case where we have a zeroth order approximation; just return the H1 norm.
  if (dim == 0) {
    //@TODO THIS CONTAINS AN ERROR
    return _sol.squaredH1Norm(elt);
  } else if( dim == 1) {
    // We solve the 2x2 system
    // [[<1,1>_K, int_K 1]  [c       = [<w,1>_K
    //  [int_K 1, 0      ]]  lambda] =  int_K w]
    //
    // now, <1,1>_K and <w,1>_K must be zero, and int_K 1 = vol(K), so the 
    // solution is
    //     c = int_K w / vol(K).
    // Now, we know that phi_1 + phi_2 + phi_3 = 1 (partition of unity), so that
    // setting the vector solution c = [c c c] has the desired property.
    setAreaRecursive(elt);
    scalar c = getArea(elt)/elt->vol();
    x = Vector::Zero(3);
    x[0] = c;
    x[1] = c;
    x[2] = c;
  } else {
    // If we have a high order approximation of a polynomial, the error is zero.
    if (_sol.has(elt) && dim >= _sol.local_dim(elt)) {
      return 0.0;
    }
    assert(dim <= elt->basis()->p());

    Vector b = approximationRhsRecursive(elt, dim);
    Matrix A = approximationMatrix(elt, dim);
    x = solveSystem(A, b);
  }

  // Transfer approximation to the elements the solution is available on.
  PiecewisePolynomial appr;
  appr.insert_vector(elt, x, true);
  // @TODO THIS CONTAINS AN ERROR
  scalar res = appr.squaredH1NormOfDifferenceWith(elt, _sol);

  if( res < -numeric_limits<scalar>::epsilon()) {
    cerr << res << endl;
    assert(false);
  }
  return fmax(0.0, res);
}

Approximator::Approximator(PiecewisePolynomial &sol) : _sol(sol) {}
