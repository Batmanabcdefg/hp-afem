#include <iostream>

#include "elementmatrices.h"
#include "element.h"

using namespace std;

Matrix &ElementMatrices::get(int tt, int tc) {
  if( _eltmats[tt][tc].rows() > 0) return _eltmats[tt][tc];

  Element *elt = _root;
  if( tc == 1) elt = _root->left();
  if( tc == 2) elt = _root->right();
  if( tc == 3) elt = _root->left()->left();
  assert(elt != nullptr);

  scalar D = 2.0L*elt->vol();
  scalar x2mx0 = (elt->v(2)->x - elt->v(0)->x), 
         x1mx0 = (elt->v(1)->x - elt->v(0)->x),
         y2my0 = (elt->v(2)->y - elt->v(0)->y), 
         y1my0 = (elt->v(1)->y - elt->v(0)->y);
  scalar E0 = x2mx0*x2mx0 + y2my0*y2my0,
         E1 = x1mx0*x2mx0 + y1my0*y2my0,
         E2 = x1mx0*x1mx0 + y1my0*y1my0;
  Matrix eltmat = (_bases.basis(tt).eltmat(0)*E0 -
                   _bases.basis(tt).eltmat(1)*E1 +
                   _bases.basis(tt).eltmat(2)*E2)/D;

  Vector v = Vector::Random(eltmat.rows());
  return _eltmats[tt][tc] = eltmat;
}
