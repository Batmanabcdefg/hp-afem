#include <cassert>
#include <algorithm>

#include "poly.h"

using namespace std;

Vector Poly::onChild(Element *elt, const Vector &v, bool rightChild) {
  assert(!elt->isLeaf());
  Vector res = elt->transferMatrix(rightChild, v.rows()) * v;
  return res;
}

Vector Poly::minus(const Vector &v1, const Vector &v2) {
  int vl = v1.rows(), ol = v2.rows();
  Vector res = Vector::Zero(max(vl, ol));

  if( vl < ol) {
    res.tail(ol - vl) = -v2.tail(ol - vl);
    res.head(vl) = v1 - v2.head(vl);
  } else if( vl > ol) {
    res.tail(vl - ol) = v1.tail( vl - ol);
    res.head(ol) = v1.head(ol) - v2;
  } else {
    res = v1 - v2;
  }

  return res;
}

Vector Poly::plus(const Vector &v1, const Vector &v2) {
  int vl = v1.rows(), ol = v2.rows();
  Vector res = Vector::Zero(max(vl, ol));

  if( vl < ol) {
    res.tail(ol - vl) = v2.tail(ol - vl);
    res.head(vl) = v1 + v2.head(vl);
  } else if( vl > ol) {
    res.tail(vl - ol) = v1.tail( vl - ol);
    res.head(ol) = v1.head(ol) + v2;
  } else {
    res = v1 - v2;
  }

  return res;
}

scalar Poly::squaredH1Norm(Element *elt, const Vector &v) {
  return v.dot( elt->elementMatrix(v.rows()) * v);
}

scalar Poly::squaredL2Norm(Element *elt, const Vector &v) {
  return v.dot( elt->massMatrix(v.rows()) * v);
}

scalar Poly::area(Element *elt, const Vector &v) {
  return v.dot( elt->elementVector(v.rows()));
}
