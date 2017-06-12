#include <iostream>
#include <cassert>

#include "poly.h"
#include "print.h"
#include "piecewisepolynomial.h"
#include "poly.h"

using namespace std;

void PiecewisePolynomial::insert_vector(Element *elt, Vector local,
    bool definedOn) {
  assert(l2g.count(elt) == 0);
  l2g.insert(make_pair(elt, local));
  if (definedOn) {
    _definedOn.insert(elt);
  }
  _availableOn.insert(elt);
}

void PiecewisePolynomial::erase(Element *elt) {
  auto it = l2g.find(elt);
  assert(it != l2g.end());

  _availableOn.erase(elt);

  // we aren't _defined_ on all nodes that we are available on
  _definedOn.erase(elt);
  l2g.erase(it);
}

const Vector &PiecewisePolynomial::locallyAt(Element *elt) const {
  auto it = l2g.find(elt);
  assert(it != l2g.end());
  return it->second;
}

void PiecewisePolynomial::copyToChildren(Element *parent) {
  assert(!parent->isLeaf());
  Vector poly = locallyAt(parent);

  Element *left = parent->left();
  Element *right = parent->right();

  insert_vector(left, Poly::onChild(parent, poly, 0), false);
  insert_vector(right, Poly::onChild(parent, poly, 1), false);
}

ElementSet PiecewisePolynomial::copyToRecursive(const PiecewisePolynomial &other, Element *elt) {
  ElementSet otherIsAvailableOn;
  if (!other.has(elt)) {
    assert(!elt->isLeaf());
    // our polynomial is available on our children, skip this
    if (!has(elt->left())) {
      copyToChildren(elt);
    }
    otherIsAvailableOn = TriangleSet::unionSets(
        copyToRecursive(other, elt->left()), 
        copyToRecursive(other, elt->right())
    );
  } else {
    otherIsAvailableOn.insert(elt);
  }
  return otherIsAvailableOn;
}

ElementSet PiecewisePolynomial::copyTo(const PiecewisePolynomial &other) {
  ElementSet otherIsAvailableOn;
  for (auto &elt : _definedOn) {
    otherIsAvailableOn =
        TriangleSet::unionSets(otherIsAvailableOn, copyToRecursive(other, elt));
  }

  return otherIsAvailableOn;
}

ElementScalarSet PiecewisePolynomial::squaredH1NormsOfDifferenceWith(Element *elt, PiecewisePolynomial &other) {
  ElementScalarSet ret;
  ElementSet elts = copyToRecursive(other, elt);

  for (auto &elt : elts) {
#if GALERKIN_ORTH
     auto error = Poly::squaredH1Norm(elt, other.locallyAt(elt))
                - Poly::squaredH1Norm(elt, locallyAt(elt));
#else
    auto curpoly = Poly::minus(locallyAt(elt), other.locallyAt(elt));
    auto error = Poly::squaredH1Norm(elt, curpoly);
#endif
    ret.insert(make_pair(elt, error));
  }

  return ret;
}

ElementScalarSet PiecewisePolynomial::squaredL2Norms(Element *elt) {
  ElementScalarSet ret;
  if(has(elt)) {
    ret.insert(make_pair(elt, Poly::squaredL2Norm(elt, locallyAt(elt))));
  } else {
    assert(!elt->isLeaf());

    //recurse left
    ElementScalarSet leftret = squaredL2Norms(elt->left());
    ret.insert(leftret.begin(), leftret.end());

    //recurse right
    ElementScalarSet rightret = squaredL2Norms(elt->right());
    ret.insert(rightret.begin(), rightret.end());
  }
  return ret;
}

ElementScalarSet PiecewisePolynomial::squaredH1Norms(Element *elt) {
  ElementScalarSet ret;
  if(has(elt)) {
    ret.insert(make_pair(elt, Poly::squaredH1Norm(elt, locallyAt(elt))));
  } else {
    assert(!elt->isLeaf());

    //recurse left
    ElementScalarSet leftret = squaredH1Norms(elt->left());
    ret.insert(leftret.begin(), leftret.end());

    //recurse right
    ElementScalarSet rightret = squaredH1Norms(elt->right());
    ret.insert(rightret.begin(), rightret.end());
  }
  return ret;
}

scalar PiecewisePolynomial::squaredH1NormOfDifferenceWith(Element *elt,
    PiecewisePolynomial &other) {
#if GALERKIN_ORTH
  return other.squaredH1Norm(elt) - squaredH1Norm(elt);
#else
  scalar res = 0;
  for (auto &child : squaredH1NormsOfDifferenceWith(elt, other)) {
    res += child.second;
  }
  return res;
#endif
}

scalar PiecewisePolynomial::squaredL2Norm(Element *elt) {
  scalar res = 0;
  for (auto &child : squaredL2Norms(elt)) {
    res += child.second;
  }
  return res;
}

scalar PiecewisePolynomial::squaredH1Norm(Element *elt) {
  scalar res = 0;
  for (auto &child : squaredH1Norms(elt)) {
    res += child.second;
  }
  return res;
}



ElementScalarSet PiecewisePolynomial::squaredH1NormsOfDifferenceWith(PiecewisePolynomial &other) {
  ElementScalarSet ret;
  ElementSet elts = copyTo(other);
  if (elts.size()) {
    // compute the error on this approximation

    for (auto &child : elts) {
#if GALERKIN_ORTH
       auto error = Poly::squaredH1Norm(child, other.locallyAt(child))
                  - Poly::squaredH1Norm(child, locallyAt(child));
#else
      auto curpoly = Poly::minus(locallyAt(child), other.locallyAt(child));
      auto error = Poly::squaredH1Norm(child, curpoly);
#endif
      ret.insert(make_pair(child, error));
    }
  } else {
    for (auto &child : other.definedOn()) {
      ret.insert(make_pair(child, Poly::squaredH1Norm(child,
                                                      other.locallyAt(child))));
    }
  }

  return ret;
}

ElementScalarSet PiecewisePolynomial::squaredL2Norms() {
  ElementScalarSet ret;
  for (auto &child : _definedOn) {
    auto curnorm = Poly::squaredL2Norm(child, locallyAt(child));
    ret.insert(make_pair(child, curnorm));
  }
  return ret;
}

ElementScalarSet PiecewisePolynomial::squaredH1Norms() {
  ElementScalarSet ret;
  for (auto &child : _definedOn) {
    ret.insert(make_pair(child, Poly::squaredH1Norm(child, locallyAt(child))));
  }
  return ret;
}

scalar PiecewisePolynomial::squaredH1NormOfDifferenceWith(PiecewisePolynomial &other) {
#if GALERKIN_ORTH
  return other.squaredH1Norm() - squaredH1Norm();
#else
  scalar res = 0;
  for (auto &child : squaredH1NormsOfDifferenceWith(other)) {
    res += child.second;
  }
  return res;
#endif
}

scalar PiecewisePolynomial::squaredL2Norm() {
  scalar res = 0;
  for( auto &child : squaredL2Norms()) res += child.second;
  return res;
}

scalar PiecewisePolynomial::squaredH1Norm() {
  scalar res = 0;
  for( auto &child : squaredH1Norms()) res += child.second;
  return res;
}

ostream &PiecewisePolynomial::printForFile(const ElementSet &on, ostream &os) {
  os << Print::formatted("%lu %d", on.size(), maximum_dim(on)) << endl;
  for (auto &elt : on) {
    Vector vec = locallyAt(elt);
    for (int i = 0; i < vec.rows(); i++) {
      if (i > 0) os << " ";
      os << Print::formatted("%.30Lf", (long double) vec[i]);
    }
    os << endl;
  }

  return os;
}

ostream &PiecewisePolynomial::print(const ElementSet &on, ostream &os) {
  os << Print::formatted("%lu", on.size()) << endl;
  for (auto &elt : on) {
    Vector vec = locallyAt(elt);
    os << Print::formatted("%lu %d %d %d %d", vec.rows(), elt->i(0), elt->i(1), elt->i(2), elt->type().toInt());
    for (int i = 0; i < vec.rows(); i++) {
      os << Print::formatted(" %.16f", (double) vec[i]);
    }
    os << endl;
  }

  return os;
}
