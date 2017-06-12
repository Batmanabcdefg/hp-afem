#pragma once
#include <map>
#include <fstream>
#include <iostream>

#include "triangleset.h"
#include "element.h"

class PiecewisePolynomial {
public:
  virtual void insert_vector(Element *elt, Vector local, bool definedOn);
  virtual void erase(Element *elt);
  virtual const Vector &locallyAt(Element *elt) const;

  const ElementSet &definedOn() const { return _definedOn; }
  const ElementSet &availableOn() const { return _availableOn; }
  void copyToChildren(Element *parent);

  bool has(Element *elt) const { return l2g.count(elt) > 0; }
  ElementPairMap<Vector> values() { return l2g; }

  virtual int local_dim(Element *elt) {
    assert(has(elt));
    return l2g.find(elt)->second.rows();
  }

  int maximum_dim(const ElementSet& on) {
    int maximum = 0;
    for (const auto elt : on) {
      if (local_dim(elt) > maximum) {
        maximum = local_dim(elt);
      }
    }
    return maximum;
  }

  /**
   *  If we know the partition on which `other' is defined, is refined relative
   *  to our partition, this works. Otherwise, bad stuff may happen.
   */
  ElementSet copyTo(const PiecewisePolynomial &other);

  // compute the norm of this polynomial on a specific triangle
  ElementScalarSet squaredH1NormsOfDifferenceWith(Element *elt, PiecewisePolynomial &other);
  ElementScalarSet squaredL2Norms(Element *elt);
  ElementScalarSet squaredH1Norms(Element *elt);

  // sum the above values
  scalar squaredH1NormOfDifferenceWith(Element *elt, PiecewisePolynomial &other);
  scalar squaredL2Norm(Element *elt);
  scalar squaredH1Norm(Element *elt);



  // compute the norm of this polynomial on all triangles it is defined on
  ElementScalarSet squaredH1NormsOfDifferenceWith(PiecewisePolynomial &other);
  ElementScalarSet squaredL2Norms();
  ElementScalarSet squaredH1Norms();

  // sum the above values
  scalar squaredH1NormOfDifferenceWith(PiecewisePolynomial &other);
  scalar squaredL2Norm();
  scalar squaredH1Norm();

  std::ostream &print(std::ostream &os = std::cout) { return print(_definedOn, os); }
  std::ostream &print(const ElementSet &on, std::ostream &os = std::cout);
  std::ostream &printForFile(const ElementSet &on, std::ostream &os = std::cout);

protected:
  /* local to global converter */
  ElementPairMap<Vector> l2g;

  /* set of elements we are defined on */
  ElementSet _definedOn;

  ElementSet _availableOn;

  ElementSet copyToRecursive(const PiecewisePolynomial &other, Element *elt);
};
