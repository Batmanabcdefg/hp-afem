#pragma once

#include <fstream>
#include <iostream>
#include <map>

#include "elementfinder.h"
#include "element.h"
#include "triangleset.h"
#include "vertex.h"

class DOFHandler {
public:
  DOFHandler(ElementFinder &finder)
      : finder(finder), _valid(false), _numDOFs(0) {}

  ElementDofsSet asSet();
  const ElementDofsMap &map() const { return _g; }

  bool valid() { return _valid; }
  int recomputeNumDOFs();
  int maxDegree();

  void set(Element *elt, int dim);
  void reset(Element *elt, int dim);
  void construct(Element *elt, int dim);
  void erase(Element *elt);
  int detectDim(Element *elt);

  bool has(Element *elt);
  Dofs find(Element *elt);
  int find(Vertex *v);
  void determine(const ElementSet &elts, int dof);
  void determine(ElementDimsSet &eltdims);
  void redetermine();
  void redetermineOn(const ElementSet &elts);
  void transferToChildren(Element *parent);

  void increaseBy(Element *parent, int dof);
  void increaseTo(Element *parent, int dof);
  void increaseDegreeBy(Element *parent, int deg);
  void increaseBy(const ElementSet &elts, int dof);
  void increaseTo(const ElementSet &elts, int dof);
  void increaseDegreeBy(const ElementSet &elts, int deg);

  std::ostream &print(std::ostream &os = std::cout);

protected:
  ElementFinder &finder;
  ElementDofsMap _g;
  std::map<Vertex *, int> _vec;
  bool _valid;

  int _numDOFs;
  int increaseNumDOFs() { _numDOFs++; return _numDOFs-1; }
  void copyToChildren(Element *parent);
};
