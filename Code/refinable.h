#pragma once

#include <vector>

#include "solvable.h"
#include "elementfinder.h"

class Refinable : public Solvable {
 public:
  void trim(Element *elt);

  const Bases &bases() { return _bases; }

  ElementSet refineElement(Element *elt, DOFHandler &handler);
  ElementSet refineElement(Element *elt) { return refineElement(elt, _handler); }

  ElementSet refine(ElementSet elts, DOFHandler &handler);
  ElementSet refine(ElementSet elts) { return refine(elts, _handler); }

  /**
   *  Creates a uniform refinement of the partition, so that we have four
   *  times as many leaves.
   */
  ElementSet refineLeavesUniformly(DOFHandler &handler);
  ElementSet refineLeavesUniformly() { return refineLeavesUniformly(_handler); }

  // determines if we are conforming.
  bool isConform(bool force = false);     

  // Construct the smallest conforming refinement from a (possibly)
  // nonconforming partition tree.
  void makeConform();                     
  void determineBoundaryVertices();
  std::set<Vertex *> vertsInside(Element *elt);
  void linearizeSolution();

  virtual void resetRoots() override;

 protected:
  // -1: no, 0: not initialized, 1: yes
  int _isConform = -1;

  Bases _bases;

  virtual ElementSet bisect(Element *elt, bool checkForBoundary = false);

  bool containsHangingVertex(Element *elt);

  Refinable(std::string basisdir) : _bases(std::move(basisdir)) {}
};
