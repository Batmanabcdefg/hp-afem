#pragma once
#include <utility>
#include <fstream>
#include <iostream>
#include <set>
#include <map>

#include "element.h"
#include "elementfinder.h"

class ElementTree {
protected:
  std::vector<Element *> _elts;
  std::vector<Vertex *> _verts;
  ElementSet _roots;
  ElementSet _leaves;

  ElementFinder _finder;

  void addRoot(Element *root) { _roots.insert(root); root->setRoot(); }
  void addLeaf(Element *leaf) { _leaves.insert(leaf); leaf->_isLeaf = true; }
  void addVertex(Vertex *vert);
  void addElement(Element *elt);
  void removeLeaf(Element *leaf) { _leaves.erase(_leaves.find(leaf)); leaf->_isLeaf = false; }
  Vertex *vertexAt(scalar x, scalar y);

public:
  std::ostream &printElements( std::ostream &os = std::cout);

  Element *element(long long i);
  int numElements() const { return _elts.size(); }

  ElementFinder &finder() { return _finder; }
  ElementSet &roots() { return _roots; }
  ElementSet &leaves() { return _leaves; }
  ElementSet subtreeLeaves(Element *elt);
  void rebuildElementFinder() { finder().rebuild(roots()); }

  virtual void resetRoots();

  bool hasOverlap(ElementSet &s) const;

  ~ElementTree();
};
