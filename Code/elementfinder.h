#pragma once

#include <map>
#include <utility>

#include "triangleset.h"
#include "vertex.h"
#include "element.h"
#include "edge.h"

typedef struct std::pair<Vertex *,Vertex *> VertPair;
typedef struct std::pair<Element *,int> EltEdgePair;

class ElementFinder {
  ElementSet roots;
  std::map<VertPair, EltEdgePair> _edges;
  void rebuildRecursively(Element *elt);

public:
  ElementFinder() {}
  ElementFinder(ElementSet &roots) : roots(roots) { rebuild(); }

  EltEdgePair elementAlongEdge(Edge e);
  EltEdgePair elementAlongEdge(Vertex *v0, Vertex *v1);
  EltEdgePair elementOppositeEdge(Edge e);
  EltEdgePair elementOppositeBisectionEdge(Element *elt);

  void remove(Element *elt);
  void rebuild(Element *elt);
  void rebuild(ElementSet &roots);
  void rebuild();
};
