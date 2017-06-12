#include "elementfinder.h"

using namespace std;

EltEdgePair ElementFinder::elementAlongEdge(Edge e) {
  VertPair vp = make_pair(e.v(0), e.v(1));
  return (_edges.count(vp) > 0) ? _edges[vp] : make_pair(nullptr, -1);
}

EltEdgePair ElementFinder::elementAlongEdge(Vertex *v0, Vertex *v1) {
  return elementAlongEdge(Edge(v0, v1));
}

EltEdgePair ElementFinder::elementOppositeEdge(Edge e) {
  return elementAlongEdge(e.v(1), e.v(0));
}

EltEdgePair ElementFinder::elementOppositeBisectionEdge(Element *elt) {
  return elementOppositeEdge(elt->bisectionEdge());
}

void ElementFinder::rebuildRecursively(Element *elt) {
  rebuild(elt);
  if( !elt->isLeaf()) {
    rebuildRecursively(elt->left());
    rebuildRecursively(elt->right());
  }
}

void ElementFinder::remove(Element *elt) {
  _edges.erase(make_pair(elt->v(0), elt->v(1)));
  _edges.erase(make_pair(elt->v(1), elt->v(2)));
  _edges.erase(make_pair(elt->v(2), elt->v(0)));
}

void ElementFinder::rebuild(Element *elt) {
  VertPair vp0 = make_pair(elt->v(0), elt->v(1)); _edges[vp0] = make_pair(elt, 0);
  VertPair vp1 = make_pair(elt->v(1), elt->v(2)); _edges[vp1] = make_pair(elt, 1);
  VertPair vp2 = make_pair(elt->v(2), elt->v(0)); _edges[vp2] = make_pair(elt, 2);
}

void ElementFinder::rebuild(ElementSet &_roots) {
  roots = _roots; rebuild();
}

void ElementFinder::rebuild() {
  for( auto root : roots) rebuildRecursively(root);
}
