#include "elementtree.h"

using namespace std;

ostream &ElementTree::printElements( ostream &os) {
  for(auto elt : _elts) {
    os << elt->index() << ": " << *elt << endl; 
  }
  return os;
}

void ElementTree::addVertex(Vertex *vert) {
  _verts.push_back( vert); 
  vert->_index = _verts.size()-1;
}

void ElementTree::addElement(Element *elt) { 
  // case: adding a root
  if( elt->isRoot()) {
    assert( (_roots.size() == _leaves.size()) && 
        "adding roots after bisection is prohibited");
    elt->_index = _elts.size();
  } else {
    elt->_index = _elts.size();
    /*
    // parent index
    long long pi = elt->parent()->index();

    // root the parent (and hence we) are in
    long long root = pi % _roots.size();   

    // parent index as if there were 1 root (single root parent index)
    long long srpi = (pi - root)/_roots.size();

    // single-root position of left child of parent is 2(srpi+1)-1
    elt->_index = _roots.size()*(2L*(srpi + 1L) - 1L) + root;

    // if we are the right child, add one for each root
    if( elt == elt->parent()->right()) elt->_index += _roots.size();
    */
  }

  _elts.push_back(elt);
}

Vertex *ElementTree::vertexAt(scalar x, scalar y) {
  for( auto &v : _verts) {
    if(v->x == x && v->y == y) {
      return v;
    }
  }
  return nullptr;
}

Element *ElementTree::element(long long i) {
  return _elts[i];
}

void ElementTree::resetRoots() {
  _roots.clear();
  for( auto &elt : _leaves) {
    addRoot(elt);
  }
}

ElementSet ElementTree::subtreeLeaves(Element *elt) {
  if(elt->isLeaf()) return ElementSet({elt});
  else {
    ElementSet v1 = subtreeLeaves(elt->left());
    ElementSet v2 = subtreeLeaves(elt->right());
    v1.insert(v2.begin(), v2.end());
    return v1;
  }
}

bool ElementTree::hasOverlap(ElementSet &s) const {
  ElementSet seen;

  for( auto &elt : s) {
    Element *cur = elt;
    if(seen.find(cur) != seen.end()) return true;
    while( !cur->isRoot()) {
      seen.insert(cur);
      cur = cur->parent();
    }
  }

  return false;
}

ElementTree::~ElementTree() {
  for( auto &v : _verts) delete v;
  for( auto &e : _elts)  delete e;
}
