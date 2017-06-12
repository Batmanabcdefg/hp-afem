#include <cassert>
#include <set>
#include <fstream>
#include <cmath>

#include "refinable.h"
#include "poly.h"

using namespace std;

/**
 *  Newest Vertex Bisection takes the newest vertex of the current element
 *  (denoted by the first vertex in its list of vertices), and draws a line from
 *  this newest vertex to the middle of the opposing edge. The two triangles
 *  that arise from this have their newest vertex set to the vertex that was
 *  just created.
 *
 *  In our case, after bisecting, we update: the set of leaves (remove current 
 *  element, add its children); the corresponding parts of the edge matrix for
 *  quickly finding neighbours; the local part of the right-hand side of the 
 *  PDE; if the solution is a polynomial locally on this element, we copy it
 *  to the children. Finally, we set isConform to 0, which means "we don't know".
 *
 *  @param  elt               the element to bisect
 *  @param  checkForBoundary  whether to check the new vertex for being on bdr
 *  @return ElementSet        the set of newly created leaves.
 */
ElementSet Refinable::bisect( Element *elt, bool checkForBoundary) {
  // this of course only works on leaves
  // force it to be a _real_ leaf
  if(elt->isLeaf(true)) {
    EltEdgePair nbredge = finder().elementOppositeBisectionEdge(elt);
    Vertex *newvert = nullptr;

    // find out if the new vertex exists already
    if(nbredge.first != nullptr) {
      if(!nbredge.first->isLeaf()) {
        // we are opposite the nbr's bisection edge
        if( nbredge.second == 1) {
          // copy vertex
          newvert = nbredge.first->left()->v(0);
        } else {
          if(nbredge.second == 0 && !nbredge.first->left()->isLeaf()) {
            //copy vertex
            newvert = nbredge.first->left()->left()->v(0);
          } else if(nbredge.second == 2 && !nbredge.first->right()->isLeaf()) {
            //copy vertex
            newvert = nbredge.first->right()->left()->v(0);
          }
        }
      }
    }

    // TODO: it is also possible that the vertex _does_ exist already, but the triangle was trimmed; you cannot find it using the element finder any more, but it is
    // still there: look one or two generations earlier.
    if( newvert == nullptr) {
      //new vertex
      //find newest vertex
      Edge edge = elt->bisectionEdge();
      Vertex *v0 = edge.v(0), *v1 = edge.v(1);

      // make sure we sum the two in the same order for both sides of the same edge
      if( v0->index() > v1->index()) {
        swap(v0, v1);
      }

      scalar x = (v0->x + v1->x)/2;
      scalar y = (v0->y + v1->y)/2;

      //add newest vertex
      newvert = vertexAt(x,y);
      if(newvert == nullptr) {
        newvert = new Vertex(x,y);

        //boundary check
        if( checkForBoundary && nbredge.first == nullptr) {
          newvert->_isBoundary = true;
        }

        addVertex(newvert);
      }
    } else {
      //cout << "copied vertex!" << endl;
    }

    // create new elements: correct vertices, tritype
    Vertex *lVerts[] = {_verts[newvert->index()], _verts[elt->i(0)], _verts[elt->i(1)]},
           *rVerts[] = {_verts[newvert->index()], _verts[elt->i(2)], _verts[elt->i(0)]};
    Element *lElt = new Element(lVerts, elt, nullptr, nullptr);
    Element *rElt = new Element(rVerts, elt, nullptr, nullptr);
    if( elt->type()._isset) {
      lElt->_basis = &_bases.basis(elt->type().left());
      rElt->_basis = &_bases.basis(elt->type().right());
    }
    elt->_left = lElt;
    elt->_right = rElt;

    // add new elements to the partition
    addElement(lElt);
    addElement(rElt);
  }

  // update leaves
  removeLeaf(elt);
  addLeaf(elt->left());
  addLeaf(elt->right());

  // update edge matrix
  finder().rebuild(elt->left());
  finder().rebuild(elt->right());

  // possibly update solution
  if( _sol.has(elt)) _sol.copyToChildren(elt);

  // update right hand side
  _rhs.copyToChildren(elt);

  //we don't know if we are still conform
  _isConform = 0;               

  // return the two elements we just created
  return ElementSet({elt->left(), elt->right()});
}

/**
 *  Trims the given element, removing all its children from the tree. This makes
 *  a recursive call.
 *
 *  This messes up conformity for certain; correctness is probably still good. 
 *  As for matchingness, I honesly don't know at this point. TODO. Not sure if
 *  we should care, either.. ElementFinder edges are messed up as well, but we 
 *  might be able to patch those while trimming.
 *
 *  @param  elt   the element
 */
void Refinable::trim(Element *elt) {
  // this only works if we are no leaf
  assert(!elt->isLeaf());
  Element *children[] = {elt->left(), elt->right()};
  //elt->_left = nullptr;
  //elt->_right = nullptr;

  for( auto &child : children) {
    // recursively trim the subtree first
    if( !child->isLeaf()) trim(child);

    // remove this element from the set of leaves
    removeLeaf(child);

    // remove this element from the piecewise polynomials rhs and sol
    _rhs.erase(child);
    if( _sol.has(child)) _sol.erase(child);

    // remove it from the DoFHandler
    if( handler().has(child)) handler().erase(child);

    // update the ElementFinder
    finder().remove(child);
  }
  finder().rebuild(elt);

  // this element is now a leaf
  addLeaf(elt);

  //we don't know if we're still conforming
  _isConform = 0;
}

ElementSet Refinable::refineElement(Element *elt, DOFHandler &handler) {
  // this only works if we're a leaf
  assert(elt->isLeaf());

  // the resulting ElementSet
  ElementSet res;

  // find neighbour
  Edge edge = elt->bisectionEdge();
  EltEdgePair nbrEdge = finder().elementOppositeEdge(edge);
  Element *nbr = nbrEdge.first;

  // case one: we are on boundary
  if( nbr == nullptr) {
    // fill res
    TriangleSet::unionSetsInto(res, bisect(elt, true));

    // transfer DOFs
    handler.transferToChildren(elt);
    return res;
  }

  // this neighbour is non-null; find its bisection edge
  Edge nbrBisectEdge = nbr->bisectionEdge();

  // case two: our bisection edge is their bisection edge
  if( edge.isOpposite(nbrBisectEdge)) {
    // bisect both this element and its neighbour
    TriangleSet::unionSetsInto(res, bisect(elt, true));
    TriangleSet::unionSetsInto(res, bisect(nbr, true));

    // tell the handler to allocate a DOF vector for both children of the neighbour;
    // this is necessary as the subsequent call to transferToChildren(elt) will
    // fill values there.
    handler.construct(nbr->left(), handler.detectDim(nbr));
    handler.construct(nbr->right(), handler.detectDim(nbr));

    // transfer the DOFs
    handler.transferToChildren(elt);
    handler.transferToChildren(nbr);
    return res;
  }

  // case three: harder. recursively call refineElement on the neighbour first
  res = refineElement(nbr, handler);
  // then bisect ourself
  TriangleSet::unionSetsInto(res, bisect(elt, true));

  // is our bisection edge opposite the neighbours left child or no?
  Edge leftEdge = nbr->left()->bisectionEdge();
  Element *nbrchild = edge.isOpposite(leftEdge) ? nbr->left() : nbr->right();

  // bisect this child as well
  TriangleSet::unionSetsInto(res, bisect(nbrchild, true));

  // tell the handler to allocate a DOF vector
  handler.construct(nbrchild->left(), handler.detectDim(nbrchild));
  handler.construct(nbrchild->right(), handler.detectDim(nbrchild));

  // transfer the DOFs
  handler.transferToChildren(elt);
  handler.transferToChildren(nbrchild);

  return res;
}

ElementSet Refinable::refine(ElementSet elts, DOFHandler &handler) {
  // this only works if we are conform
  assert(isConform());

  // the resulting ElementSet
  ElementSet res;

  // refine each element in the set. As it's possible that a call to refineElement
  // already refines some other element in this set, first we check if it is
  // still a leaf.
  for( auto &elt : elts) {
    if( elt->isLeaf()) {
      // refine and save the resulting elements into res
      TriangleSet::unionSetsInto(res, refineElement(elt, handler));
    }
  }

  cout << "ik ben hier" << endl;

  // clean-up phase; remove those elements that've become non-leaves a 
  // result of this algorithm
  // TODO: try and do this a little nicer
  ElementSet ret;
  for( auto &elt : res) {
    if( elt->isLeaf()) ret.insert(elt);
  }

  // we know we're still conforming
  _isConform = 1;

  return ret;
}

/**
 *  Bisect each leaf twice, yielding your leaves per original leaf.
 */
ElementSet Refinable::refineLeavesUniformly(DOFHandler &handler) {
  ElementSet ccurleaves = leaves();
  cout << "first refine" << endl;
  refine(ccurleaves, handler);
  ElementSet newleaves;
  cout << "nu hier" << endl;
  for( auto &elt : ccurleaves) {
    newleaves.insert(elt->left());
    newleaves.insert(elt->right());
  }
  cout << "second refine" << endl;
  return refine(newleaves, handler);
}

set<Vertex *> Refinable::vertsInside(Element *elt) {
  if(elt->isLeaf()) {
    return set<Vertex *>({elt->v(0), elt->v(1), elt->v(2)});
  } else {
    set<Vertex *> s1 = vertsInside(elt->left());
    set<Vertex *> s2 = vertsInside(elt->right());
    s1.insert(s2.begin(), s2.end());
    return s1;
  }
}

void Refinable::linearizeSolution() {
  while(true) { 
    ElementSet toRefine;
    for( auto &elt : _sol.definedOn()) {
      //auto size = vertsInside(elt).size();
      ElementSet eltleaves = subtreeLeaves(elt);
      scalar factor = pow(0.5L, elt->gen());
      if(eltleaves.size() < factor * _sol.locallyAt(elt).size()) {
        toRefine.insert(eltleaves.begin(), eltleaves.end());
      }
    }
    if(toRefine.size()) {
      refine(toRefine);
    } else {
      break;
    }
  }
}

/**
 *  This method finds if an element contains a hanging vertex, i.e., if some
 *  vertex from another element is inside this one. Very crappy drawing where
 *  the left element contains a hanging vertex:
 *   /|\
 *  < |->
 *   \|/
 *  I believe this is also the _only_ way we can have a hanging vertex.
 *    
 *  Old way -- O(N):
 *    for each leaf:
 *      for each edge of current element:
 *        for each vertex of leaf:
 *          if our edge contains vertex:
 *            continue;
 *          else:
 *            if vertex is on edge:
 *              return true
 *    return false;
 *  New way -- O(1):
 *    for each edge of element:
 *      find neighbour along edge;
 *      if nbr == nullptr:
 *        we are on the edge of the domain, so continue;
 *      if nbr is a leaf:
 *        we share this edge with a leaf: no hanging vertex possible, so continue;
 *      if edge is opposite the bisection edge of nbr:
 *        nbr has children with newest vertex inside element, so return true;
 *    return false;
 */
bool Refinable::containsHangingVertex(Element *elt) {
  assert(elt->isLeaf());

  for( auto edge : elt->edges()) {
    // get the neighbour opposite the current edge
    EltEdgePair nbrEdge = finder().elementOppositeEdge(edge);
    Element *nbr = nbrEdge.first;

    // are we on the edge of the domain?
    if( nbr == nullptr) continue; 

    // not possible to have a hanging vertex here
    if( nbr->isLeaf()) continue;

    // see above
    if( edge.isOpposite(nbr->bisectionEdge())) return true;
  }

  return false;
}

/**
 *  Check if we are conforming.  If force is false, look it up in the private
 *  member first; if it is true, always do a "full check".
 *
 *  Full check is done by just looping over all leaves, looking if they have a
 *  hanging vertex.
 */
bool Refinable::isConform( bool force) {
  // if no force, just look it up
  if( !force && _isConform != 0) return _isConform > 0;

  // loop over all leaves
  for( auto &elt : _leaves) {
    if( containsHangingVertex(elt)) {
      _isConform = -1;
      return false;
    }
  }

  // no leaf has a hanging vertex
  _isConform = 1;
  return true;
}

/**
 *  See "Optimality of a standard adaptive finite element method" by R.P.
 *  Stevenson, page 5, for this algorithm. In pseudocode:
 *  let M be an empty set;
 *  for each node in the partition:
 *    if node contains a hanging vertex, put node in M;
 *
 *  while M is not empty:
 *    extract node from M;
 *    bisect node;
 *    for both children of node:
 *      if child contains a hanging vertex, put child into M;
 *    for each neighbour of node:
 *      if nbr contains a hanging vertex:
*         put nbr in M;
 */
void Refinable::makeConform() {
  if(isConform()) return;

  ElementSet M;
  for( auto &elt : _leaves) if( containsHangingVertex(elt)) M.insert(elt);

  while(M.size()) {
    auto it = M.begin();
    Element *elt = *it;
    M.erase(it);

    bisect(elt);

    if( containsHangingVertex(elt->left())) {
      M.insert(elt->left());
    }
    if( containsHangingVertex(elt->right())) {
      M.insert(elt->right());
    }

    for( auto edge : elt->edges()) {
      EltEdgePair nbrEdge = finder().elementOppositeEdge(edge);
      Element *nbr = nbrEdge.first;
      if( nbr == nullptr) continue;
      if( !nbr->isLeaf()) continue;
      if( containsHangingVertex(nbr)) {
        M.insert(nbr);
      }
    }
  }
  _isConform = 1;
}

/**
 *  A vertex v is a boundary vertex iff there is a vertex w such that
 *  along (v, w), there is an element, but along (w, v), there is none.
 *
 *  So: for each vertex, loop over all other vertices.  If the above condition
 *  holds, we can break.
 */
void Refinable::determineBoundaryVertices() {
  makeConform();
  for( auto verti : _verts ) {
    verti->_isBoundary = false;
    for( auto vertj : _verts ) {
      if( finder().elementAlongEdge(verti, vertj).first != nullptr && 
          finder().elementAlongEdge(vertj, verti).first == nullptr) {
        verti->_isBoundary = true; 
        break;
      }
    }
  }
}

void Refinable::resetRoots() {
  ElementTree::resetRoots();
  cerr << "RESET ROOTS"  << endl;
  for(auto &elt : _roots) {
    elt->_eltmats = new ElementMatrices(_bases, elt);
  }
}
