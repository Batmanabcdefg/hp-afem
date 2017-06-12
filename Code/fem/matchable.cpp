#include "matchable.h"

using namespace std;

bool Matchable::isMatching( bool force) {
  if( !isConform(force)) {
    _isMatching = -1;
    if( !force) return false;
  }
  if( !force && _isMatching != 0) return _isMatching > 0;
  //for each leaf, we want the neighbour along its bisection edge to have its
  //bisection edge equal to ours
  for( auto elt : _leaves) {
    EltEdgePair nbrEdge = finder().elementOppositeBisectionEdge(elt);
    if( nbrEdge.first == nullptr) continue;
    if( !nbrEdge.first->isLeaf()) {
      _isMatching = -1;
      return false;
    }
    EltEdgePair nbrnbr = finder().elementOppositeBisectionEdge(nbrEdge.first);
    if( nbrnbr.first == nullptr || nbrnbr.first->index() != elt->index()) {
      _isMatching = -1;
      return false;
    }
  }
  _isMatching = 1;
  return true;
}

void Matchable::match( Element *elt) {
  //bisect
  bisect(elt);

  Element *l = elt->left();
  Element *r = elt->right();

  Vertex *vLeft[] = {l->v(0), l->v(1), l->v(2)};
  Vertex *vRight[] = {r->v(0), r->v(1), r->v(2)};

  l->_v[0] = vLeft[2]; l->_v[1] = vLeft[0]; l->_v[2] = vLeft[1];
  r->_v[0] = vRight[1]; r->_v[1] = vRight[2]; r->_v[2] = vRight[0];
}

void Matchable::makeMatching() {
  if(isMatching()) return;
  if( !isConform()) makeConform();

  //create subdivision
  ElementSet leaves(_leaves); //copy leaves to prevent infinite loops and such
  for( auto elt : leaves) {
    bisect(elt);
    match(elt->left());
    match(elt->right());
  }

  resetRoots();

  _isMatching = 1;

  cout << "made matching" << endl;
}

void Matchable::resetRoots() {
  Refinable::resetRoots();

  //rebuild edge matrix
  finder().rebuild(roots());

  //reset rhs
  rhs().retainOnly(roots());
}

bool Matchable::isTriTypesCorrect(bool force) {
  if( !isConform(force)) {
    _isTriTypesCorrect = -1;
    if( !force) return false;
  }

  if( !force && _isTriTypesCorrect != 0) return _isTriTypesCorrect > 0;
  for( auto elt : _leaves) {
    for( auto edge : elt->edges()) {
      EltEdgePair eltEdge = finder().elementAlongEdge(edge);
      EltEdgePair nbrEdge = finder().elementOppositeEdge(edge);
      Element *nbrElt = nbrEdge.first;
      if( nbrElt == nullptr) continue;

      assert(nbrElt->type()._isset); //not yet implemented
      if( nbrElt->type().i(nbrEdge.second) == elt->type().i(eltEdge.second)) return false;
    }
  }
  return true;
}

void Matchable::makeTriTypesCorrect() {
  if( isTriTypesCorrect()) return;
  resetRoots(); finder().rebuild(roots());
  assert(_roots.size() == _leaves.size()); //not yet implemented
  for( auto elt : _roots) {
    bool tt[3] = {false, false, false};
    for(auto edge : elt->edges()) {
      EltEdgePair eltEdge = finder().elementAlongEdge(edge);
      EltEdgePair nbrEdge = finder().elementOppositeEdge(edge);
      Element *nbrElt = nbrEdge.first;
      if( nbrElt == nullptr) continue;

      if( nbrElt->index() < elt->index()) { //we already assigned this neighbour a tritype
        assert(nbrElt->type()._isset); //not yet implemented
        TriType nbrTt = nbrElt->type();
        tt[eltEdge.second] = !nbrTt.i(nbrEdge.second);
      }
    }
    elt->_basis = &_bases.basis(TriType(tt[0], tt[1], tt[2]));
  }
  _isTriTypesCorrect = 1;
}
