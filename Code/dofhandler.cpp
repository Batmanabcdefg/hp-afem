#include "dofhandler.h"

#include <cassert>
#include <set>

#include "degree.h"
#include "print.h"

using namespace std;

ElementDofsSet DOFHandler::asSet() {
  ElementDofsSet current;
  for (const auto& ed : _g) {
    current.insert(ed);
  }
  return current;
}

int DOFHandler::recomputeNumDOFs() {
  std::set<int> vals;
  int maxDof = -1;
  for (const auto& ed : _g) {
    cout << *(ed.first);
    ed.second.print();
    for (int gi : ed.second.values()) {
      if (gi != -1) {
        maxDof = max(maxDof, gi);
        vals.insert(gi);
      }
    }
  }
  //for( int dof : vals) cout << dof << " "; cout << endl;
  //cout << maxDof + 1 << " " << vals.size() << endl;
  assert(vals.size() == maxDof+1);

  return 1 + maxDof;
}

void DOFHandler::redetermine() {
  ElementDimsSet current;
  for (const auto& ed : _g) {
    current.insert(make_pair(ed.first, ed.second.dof()));
  }
  determine(current);
}

int DOFHandler::maxDegree() {
  int max_degree = -1;
  for (const auto &ed : _g) {
    if (max_degree < ed.second.degree()) {
      max_degree = ed.second.degree();
    }
  }
  return max_degree;
}

int DOFHandler::detectDim(Element *elt) {
  if (_g.count(elt) > 0) {
    return find(elt).size();
  }
  assert(!elt->isRoot());
  return detectDim(elt->parent());
}

void DOFHandler::redetermineOn(const ElementSet &elts) {
  ElementDimsSet cur;
  for (const auto &e : elts) {
    cur.insert(make_pair(e, detectDim(e)));
  }
  determine(cur);
}

void DOFHandler::determine(const ElementSet &elts, int dof) {
  ElementDimsSet current;
  for (const auto &elt : elts) {
    current.insert(make_pair(elt, dof));
  }
  determine(current);
}

void DOFHandler::determine(ElementDimsSet &eltdims) {
  // reset everything
  _numDOFs = 0;
  _vec.clear();
  _g.clear();

  // reset the dof vectors
  for (const auto &ed : eltdims) {
    Element *elt = ed.first;
    int dim = Degree::dofToDim(ed.second);
    _g.insert(make_pair(elt, Dofs(dim)));
  }

  for (const auto &eg : _g) {
    Element *elt = eg.first;

    // This elements DOF vector, possibly partially filled already.
    Dofs g = eg.second;

    for (int i = 0; i < 3; i++) {
      // If we haven't seen this vertex before, add it to our vertex map.
      if (_vec.count(elt->v(i)) == 0) {
        _vec[elt->v(i)] = -1;
      }

      // If it is on the boundary, it should not get a DOF.
      if (elt->v(i)->isBoundary()) {
        continue;
      }

      // In all other cases, it should.
      if (_vec[elt->v(i)] == -1) {
        _vec[elt->v(i)] = increaseNumDOFs();
      }

      // Propagate this DOF change to this element.
      g.setVertex(i, _vec[elt->v(i)]);
    }

    if (g.degree() >= 2) {
      for (int k = 1; k <= g.degree(); k++) {
        // face nodes
        if (k >= 3) {
          for (int r1 = 0; r1 <= k-3; r1++) {
            g.setFace(k, r1, increaseNumDOFs());
          }
        }

        // edge nodes
        if (k < g.degree()) {
          for (const auto& edge : elt->edges()) {
            EltEdgePair leafEdge = finder.elementAlongEdge(edge);
            EltEdgePair nbrEdge = finder.elementOppositeEdge(edge);

            //on the bdr
            if (nbrEdge.first == nullptr) {
              continue;
            }

            Element *nbr = nbrEdge.first;
            Dofs nbrg = find(nbr);
            if (!(k < nbrg.degree())) {
              continue;
            }

            // found opposing edge; either create or copy DOF
            if (nbrg.getEdge(k, nbrEdge.second) == -1) {
              g.setEdge(k, leafEdge.second, increaseNumDOFs());
            } else { // We share a global index with this guy.
              g.setEdge(k, leafEdge.second, nbrg.getEdge(k, nbrEdge.second));
            }

            _g[nbr] = nbrg;
          }
        }
      }
    }

    _g[elt] = g;
  }

  // The values in this DOFHandler are valid.
  _valid = true;
}

void DOFHandler::set(Element *elt, int dim) {
  Dofs g;
  g.reset(dim);
  _g[elt] = g;
  _valid = false;
}

void DOFHandler::reset(Element *elt, int dim) {
  Dofs g = find(elt);
  g.reset(dim); 
  _g[elt] = g;
  _valid = false;
}

void DOFHandler::construct(Element *elt, int dim) {
  assert(_g.count(elt) == 0);
  Dofs g;
  g.reset(dim);
  _g[elt] = g;
}

void DOFHandler::erase(Element *elt) {
  auto it = _g.find(elt);
  assert(it != _g.end());
  _g.erase(it);
  _valid = 0;
}

void DOFHandler::transferToChildren(Element *parent) {
  copyToChildren(parent);
  _g.erase(parent);
}

void DOFHandler::copyToChildren(Element *parent) {
  assert(valid());
  assert(!parent->isLeaf());
  Dofs g = find(parent);

  Element *left = parent->left();
  Element *right = parent->right();

  // Reset all.
  Dofs leftg, rightg;
  leftg.reset(g.dim());
  rightg.reset(g.dim());

  // If newest vertex is not on boundary, assign it a DOF.
  if (_vec.find(left->v(0)) == _vec.end()) {
    _vec[left->v(0)] = -1;
  }
  if (!left->v(0)->isBoundary()) {
    if (_vec[left->v(0)] == -1) {
      _vec[left->v(0)] = increaseNumDOFs();
    }
    leftg.setVertex(0, _vec[left->v(0)]);
    rightg.setVertex(0, _vec[left->v(0)]);
  }

  // Copy vertex DOFs from parent.
  leftg.setVertex(1, _vec[left->v(1)]); 
  rightg.setVertex(1, _vec[right->v(1)]);
  leftg.setVertex(2, _vec[left->v(2)]);
  rightg.setVertex(2, _vec[right->v(2)]);

  // Find neighbour along bisection edge.
  auto nbrEdge = finder.elementOppositeBisectionEdge(parent);
  Element *nbr = nbrEdge.first;

  // they must share a refinement edge, and it needs to have children as well.
  assert(nbr == nullptr || (nbr->v(1) == parent->v(2) 
                         && nbr->v(2) == parent->v(1) 
                         && !nbr->isLeaf()));

  if (g.degree() >= 2) {
    for (int k = 1; k <= g.degree(); k++) {
      // Add face DOFs; these are easy.
      if (k >= 3) {
        for (int r1 = 0; r1 <= k-3; r1++) {
          leftg.setFace(k, r1, g.getFace(k, r1));
          rightg.setFace(k, r1, increaseNumDOFs());
        }
      }

      if (k < g.degree()) {
        // Copy parent v0-v1 edge DOFs to left v1-v2 edge DOFs.
        leftg.setEdge(k, 1, g.getEdge(k, 0));
        // Copy parent v2-v0 edge DOFs to right v1-v2 edge DOFs.
        rightg.setEdge(k, 1, g.getEdge(k, 2));

        // Give the shared edge a new DOF.
        leftg.setEdge(k, 0, increaseNumDOFs());
        rightg.setEdge(k, 2, leftg.getEdge(k, 0));

        // Bisection edge of parent was on boundary, so no DOFs for children 
        // here.
        if( nbr == nullptr) {
          continue;
        }

        // Get the neighbours and its children's DOF vectors.
        Dofs nbrleftg = find(nbr->left());
        Dofs nbrrightg = find(nbr->right());

        // Have a DOF along this edge?
        // TODO: why > k and not >= k?
        if (nbrrightg.degree() > k) {
          // If this is -1, our right neighbour will get its DOFs copied in a
          // bit.
          if (nbrrightg.getEdge(k, 0) == -1) {
            // So just copy from parent to here.
            leftg.setEdge(k, 2, g.getEdge(k, 1));
          } else {
            // Else, we copy.
            leftg.setEdge(k, 2, nbrrightg.getEdge(k, 0));
          }
        }

        // Have a DOF along this edge?
        if (nbrleftg.degree() > k) {
          // If this is -1, our right neighbour will get its DOFs copied in a
          // bit.
          if (nbrleftg.getEdge(k, 2) == -1) {
            // So just copy create a new DOF.
            rightg.setEdge(k, 0, increaseNumDOFs());
          } else {
            // Else, we copy.
            rightg.setEdge(k, 0, nbrleftg.getEdge(k, 2));
          }
        }

        _g[nbr->left()] = nbrleftg;
        _g[nbr->right()] = nbrrightg;

        // TODO: why is this commented out?
        /*
        if( nbr->index() > parent->index()) {
          if( nbrrightg.degree() >= k) leftg.setEdge(k, 2, g.getEdge(k, 1));
          if( nbrleftg.degree() >= k) rightg.setEdge(k, 0, increaseNumDOFs());
        } else {
          // our left v2-v0 bdr is their right v0-v1 bdr and vice versa
          if( nbrrightg.degree() >= k) leftg.setEdge(k, 2, nbrrightg.getEdge(k, 0));
          if( nbrleftg.degree() >= k) rightg.setEdge(k, 0, nbrleftg.getEdge(k, 2));
        }
        */
      }
    }
  }

  _g[left] = leftg;
  _g[right] = rightg;

  _valid = true;
}

void DOFHandler::increaseBy(const ElementSet &elts, int dof) {
  for (auto &elt : elts) {
    increaseBy(elt, dof);
  }
}

void DOFHandler::increaseDegreeBy(Element *elt, int deg) {
  int residue = find(elt).dof() - find(elt).dim();
  increaseBy(elt, Degree::degreeToDim(find(elt).degree() + deg) - Degree::degreeToDim(find(elt).degree()) + residue);
}

void DOFHandler::increaseTo(Element *elt, int dof) {
  assert(find(elt).dof() <= dof);
  return increaseBy(elt, dof - find(elt).dof());
}

void DOFHandler::increaseDegreeBy(const ElementSet &elts, int deg) {
  for (auto &elt : elts) {
    increaseDegreeBy(elt, deg);
  }
}

void DOFHandler::increaseTo(const ElementSet &elts, int dof) {
  for (auto &elt : elts) {
    increaseTo(elt, dof);
  }
}

void DOFHandler::increaseBy(Element *parent, int dof) {
  cout << parent->index() << " " << dof << endl;
  assert(valid());
  Dofs g = find(parent);

  int stopdim = min(g.dof() + dof, parent->basis()->dim());

  if(Degree::dofToDim(g.dof() + dof) > parent->basis()->dim()) {
    cout << "Truncating degree!!" << endl;
  }
  
  int startdegree = Degree::dofToDegree(g.dof())+1;
  int stopdegree = Degree::dofToDegree(stopdim);
  g.resize(stopdim);

  if (stopdegree >= 2) {
    for (int k = startdegree; k <= stopdegree; k++) {
      // Face DOFs are easy.
      if (k >= 3) {
        for (int r1 = 0; r1 <= k-3; r1++)  {
          g.setFace(k, r1, increaseNumDOFs());
        }
      }

      // We require k-1 instead of k here because the edge DOFs will otherwise
      // not be given properly.
      // TODO: what the hell man
      if (k-1 < stopdegree) {
        for (auto edge : parent->edges()) {
          EltEdgePair leafEdge = finder.elementAlongEdge(edge);
          EltEdgePair nbrEdge = finder.elementOppositeEdge(edge);

          // On the bdr.
          if (nbrEdge.first == nullptr) {
            continue;
          }
          Element *nbr = nbrEdge.first;
          Dofs nbrg = find(nbr);
          if (!(k-1 < nbrg.degree())) {
            continue;
          }

          assert(nbrg.getEdge(k-1, nbrEdge.second) == -1);
          g.setEdge(k-1, leafEdge.second, increaseNumDOFs());
          nbrg.setEdge(k-1, nbrEdge.second, g.getEdge(k-1, leafEdge.second));

          // Update the DOFs again.
          _g[nbr] = nbrg;
        }
      }
    }
  }

  _g[parent] = g;
}

bool DOFHandler::has(Element *elt) {
  return _g.find(elt) != _g.end();
}

Dofs DOFHandler::find(Element *elt) {
  auto it = _g.find(elt);
  assert(it != _g.end());
  return it->second;
}

int DOFHandler::find(Vertex *v) {
  auto it = _vec.find(v);
  if(it != _vec.end()) return it->second;
  return -1;
}

ostream &DOFHandler::print( ostream &os) {
  os << Print::formatted("%lu", _g.size()) << endl;
  for (const auto &ed : _g) {
    Element *elt = ed.first;
    Dofs g = ed.second;
    os << Print::formatted("%lu %d %d %d %d", g.dim(), elt->i(0), elt->i(1),
                           elt->i(2), elt->type().toInt());
    for (const auto gi : g.values()) {
      os << " " << gi;
    }
    os << endl;
  }

  return os;
}
