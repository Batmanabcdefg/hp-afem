#include <math.h>
#include <cfloat>
#include <queue>
#include <fstream>

#include "print.h"
#include "nearbest.h"
#include "approximator.h"
#include "degree.h"
#include "math.h"

#define DEBUG 0

using namespace std;

enum class ElementType { Combined, Element, Virtual};

void NearBest::combineRoots() {
  queue<Element *> push;
  queue<Element *> pull;

  for( auto &root : partition.roots()) {
    push.push(root);
  }

  int totalAdded = 0;

  while(push.size() >= 2) {
    //swap push and pull
    swap(push, pull);
    //while pull has two elements:
    while( pull.size() >= 2) {
      //take two out, combine, push it
      Element *one = pull.front(); pull.pop();
      Element *two = pull.front(); pull.pop();
      //this is required to make the elements still sortable
      Element *newElt = new Element(one, two, -(++totalAdded));
      if (DEBUG) {
        cout << "inserting a combination " << newElt->index()
             << " as parent of " << one->index() << " and " << two->index()
             << " " << totalAdded << endl;
      }
      _combinations.insert(newElt);
      one->_parent = newElt;
      two->_parent = newElt;

      push.push(newElt);
    }

    //if pull has element, create a virtual
    if( pull.size()) {
      //create virtual
      //this is required to make the elements still sortable
      Element *virt = new Element(nullptr, nullptr, -(++totalAdded));
      _virtuals.insert(virt);

      //take out last
      Element *last = pull.front(); pull.pop();

      //combine, push it
      //this is required to make the elements still sortable
      Element *newElt = new Element(last, virt, -(++totalAdded));
      if (DEBUG) {
        cout << "inserting virtual under " << newElt->index() << endl;
      }
      _combinations.insert(newElt);
      last->_parent = newElt;
      virt->_parent = newElt;
      push.push(newElt);
    }
  }
  assert(push.size() == 1);
  _root = push.front(); push.pop();
}

void NearBest::removeCombinedRoots() {
  removeCombinedRootsRecursive(_root);
  _virtuals.clear();
  _combinations.clear();
  _root = nullptr;
}

void NearBest::removeCombinedRootsRecursive(Element *elt) {
  if (DEBUG) {
    cout << "now looking at " << elt->index() << endl;
  }
  // virtuals are leaves
  if( _virtuals.find(elt) != _virtuals.end()) {
    //cout << "erasing virtual " << elt->index() << endl;
    if(_nearbestleaves.find(elt) != _nearbestleaves.end()) 
      _nearbestleaves.erase(elt);
    _virtuals.erase(elt);
    assert(elt->isLeaf());

    // remove every trace
    _e.erase(elt);
    _ehp.erase(elt);
    _ehpTilde.erase(elt);
    _eTilde.erase(elt);
    _r.erase(elt);
    _t.erase(elt);
    _q.erase(elt);
    delete elt;
    return;
  }

  if( _combinations.find(elt) != _combinations.end()) {
    if (DEBUG) {
      cout << "erasing combination " << elt->index() << " with children "
           << elt->left()->index() << " " << elt->right()->index() << endl;
    }
    assert(!elt->isLeaf());

    // if we are a nearbest leaf, make its children nearbest leaves.
    // recursive call will take care of possible hiccups
    if(_nearbestleaves.find(elt) != _nearbestleaves.end()) {
      if (DEBUG) {
        cout << "we are a nearbest leaf; inserting to children" << endl;
      }
      assert(!elt->isLeaf());
      _nearbestleaves.insert(elt->left());
      _nearbestleaves.insert(elt->right());
      _nearbestleaves.erase(elt);

      //special case: a combination node is a nearbest leaf; if its left- or right
      //child is not a combination nor a virtual, we should trim that
      for(auto &child : {elt->left(), elt->right()}) {
        if(_combinations.find(child) == _combinations.end() &&
           _virtuals.find(child) == _virtuals.end()) {
          if (DEBUG) {
            cout << "special case!" << endl;
          }
          trimNode(child);
        }
      }
    }
    _combinations.erase(elt);

    elt->left()->_parent = nullptr;
    elt->right()->_parent = nullptr;
    if (DEBUG) {
      cout << "recursing into children" << endl;
    }
    removeCombinedRootsRecursive(elt->left());
    removeCombinedRootsRecursive(elt->right());

    if (DEBUG) {
      cout << "children of " << elt->index() << " are fixed" << endl;
    }
    // remove every trace
    _e.erase(elt);
    _ehp.erase(elt);
    _ehpTilde.erase(elt);
    _eTilde.erase(elt);
    _r.erase(elt);
    _t.erase(elt);
    _q.erase(elt);

    delete elt;
    return;
  }
}

scalar NearBest::error(Element *elt, int dof) {
  int degree = Degree::dofToDegree(dof);
  int dim = Degree::dofToDim(dof);

  if( !_e[elt].count(degree)) {
    if( _virtuals.find(elt) != _virtuals.end()) {
      _e[elt][degree] = 0.0;
    } else if( _combinations.find(elt) != _combinations.end()) {
      scalar minval = std::numeric_limits<scalar>::max();
      for( int d1 = 0; d1 <= dim; d1++) {
        minval = min(minval, error(elt->left(), d1) + error(elt->right(), dim-d1));
      }
      _e[elt][degree] = minval;
    } else {
      _e[elt][degree] = Approximator(poly).error(elt, dof);
    }
  }

  //cout << "error of " << elt->index() << " with " << dof << " dof (deg " << degree << ") is " << _e[elt][degree] << endl;
  return _e[elt][degree];
}

scalar NearBest::ehp(Element *elt, int dof, scalar val) {
  int degree = Degree::dofToDegree(dof);
  if(val != -1) {
    _ehp[elt][degree] = val;
  } else {
    assert(_ehp[elt].count(degree));
  }
  return _ehp[elt][degree];
}

scalar NearBest::ehpTilde(Element *elt, int dof, scalar val) {
  int degree = Degree::dofToDegree(dof);
  if(val != -1) {
    _ehpTilde[elt][degree] = val;
  } else {
    assert(_ehpTilde[elt].count(degree));
  }
  return _ehpTilde[elt][degree];
}

scalar NearBest::eTilde(Element *elt, scalar val) {
  if( val != -1) {
    _eTilde[elt] = val;
  } else {
    assert(_eTilde.count(elt));
  }
  return _eTilde[elt];
}

int NearBest::r_forced(Element *elt, int val) {
  if( val != -1) return _r[elt] = val;

  // check if we have a correctly set r(elt); if not, create it
  if( !_r.count(elt) || (_r[elt] == -1)) {
    // check if we are a leaf of the near best tree
    if( _nearbestleaves.find(elt) != _nearbestleaves.end()) {
      return _r[elt] = 1;
    } else {
      // to die if we are no leaf and none of our ancestors isnt either
      if( elt->isLeaf()) {
        cout << *elt << endl;
        assert(false);
      }
      return _r[elt] = r_forced(elt->left()) + r_forced(elt->right());
    }
  }

  // return the value
  return _r[elt];
}

int NearBest::r(Element *elt, int val) {
  return min(partition.bases().dim(), r_forced(elt, val));
}

Element *NearBest::t(Element *elt, Element *val) {
  if( val != nullptr) {
    _t[elt] = val;
  } else {
    assert(_t.count(elt));
  }
  // return the value
  return _t[elt];
}

scalar NearBest::q(Element *elt, scalar val) {
  if( val != -1) {
    _q[elt] = val;
  } else {
    assert(_q.count(elt));
  }
  return _q[elt];
}

void NearBest::eraseFromDataStructures(Element *elt) {
  // if we are _no_ nearbest leaf, recurse until we are
  if( _nearbestleaves.find(elt) == _nearbestleaves.end()) {
    if(!elt->isLeaf()) {
      eraseFromDataStructures(elt->left());
      eraseFromDataStructures(elt->right());
    } else {
      //this can happen exactly when a combination node is a nearbest leaf, and 
      //we call the trim operation on a root that is a child of this combination.
      //a recursive trimming operation then comes here
    }
    // now we are
    _nearbestleaves.insert(elt);
  }

  // and remove us
  _nearbestleaves.erase(elt);
  _e.erase(elt);
  _ehp.erase(elt);
  _ehpTilde.erase(elt);
  _eTilde.erase(elt);
  _r.erase(elt);
  _t.erase(elt);
  _q.erase(elt);
}

void NearBest::trimNode(Element *elt) {
  if(!elt->isLeaf()) {
    eraseFromDataStructures(elt->right());
    eraseFromDataStructures(elt->left());
    partition.trim(elt);
  }
  _nearbestleaves.insert(elt);
}

/**
 *  Starting from the leaves of the nearbest tree, walk upwards to the roots,
 *  trimming whenever it is better to use hp than h.
 *
 *  @param  elt   the element
 */
void NearBest::trimRecursively(Element *elt) {
  //for( int i = 0; i < elt->gen(); i++) cout << "  ";
  //cout << "trimming " << *elt << endl;
  // if we are NOT a nearbest-leaf, first traverse children
  // can I put this after the child traversal?
  if( _nearbestleaves.find(elt) == _nearbestleaves.end()) {
    trimRecursively(elt->left());
    trimRecursively(elt->right());

    if (DEBUG) {
      cout << "comparing for index " << elt->index() << "," << r(elt)
           << " values " << ehp(elt, r(elt))
           << " with " << error(elt, r(elt)) << endl;
    }
    if( Math::almost_equal(ehp(elt, r(elt)), error(elt, r(elt)), 1)) {
      trimNode(elt);
    }
  }
  r(elt, max(3, r(elt)));
}

void NearBest::setupRoot(Element *root) {
  assert(root->isRoot());
  r(root, Degree::Constant);
  eTilde(root, error(root, Degree::Constant));
  ehp(root, Degree::Constant, eTilde(root));
  ehpTilde(root, Degree::Constant, eTilde(root));

  q(root, eTilde(root));
  t(root, root);
}

void NearBest::setupLeaf(Element *leaf) {
  assert(!leaf->isRoot());
  assert(_nearbestleaves.find(leaf) != _nearbestleaves.end());

  r(leaf, Degree::Constant);
  scalar nodeErr = error(leaf, Degree::Constant);
  eTilde(leaf,  Math::recip(nodeErr, eTilde(leaf->parent())));
  ehp(leaf, Degree::Constant, error(leaf, Degree::Constant));
  ehpTilde(leaf, Degree::Constant, eTilde(leaf));
  q(leaf, eTilde(leaf));
  t(leaf, leaf);
}

void NearBest::trim() {
  //TODO some error is here; if a combined root was a nearbest leaf, we transfer
  //its nearbest-leaf state to its non-combination descendant roots. These are
  //not flagged for trimming, so no trimming is done.
  cerr << "size: " << _nearbestleaves.size() << endl;
  for( auto &root : partition.roots()) {
    trimRecursively(root);
  }
}

void NearBest::setDOFvalues() {
  for( auto &elt : _nearbestleaves) {
    partition.handler().set(elt, r(elt));
  }
}

int NearBest::countRealLeavesInSubTree(Element *e, bool parentWasLeaf) {
  if(_virtuals.find(e) != _virtuals.end()) return 0;

  bool isLeaf = parentWasLeaf || (_nearbestleaves.find(e) != _nearbestleaves.end());
  if(_combinations.find(e) != _combinations.end()) {
    return countRealLeavesInSubTree(e->left(), isLeaf) + 
           countRealLeavesInSubTree(e->right(), isLeaf);
  }
  if(isLeaf) {
    return 1;
  }
  return countRealLeavesInSubTree(e->left(), false) + 
         countRealLeavesInSubTree(e->right(), false);
}

void NearBest::print(ostream &os, Element *e) {
  //cout << "gonna print for " << e->index() << endl;
  if(_virtuals.find(e) != _virtuals.end()) return;
  if(_combinations.find(e) != _combinations.end()) {
    print(os, e->left());
    print(os, e->right());
    return;
  }
  os << Print::formatted("%d %d %d %d\n", e->i(0), e->i(1), e->i(2), e->type().toInt());//ehp(e, r(e)));
}

void NearBest::printNearBestLeafHpErrors(const std::string &filename) {
  ofstream os(filename, ofstream::trunc);
  //os << "tridim error" << endl;
  partition.printVerts(os);
  os << countRealLeavesInSubTree(_root) << endl;//partition.roots().size() << endl;
  for( auto &elt : _nearbestleaves) {//partition.roots()) {
    //os << elt->index() << endl;
    print(os, elt);
  }
  os.flush();
  os.close();
}

// returns squared error on the roots
scalar NearBest::errorOnRoots() {
  return ehp(_root, r(_root));
  scalar s = 0;
  for( auto &root: partition.roots()) {
    s += ehp(root, r(root));
  }
  //cerr << s << endl;
  return s;
}

NearBest::NearBest(Partition &partition, PiecewisePolynomial &poly, 
                   scalar epsilon, int maxN, std::ostream& os) : 
                   partition(partition), poly(poly)
{
  combineRoots();

  poly.print(std::cout);
  os << "4 " << epsilon << std::endl;
  assert(maxN > 0);
  //cout << "Nearbest step 1" << endl;
  /* step 1: setup */
  _nearbestleaves.clear();
  _nearbestleaves.insert(_root);

  //cout << "NIEUWE ROOT " << root->index() << endl;
  auto D = _root;
  setupRoot(D);
  cout << "hier dan" << endl;

  int N = 1;

  if( N < maxN && errorOnRoots() > epsilon*epsilon) {
    while(true) {
      /* step 2: create new tree */
      //cout << "NearBest Step 2" << endl;
      if (DEBUG) {
        cout << endl;
        cout << "Increasing the tree." << endl;
      }
      Element *node_N = t(D);
      if (DEBUG) {
        fprintf(stdout, "node_N = t(%lld) = %lld; total tris=%d, error=%g\n",
                D->index(), node_N->index(), partition.leaves().size(),
                (double) q(D));
      }
      //cout << *node_N << endl;
      if( node_N->isLeaf()) {
        if (DEBUG) {
          cout << "Bisecting " << *node_N << endl;
        }
        partition.bisect(node_N, true);
      }
      //printNearBestLeafHpErrors("output/floep_" + to_string(N) + ".dof" );

      _nearbestleaves.erase(node_N);
      _nearbestleaves.insert(node_N->left());
      _nearbestleaves.insert(node_N->right());

      //cout << "NearBest Step 3" << endl;
      /* step 3: update errors or new leaves */
      setupLeaf(node_N->left());
      setupLeaf(node_N->right());

      //cout << "NearBest Step 4" << endl;
      /* step 4: update current node */
      Element *node = node_N;

      //cout << "NearBest Step 5" << endl;
      /* step 5: check for termination */
      N += 1;
      //cerr << "We have N=" << N << " and total error = " << errorOnRoots() << "; need " << epsilon << endl;
      if( N > maxN || errorOnRoots() < epsilon*epsilon) {
        cerr << "hiera" << endl;
        os << "3 N=" << N << " and total error = " << sqrt(errorOnRoots()) << endl;
        break;
      }

      //cout << "NearBest Step 6" << endl;
      while(true) {
        if (DEBUG) {
          cout << endl;
        }
        /* step 6a: increase DoF */
        //cout << "Step 6a" << endl;
        r(node, r(node) + 1);
        //if(poly.has(node) && (r(node) > poly.dim(node))) cerr << "GREATER" << endl;
        if (DEBUG) {
          cout << "r(" << node->index() << ") = " << r(node) << endl;
        }
        /* step 6b: compute children of current node */
        //cout << "Step 6b" << endl;
        Element *left = node->left();
        Element *right = node->right();
        if (DEBUG) {
          cout << "left = " << left->index() << "; right = " << right->index()
               << endl;
        }

        /* step 6c: compute hp-error 
         * TODO: how does it change when r doesn't always decrease the error? */
        //cout << "Step 6c" << endl;
        scalar newError = min(ehp(left, r(left)) + ehp(right, r(right)), error(node, r(node)));
        if (DEBUG) {
          fprintf(stdout, "r(%lld) = %d; r(%lld) = %d; r(%lld) = %d\n", 
                 node->index(), r(node),
                 node->left()->index(), r(node->left()),
                 node->right()->index(), r(node->right()));
          fprintf(stdout, "e_hp^r(%lld) = min(%g + %g, %g) = %g\n", 
                 node->index(), (double) ehp(left, r(left)),
                 (double) ehp(right, r(right)), (double) error(node, r(node)),
                 (double) newError);
        }
        ehp(node, r(node), newError);

        /* step 6d: set tilde hp-error */
        //cout << "Step 6d" << endl;
        scalar newErrorTilde = Math::recip(ehp(node, r(node)), ehpTilde(node, r(node)-1));
        if (DEBUG) {
          printf("~e_hp^r(%lld) = 1/(1/%g + 1/%g) = %g\n",
              node->index(), (double) ehp(node, r(node)),
              (double) ehpTilde(node, r(node)-1), (double) newErrorTilde);
        }
        ehpTilde(node, r(node), newErrorTilde);

        /* step 6e: TODO what did this do again? */
        //cout << "Step 6e" << endl;
        Element *X = left;
        if( q(right) > q(left)) X = right;
        q(node, min(q(X), ehpTilde(node, r(node))));
        t(node, t(X));

        if (DEBUG) {
          fprintf(stdout, "X = argmax(q(%lld), q(%lld)) = argmax(%g, %g) = %lld\n",
              left->index(), right->index(), (double) q(left), (double) q(right),
              X->index());
          fprintf(stdout, "q(%lld) = min(%g, %g) = %g\n",
              node->index(), (double) q(X), (double) ehpTilde(node, r(node)),
              (double) q(node));
          fprintf(stdout, "t(%lld) = t(%lld) = %lld\n", node->index(), X->index(),
              t(node)->index());
        }

        /* step 6f: continue */
        //cout << "Step 6f" << endl;
        if( node != D) {
          if (DEBUG) {
            cout << "Setting node from " << node->index() << " to its parent "
                 << node->parent()->index() << endl;
          }
          node = node->parent();
        } else break;
      }
    }
  }

  if (DEBUG) {
    cout << "hier benm ik" << endl;
  }
  removeCombinedRoots();

  if (DEBUG) {
    cout << "hier dan" << endl;
  }
  //printNearBestLeafHpErrors("output/floep_" + to_string(N) + ".dof" );
  //assert(!partition.hasOverlap(_nearbestleaves));

  if (DEBUG) {
    cout << "hier dan" << endl;
  }
  /* trimming step: create hp-tree from large h-tree */
  //cout << "Trimming step" << endl;
  trim();
  if (DEBUG) {
    cerr << "Trimming done; total leaves " << _nearbestleaves.size() << " "
         << partition.leaves().size() << endl;
  }

  setDOFvalues();
  cout << "NearBest done" << endl;
}
