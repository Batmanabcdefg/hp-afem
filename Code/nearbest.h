#pragma once

#include <fstream>
#include <vector>
#include <set>

#include "matrix.h"
#include "partition.h"
#include "triangleset.h"

/**
 *  hp-NearBest algorithm. Given a partition, a piecewise polynomial to approximate,
 *  an error bound epsilon, and a maximum number of iterations, this routine
 *  creates a (non-conforming) refinement of the given partition such that
 *  the total error
 *    E_{partition} := \sum_{element (K,dof) in partition} e_{K,dof}(poly|_K)
 *  where 
 *    e_{K,dof}(w) := = inf |Q_dof(w) - w|^2_{H^1_0(K)}
 *  satisfies
 *    E_{partition} <= epsilon*epsilon.
 */

class NearBest {
private:

  Partition &partition;
  PiecewisePolynomial &poly;

  ElementSet _nearbestleaves;
  ElementSet _virtuals;
  ElementSet _combinations;

  Element *_root;

  ElementPairMap<std::map<int, scalar>> _e;
  ElementPairMap<std::map<int, scalar>> _ehp;
  ElementPairMap<std::map<int, scalar>> _ehpTilde;
  ElementPairMap<scalar> _eTilde;
  ElementPairMap<int> _r;
  ElementPairMap<Element *> _t;
  ElementPairMap<scalar> _q;

  int r_forced(Element *elt, int val = -1);
  void trimNode(Element *elt);
  void trimRecursively(Element *elt);
  void setupRoot(Element *root);
  void setupLeaf(Element *leaf);
  void eraseFromDataStructures(Element *elt);
  void combineRoots();
  void removeCombinedRoots();
  void removeCombinedRootsRecursive(Element *elt);

  void print(std::ostream &os, Element *e);
  int countRealLeavesInSubTree(Element *e, bool parentWasLeaf = false);

protected:
  scalar error(Element *elt, int dof);
  int r(Element *elt, int val = -1);

  scalar ehp(Element *elt, int dof, scalar val = -1);
  scalar ehpTilde(Element *elt, int dof, scalar val = -1);
  scalar eTilde(Element *elt, scalar val = -1);
  Element *t(Element *elt, Element *val = nullptr);
  scalar q(Element *elt, scalar val = -1);

  void trim();
  void setDOFvalues();

public:
  NearBest(Partition &partition, PiecewisePolynomial &poly, scalar epsilon, int maxN, std::ostream& os);
  NearBest(Partition &partition, PiecewisePolynomial &poly, scalar epsilon, int maxN) :
    NearBest(partition, poly, epsilon, maxN, std::cerr) {}

  void printNearBestLeafHpErrors(const std::string &filename);
  scalar errorOnRoots();
};
