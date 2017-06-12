#include <Eigen/Sparse>
#include <vector>
#include <set>

#include "solver.h"
#include "../element.h"

#define LDLT 1

using namespace std;

namespace FEM {

int Solver::estimateNonZeros() {
  int estimated = 0;
  for( auto &pair : _elts) {
    int curdim = pair.second.dim();
    estimated += curdim * curdim;
  }
  return estimated;
}

void Solver::computeSystemMatrix() {
  int estimated = estimateNonZeros();
  TripVec triplets;
  triplets.reserve(estimated);

  for( auto &pair : _elts) {
    Element *elt = pair.first;
    Dofs dofs = pair.second;
    //cerr << elt->index() << " has dof " << dofs.size() << endl;
    ElementMatrix eltmat = elt->elementMatrix(dofs.size());

    for( int i = 0; i < dofs.size(); i++) {
      if( !dofs.is(i)) continue;
      for( int j = 0; j < dofs.size(); j++) {
        if( !dofs.is(j)) continue;
        triplets.push_back(Trip(dofs.get(i), dofs.get(j), eltmat(i,j)));
      }
    }
  }

  _sysmat = StiffnessMatrix(numDOFs(), numDOFs());
  _sysmat.setFromTriplets(triplets.begin(), triplets.end());

  cerr << "total dofs: " << numDOFs()
       << "; total nonzeros: " << _sysmat.nonZeros()
       << "; estimated: " << estimated << endl;
}

void Solver::computeSystemRhs() {
  _sysrhs = LoadVector::Zero(numDOFs());

  /* use each element's information to build the system */
  for( auto pairr : _elts) {
    Element *elt = pairr.first;
    Dofs dofs = pairr.second;
    int curdim = min(dofs.size(), _femrhs.dim());
    ElementVector eltrhs =
      elt->massMatrix(dofs.size(), curdim)*_femrhs.locallyAt(elt).head(curdim);

    for( int i = 0; i < dofs.size(); i++) {
      if( !dofs.is(i)) continue;
      _sysrhs[dofs.get(i)] += eltrhs[i];
    }
  }
}

int Solver::numDOFs() {
  if( _numDOFs > -1) return _numDOFs;

  set<int> values;
  int maxDof = -1;
  for( auto p : _elts) {
    for( int gi : p.second.values()) {
      if( gi != -1) {
        maxDof = max(maxDof, gi);
        values.insert(gi);
      }
    }
  }
  assert(values.size() == maxDof+1);

  return _numDOFs = 1 + maxDof;
}

void Solver::solveSystem() {
#if LDLT
  Eigen::SimplicialLDLT<StiffnessMatrix> sldlt;
  sldlt.compute(_sysmat);
  Vector sol = sldlt.solve(_sysrhs);
#else
  /* solve the system using PCG */
  Eigen::ConjugateGradient<StiffnessMatrix> cg;
  cerr << cg.maxIterations() << endl;
  cerr << cg.tolerance() << endl;
  cg.compute(_sysmat);
  Vector sol = cg.solve(_sysrhs);
  cerr << cg.iterations() << endl;
#endif

  _sol = Solution(sol, _elts, _femrhs);

  _sol._systemMat = _sysmat;
  _sol._systemRhs = _sysrhs;
}

Solver::Solver(const ElementDofsMap &elts, Rhs &rhs)
    : _elts(elts), _femrhs(rhs) {
  //cout << "In solver with " << numDOFs() << " dofs" << endl;
  if( numDOFs() == 0) {
    _sol = Solution(_sysrhs, _elts, _femrhs);
  } else {
    cout << "gonna compute matrix" << endl;
    computeSystemMatrix();
    cout << "gonna compute rhs" << endl;
    computeSystemRhs();
    cout << "gonna solve" << endl;
    solveSystem();
  }
}
}
