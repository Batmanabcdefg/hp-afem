#include <queue>

#include "fem/solution.h"
#include "math.h"
#include "print.h"
#include "solvable.h"
#include "poly.h"

using namespace std;

/**
 * This routine solves the sparse linear system Ax=b, where A is the stiffness 
 * matrix, and b the load vector.
 *
 * We do this by assembling A and b from the element matrix and 
 * element vector on each element.  This involves converting local degrees of
 * freedom to the global ones.
 */
std::unique_ptr<FEM::Solver> Solvable::solve() {
  assert(handler().valid());
  auto solver = std::make_unique<FEM::Solver>(handler(), _rhs);
  _sol = solver->sol();
  return solver;
}

ostream &Solvable::printLinearInterpolant( ostream &os) {
  os << Print::formatted("tridim linsol") << endl;
  os << Print::formatted("%lu", _verts.size()) << endl;
  for(auto v : _verts) os << Print::formatted("%.16f %.16f", (double) v->x, (double) v->y) << endl;

  os << Print::formatted("%lu", _sol.definedOn().size()) << endl;
  for(Element *elt : _sol.definedOn()) {
    os << Print::formatted("%lu %d %d %d",
        _sol.local_dim(elt), elt->i(0), elt->i(1), elt->i(2)) << endl;
  }

  os << Print::formatted("%lu", leaves().size()) << endl;
  for(Element *elt : leaves()) {
    Vector local = _sol.locallyAt(elt).head(3);
    os << Print::formatted("%lu %d %d %d",
        _sol.local_dim(elt), elt->i(0), elt->i(1), elt->i(2));
    for( int i = 0; i < 3; i++) {
      os << Print::formatted(" %.16f", (double) local[i]);
    }
    os << endl;
  }

  return os;
}

ostream &Solvable::printDOFs( ostream &os) {
  os << "tridim tritype dof" << endl;
  os << _verts.size() << endl;
  for(auto v : _verts) {
    os << v->x << " " << v->y << " " << handler().find(v) << endl;
  }
  handler().print(os);
  return os;
}

/**
 * TODO!
 */
ostream &Solvable::printSolution( ostream &os) {
  os << Print::formatted("tridim tritype sol") << endl;
  os << Print::formatted("%lu", _verts.size()) << endl;
  for(auto v : _verts) {
    os << Print::formatted("%.16f %.16f %d",
        (double) v->x, (double) v->y, v->isBoundary()) << endl;
  }

  _sol.print(os);
  return os;
}

ostream &Solvable::printRhsMesh( ostream &os) {
  os << Print::formatted("tridim tritype sol") << endl;
  os << Print::formatted("%lu", _verts.size()) << endl;
  for(auto v : _verts) {
    os << Print::formatted("%.16f %.16f %d",
        (double) v->x, (double) v->y, v->isBoundary()) << endl;
  }

  _rhs.print(leaves(), os);
  return os;
}

ostream &Solvable::printElementScalarSet(const ElementScalarSet &on, ostream &os) {
  os << Print::formatted("tridim error") << endl;
  os << Print::formatted("%lu", _verts.size()) << endl;
  for(auto v : _verts) {
    os << Print::formatted("%.16f %.16f %d",
        (double) v->x, (double) v->y, v->isBoundary()) << endl;
  }


  os << Print::formatted("%lu", on.size()) << endl;
  for( auto &eltdbl : on) {
    Element *elt = eltdbl.first;
    os << Print::formatted("%d %d %d %d %g\n",
        _sol.local_dim(elt), elt->i(0), elt->i(1), elt->i(2),
        (double) eltdbl.second);
  }
  return os;
}

ostream &Solvable::printVerts( ostream &os) {
  os << _verts.size() << endl;
  for(auto v : _verts) {
    os << Print::formatted("%.16f %.16f", (double) v->x, (double) v->y) << endl;
  }
  return os;
}

ostream &Solvable::printLeaves( ostream &os) {
  os << "tritype" << endl;
  printVerts(os);
  os << Print::formatted("%lu", _leaves.size()) << endl;
  for(auto &elt : _leaves) {
    os << elt->i(0) << " " << elt->i(1) << " " << elt->i(2) << " "
       << elt->type().toInt() << endl;
  }
  return os;
}

ostream &Solvable::printHpMesh( ostream &os) {
  os << "tritype tridim" << endl;
  printVerts(os);
  os << Print::formatted("%lu", _leaves.size()) << endl;
  for(auto &elt : _leaves) {
    os << handler().find(elt).dim() 
       << " " << elt->i(0) 
       << " " << elt->i(1) 
       << " " << elt->i(2) 
       << " " << elt->type().toInt() 
       << endl;
  }

  return os;
}

void Solvable::printHpMesh(string filename) {
  ofstream os(filename, ofstream::trunc);
  printHpMesh(os);
  os.flush();
  os.close();
}

void Solvable::printLinearInterpolant(string filename) {
  ofstream os(filename, ofstream::trunc);
  printLinearInterpolant(os);
  os.flush();
  os.close();
}

void Solvable::printLeaves(string filename) {
  ofstream os(filename, ofstream::trunc);
  printLeaves(os);
  os.flush();
  os.close();
}

void Solvable::printRhsMesh(string filename) {
  ofstream os(filename, ofstream::trunc);
  printRhsMesh(os);
  os.flush();
  os.close();
}

void Solvable::printDOFs(string filename) {
  ofstream os(filename, ofstream::trunc);
  printDOFs(os);
  os.flush();
  os.close();
}

void Solvable::printSolution(string filename) {
  ofstream os(filename, ofstream::trunc);
  printSolution(os);
  os.flush();
  os.close();
}

void Solvable::printElementScalarSet(
    const ElementScalarSet &on, string filename) {
  ofstream os(filename, ofstream::trunc);
  printElementScalarSet(on, os);
  os.flush();
  os.close();
}
