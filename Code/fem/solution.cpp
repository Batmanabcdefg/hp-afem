#include "solution.h"

#include <algorithm>
#include <fstream>
#include <vector>

#include "../matrix.h"

using namespace std;

std::vector<std::string> split(const std::string &text, char sep) {
  std::vector<std::string> tokens;
  std::size_t start = 0, end = 0;
  while ((end = text.find(sep, start)) != std::string::npos) {
    tokens.push_back(text.substr(start, end - start));
    start = end + 1;
  }
  tokens.push_back(text.substr(start));
  return tokens;
}

namespace FEM {
Solution::Solution(Vector &sol, const ElementDofsSet &elts, const Rhs &rhs) :
    _global(sol), _rhs(rhs) {
  for( auto &p : elts) insert_dofs(p.first, p.second, true);
}

Solution::Solution(Vector &sol, const ElementDofsMap &elts, const Rhs &rhs) :
  _global(sol), _g(elts), _rhs(rhs)
{
  for( auto &eltdof : elts) 
    insert_vector(eltdof.first, indexSol(eltdof.second), true);
}

Vector Solution::indexSol(const Dofs &g) {
  Vector vec(g.size());

  for( int i = 0; i < g.size(); i++) {
    vec[i] = (g.get(i) == -1) ? 0.0 : _global[g.get(i)];
  }
  return vec;
}

void Solution::erase(Element *elt) {
  PiecewisePolynomial::erase(elt);
  
  // it's possible we don't have a set of dofs for this element as it was copied
  auto it = _g.find(elt);
  if( it != _g.end()) _g.erase(it);
}

void Solution::insert_dofs(Element *elt, Dofs local, bool definedOn) {
  assert(_g.count(elt) == 0);
  _g.insert(make_pair(elt, local));

  insert_vector(elt, indexSol(local), definedOn);
}

Solution Solution::fromSolFile(string filename, ElementDofsSet elts, Rhs &rhs) {
  ifstream solfile(filename);
  assert(solfile.good());

  string settingString;
  getline(solfile, settingString);
  vector<string> settings = ::split(settingString, ' ');
  auto it = find(settings.begin(), settings.end(), "sol");
  assert(it != settings.end());

  int nVerts, nElts, on_boundary;
  scalar xVert, yVert;
  solfile >> nVerts;
  for( int i = 0; i < nVerts; i++) solfile >> xVert >> yVert >> on_boundary;
  solfile >> nElts;
  assert(nElts == elts.size());

  int maxDof = 0;
  for( auto p : elts) for( int gi : p.second.values()) maxDof = max(maxDof, gi);
  Vector sol = Vector::Zero(maxDof + 1);

  for( auto eltdofs : elts) {
    Element *elt = eltdofs.first;
    Dofs dofs = eltdofs.second;
    int nDofs;
    solfile >> nDofs;
    assert(nDofs == dofs.size());
    int v1, v2, v3, tt;
    solfile >> v1 >> v2 >> v3 >> tt;
    assert( v1 == elt->i(0) && v2 == elt->i(1) && v3 == elt->i(2) 
         && tt == elt->type().toInt());

    for( auto dof : dofs.values()) {
      scalar solval;
      solfile >> solval;
      if( dof == -1) {
        assert(solval == 0.0);
      } else {
        sol[dof] = solval;
      }
    }
  }

  return Solution(sol, elts, rhs);
}

Vector Solution::residualVector() { 
  cout << _systemMat.rows() << " " << _global.rows() << " " << _systemRhs.rows() << endl;
  return _systemMat * _global - _systemRhs; 
}
}
