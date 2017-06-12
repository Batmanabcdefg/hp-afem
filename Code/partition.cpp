#include <fstream>

#include "degree.h"
#include "partition.h"
#include "fem/rhs.h"

using namespace std;

std::vector<std::string> splitt(const std::string &text, char sep) {
  std::vector<std::string> tokens;
  std::size_t start = 0, end = 0;
  while ((end = text.find(sep, start)) != std::string::npos) {
    tokens.push_back(text.substr(start, end - start));
    start = end + 1;
  }
  tokens.push_back(text.substr(start));
  return tokens;
}

Partition::Partition(string basisdir, string meshfn)
  : Matchable(std::move(basisdir))
{
  bool reading_solution = false;

  ifstream meshfile(meshfn);
  assert(meshfile.good());

  if(meshfn.compare(meshfn.size()-4, 4, ".sol") == 0) {
    reading_solution = true;
  }

  if (reading_solution) {
    string settingString;
    getline(meshfile, settingString);
    vector<string> settings = ::splitt(settingString, ' ');
    auto it = find(settings.begin(), settings.end(), "sol");
    assert(it != settings.end());
  }

  int nVerts, nElements;
  meshfile >> nVerts;
  for (auto i = 0; i < nVerts; i++) {
    scalar x, y;
    meshfile >> x >> y;
    auto *vert = new Vertex(x,y);
    if(reading_solution) {
      bool on_boundary;
      meshfile >> on_boundary;
      vert->_isBoundary = on_boundary;
    }
    addVertex(vert);
  }

  meshfile >> nElements;
  ElementDimsSet current;
  for (auto i = 0; i < nElements; i++) {
    int ndof;
    if (reading_solution) {
      meshfile >> ndof;
    }
    int i0, i1, i2, type;
    meshfile >> i0 >> i1 >> i2 >> type;
    Vertex *v[3] = {_verts[i0], _verts[i1], _verts[i2]};
    Element *root = new Element(v, &_bases.basis(type));
    addElement(root);
    addRoot(root);
    addLeaf(root);
    if (reading_solution) {
      current.insert(make_pair(root, ndof));
      scalar value;
      for (int j = 0; j < Degree::dofToDim(ndof); j++) {
        meshfile >> value;
      }
    }
  }

  if (reading_solution) {
    handler().determine(current);
  }

  meshfile.close();
}

Partition::Partition(string basisdir, string meshfn, string rhsfn)
  : Partition(std::move(basisdir), std::move(meshfn))
{
  ifstream rhsfile(rhsfn);
  assert(rhsfile.good());

  int nRhsElements, nRhs;
  rhsfile >> nRhsElements >> nRhs;
  assert(nRhsElements == numElements());
  assert(nRhs <= _bases.dim());

  _rhs = FEM::Rhs(nRhs);

  for( int i = 0; i < nRhsElements; i++) {
    //create vector of largest size a RHS can have; initialize to zero to null
    //out every component after nRhs
    Vector rhs = Vector::Zero(nRhs);
    for( int j = 0; j < nRhs; j++) {
      rhsfile >> rhs[j];
    }
    _rhs.insert_vector(element(i), rhs, true);
  }

  resetRoots();
  rebuildElementFinder();
  makeConform();
  makeMatching();
  determineBoundaryVertices();
  makeTriTypesCorrect();
  handler().determine(leaves(), Degree::Linear);
}
