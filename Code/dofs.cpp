#include <iostream>

#include "dofs.h"
#include "math.h"

using namespace std;

void Dofs::set(int local, int global) {
  assert(local < _g.size());
  _g[local] = global;
}

bool Dofs::has() const {
  for( auto i : _g) if(i != -1) return true;
  return false;
}

int Dofs::get( int local) const {
  assert(local < _g.size());
  return _g[local];
}

bool Dofs::is( int local) const {
  return _g[local] != -1;
}

int Dofs::size() const {
  return _g.size();
}

int  Dofs::getVertex(int local) const {
  assert(local < 3);
  return get(local);
}

void Dofs::setVertex(int local, int global) {
  assert(local < 3);
  return set(local, global);
}

int  Dofs::getEdge(int degree, int edge) const {
  return get(Math::triNum(degree) + edge);
}

void Dofs::setEdge(int degree, int edge, int global) {
  return set(Math::triNum(degree) + edge, global);
}

int  Dofs::getFace(int degree, int r1) const {
  return get(Math::triNum(degree) - (degree - 2) + r1);
}

void Dofs::setFace(int degree, int r1, int global) {
  return set(Math::triNum(degree) - (degree - 2) + r1, global);
}

void Dofs::print() const {
  cout << "  .   .   .";
  for( int k = 1; k < degree(); k++) {
    cout << " - |   |   |";
    if( k >= 2) {
      for( int j = 0; j < k-1; j++) {
        cout << "   O";
      }
    }
  }
  std::cout << std::endl;
  for( auto gi : _g) {
    printf("%3d ", gi);
  }
  std::cout << std::endl;
}
