#pragma once
#include <iostream>
#include "math.h"
#include "indexed.h"

class Vertex : public Indexed {
protected:
  int _dof = -1;
public:
  friend std::ostream & operator<<( std::ostream &os, const Vertex &v);

  bool isSameAs(scalar x, scalar y);
  bool isSameAs(Vertex *v) { return isSameAs(v->x, v->y); }

  Vertex operator - (const Vertex &v) const { return Vertex(x-v.x, y-v.y); }
  scalar norm() { return hypot(x, y); }

  bool hasDOF() { return _dof != -1; }
  int getDOF() { return _dof; }
  void setDOF(int global) { _dof = global; }
  void resetDOF() { setDOF(-1); }

  bool _isBoundary = false;
  bool isBoundary() { return _isBoundary; }

  scalar x, y;
  scalar cross(Vertex v) { return v.x*y - v.y*x; }
  Vertex( scalar _x, scalar _y) : Indexed(), x(_x), y(_y) {}
};
