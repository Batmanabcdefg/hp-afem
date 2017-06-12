#include <cmath>
#include <cassert>
#include "vertex.h"
#include "triangle.h"

using namespace std;

scalar vol(Vertex *v[3]) {
  scalar vol =(- v[1]->x*v[0]->y + v[2]->x*v[0]->y 
              + v[0]->x*v[1]->y - v[2]->x*v[1]->y 
              - v[0]->x*v[2]->y + v[1]->x*v[2]->y)/2;
  assert(vol > 0);
  return fabs(vol);
}

/*
ostream &operator<<( ostream &os, Triangle &tri) {
  os << "[" << tri.i(0) << "," << tri.i(1) << "," << tri.i(2) << "]";
  return os;
}
*/

array<Edge, 3> Triangle::edges() {
  array<Edge, 3> _edges = {
    {Edge(v(0), v(1)), Edge(v(1), v(2)), Edge(v(2), v(0))}
  };
  return _edges;
}

array<Vertex *, 3> Triangle::verts() {
  array<Vertex *, 3> _verts = {{_v[0], _v[1], _v[2]}};
  return _verts;
}

Triangle::Triangle(Vertex *v[3], scalar vol) : Indexed() {
  _v[0] = v[0];
  _v[1] = v[1];
  _v[2] = v[2];
  if( vol == -1.0) vol = ::vol(v);
  _vol = vol;
  assert(_vol > 0);
}

Triangle::Triangle( Vertex *v[3]) : Triangle(v, ::vol(v)) {}

