#include "vertex.h"

using namespace std;

ostream & operator<<( ostream &os, const Vertex &v) {
  os << "[" << v.x << "," << v.y << "] (" << v.index() << ")"; 
  return os; 
}

bool Vertex::isSameAs( scalar xx, scalar yy) {
  return Math::almost_equal(x, xx, 1) && Math::almost_equal(y, yy, 1); 
}
