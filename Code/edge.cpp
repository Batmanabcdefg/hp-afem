#include "edge.h"

using namespace std;

ostream & operator<<( ostream &os, Edge e) {
  os << *(e.v0) << "," << *(e.v1);
  return os;
}
