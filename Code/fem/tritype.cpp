#include "tritype.h"

using namespace std;

ostream & operator<<( ostream &os, TriType tt) {
  if( tt._isset) os << tt.x << tt.y << tt.z;
  else os << "(not set)";
  return os;
}
