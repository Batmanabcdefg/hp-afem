#pragma once

namespace Precomputations {
class Tritype {
  public:
    friend std::ostream & operator<<( std::ostream &os, Tritype tt) {
      os << (tt.x > 0) << " " << (tt.y > 0) << " " << (tt.z > 0);
      return os;
    }

    bool x, y, z, _isset;
    bool i(int idx) { return (toInt() >> idx) % 2; }
    int toInt()     { return 4*z + 2*y + x; }
    Tritype left()  { return _isset ? Tritype(!x, x, y) : Tritype(); }
    Tritype right() { return _isset ? Tritype(y, z, x)  : Tritype(); }
    Tritype(bool _x, bool _y, bool _z) : x(_x), y(_y), z(_z), _isset(true) {}
    Tritype(int i) : x((i >> 0) % 2), y((i >> 1) % 2), z((i >> 2) % 2), _isset(true) {}
    Tritype() : _isset(false) {}
};
}
