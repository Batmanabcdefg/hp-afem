#pragma once
#include <iostream>

/**
 *  Tritype.h
 *
 *  A triangle type or tritype defines the orientation of the edge basis functions
 *  along each edge.  The edge to the right of the newest vertex is 0; after that,
 *  the edges are numbered in ccw order.  On each edge, the orientation can either
 *  be forward (bool false) or backward/reversed (bool true).  This yields 2^3=8
 *  different possible triangle types (x,y,z).  We go to and from integers by
 *  (x,y,z) = 0bzyx, so that the orientation along edge i is equal to 
 *  ((0bzyx) >> i) % 2.
 */
class TriType {
  public:
    friend std::ostream & operator<<( std::ostream &os, TriType tt);

    bool x, y, z, _isset;
    bool i(int idx) const { return (toInt() >> idx) % 2; }
    int toInt()     const { return 4*z + 2*y + x; }
    operator int () const { return toInt(); }
    TriType left()  const { return _isset ? TriType(!x, x, y) : TriType(); }
    TriType right() const { return _isset ? TriType(y, z, x)  : TriType(); }
    TriType(bool _x, bool _y, bool _z) : x(_x), y(_y), z(_z), _isset(true) {}
    TriType(int i) : x((i >> 0) % 2), y((i >> 1) % 2), z((i >> 2) % 2), _isset(true) {}
    TriType() : _isset(false) {}
};
