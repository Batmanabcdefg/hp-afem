#pragma once
#include <iostream>
#include "vertex.h"

class Edge {
  protected:
    Vertex *v0, *v1;

  public:
    friend std::ostream & operator<<( std::ostream &os, Edge e);

    Vertex *v(int index) { return (index == 0) ? v0 : ((index == 1) ? v1 : nullptr); }

    bool isOpposite(Edge e) { return v0->isSameAs(e.v(1)) && v1->isSameAs(e.v(0)); }
    Edge( Vertex *_v0, Vertex *_v1 ) : v0(_v0), v1(_v1) {}
};
