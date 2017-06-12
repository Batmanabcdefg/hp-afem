#pragma once
#include <iostream>
#include <array>
#include "indexed.h"
#include "vertex.h"
#include "edge.h"

class Matchable;

class Triangle : public Indexed {
  protected:
    Vertex *_v[3];
    scalar _vol;

    Triangle() = default;

  public:
    std::array<Edge, 3> edges();
    std::array<Vertex *, 3> verts();
    Vertex **v() { return _v; }

    scalar vol() const { return _vol; }
    int i(int index) const { return (index < 3 && index >= 0) ? _v[index]->index() : -1; }
    Vertex *v(int index) const { return (index < 3 && index >= 0) ? _v[index] : nullptr; }
    Triangle(Vertex *v[3], scalar _vol);
    Triangle(Vertex *v[3]);

    bool contains(const Vertex *u) const {
      scalar lambda0 = (_v[0]->y*_v[2]->x - _v[0]->x*_v[2]->y +
          (_v[2]->y - _v[0]->y)*u->x + (_v[0]->x - _v[2]->x)*u->y)/(2*vol());
      scalar lambda1 = (_v[0]->x*_v[1]->y - _v[0]->y*_v[1]->x +
          (_v[0]->y - _v[1]->y)*u->x + (_v[1]->x - _v[0]->x)*u->y)/(2*vol());
      scalar lambda2 = 1 - lambda0 - lambda1;
      //std::cout << lambda0 << " " << lambda1 << " " << lambda2 << std::endl;
      return (lambda0 >= -0.1) && (lambda1 >= -0.1) && (lambda2 >= -0.1);
    }

    bool contains(const Triangle* t) {
      return contains(t->v(0)) && contains(t->v(1)) && contains(t->v(2));
    }

    friend class Matchable;
};
