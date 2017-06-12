#pragma once

#include <iostream>
#include <cassert>
#include <vector>

#include "degree.h"
#include "hasdof.h"

class DOFHandler;

class Dofs : public HasDOF {
protected:
  std::vector<int> _g;

public:
  Dofs() = default;

  Dofs(int dim) : _g(dim, -1) {}

  int dof() const override { return _g.size(); }
  std::vector<int> values() const { return _g; }

  void reset(int dim) { _g.assign(dim, -1); }
  void resize(int dim) { assert(dim >= size()); _g.resize(dim, -1);}
  bool has() const;
  int size() const;

  void set(int local, int global);
  int get(int local) const;
  bool is(int local) const;

  int getVertex(int local) const;
  void setVertex(int local, int global);

  int getEdge(int degree, int edge) const;
  void setEdge(int degree, int edge, int global);

  int getFace(int degree, int r1) const;
  void setFace(int degree, int r1, int global);

  void print() const;

  friend class DOFHandler;
};
