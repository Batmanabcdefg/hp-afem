#pragma once

#include "basis.h"
#include "matrix.h"

class Element;

class ElementMatrices {
  Bases &_bases;
  Element *_root;
  Matrix _eltmats[8][4];
public:
  ElementMatrices(Bases &bases, Element *root) : _bases(bases), _root(root) {}
  Matrix &get(int tt, int tc);
  virtual ~ElementMatrices() = default;
};
