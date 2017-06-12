#pragma once

#include <map>

#include "triangleset.h"
#include "element.h"

class Errors {
public:
  Errors(const ElementScalarSet &errors);
  Errors() {}

  ElementScalarSet on(const ElementSet &elts);
  scalar on(Element *elt);

private:
  std::map<Element *, scalar> _errors;
};
