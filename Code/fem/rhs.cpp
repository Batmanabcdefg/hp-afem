#include "rhs.h"

#include <cassert>

using namespace std;

namespace FEM {

void Rhs::retainOnly(ElementSet &elts) {
  ElementPairMap<Vector> newl2g;
  for(auto &elt : elts) {
    assert(has(elt));
    newl2g.insert(*l2g.find(elt));
  }

  l2g = newl2g;
  _availableOn = elts;
  _definedOn = elts;
}

}
