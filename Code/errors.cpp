#include <cassert>

#include "errors.h"

using namespace std;

Errors::Errors(const ElementScalarSet &errors)
{
  for(auto p : errors) {
    _errors.insert(p);
  }
}

ElementScalarSet Errors::on(const ElementSet &elts)
{
  ElementScalarSet ret;
  for( auto elt : elts) {
    ret.insert(make_pair(elt, on(elt)));
  }

  return ret;
}

scalar Errors::on(Element *elt)
{
  if(_errors.count(elt) > 0) return _errors[elt];
  assert(!elt->isLeaf());

  scalar ret = on(elt->left()) + on(elt->right());
  _errors.insert(make_pair(elt, ret));
  
  return ret;
}
