#include "element.h"

using namespace std;

ostream & operator<<( ostream &os, Element &elt) {
  os << "tri " << elt.index() << ":  [" << *elt.v(0) << "," << *elt.v(1) 
                                 << "," << *elt.v(2) << "]; ";
  os << "parent " << ((elt.parent() == nullptr) ? -1 : elt.parent()->index());
  os << "; left " << ((elt.left() == nullptr) ? -1 : elt.left()->index());
  os << "; right " << ((elt.right() == nullptr) ? -1 : elt.right()->index());
  os << "; type " << elt.type() << endl;
  return os;
}

const TriClass &Element::triclass() {
  if( _tc.valid()) return _tc;
  if( parent()) {
    if(parent()->right() == this) {
      return _tc = parent()->triclass().right();
    } else {
      return _tc = parent()->triclass().left();
    }
  }
  return _tc = TriClass(0);
}

Matrix Element::elementMatrix( int dof) {
  int curdim = Degree::dofToDim(dof);
  return _eltmats->get(type(), triclass()).topLeftCorner(curdim, curdim);
}

ElementMatrix Element::massMatrix( int dof, int dof2) {
  int curdim = Degree::dofToDim(dof);
  if( dof2 == -1) dof2 = dof;
  int curdim2 = Degree::dofToDim(dof2);
  assert(curdim2 <= curdim);

  scalar D = 2.0L*vol();
  return basis()->massmat().topLeftCorner(curdim,curdim2)*D;
}

ElementVector Element::elementVector( int dof) {
  int curdim = Degree::dofToDim(dof);
  scalar D = 2.0L * vol();
  return basis()->eltvec().head(curdim)*D;
}

Matrix Element::transferMatrix( bool rightChild, int dof) {
  int curdim = Degree::dofToDim(dof);
  assert(basis()->dim() >= curdim);
  return basis()->transfermat((int)rightChild).topLeftCorner(curdim, curdim);
}

Element::Element(Vertex *v[3], Element *parent, Element *left, Element *right, const Basis *basis) : 
  Triangle(v, parent ? parent->vol()/2 : -1.0), 
  Node(parent, left, right),
  _basis(basis), 
  _eltmats(parent ? parent->_eltmats : nullptr)
{}

Element::Element(Vertex *v[3], Element *parent, Element *left, Element *right) : 
  Element::Element(v, parent, left, right, nullptr) {}

Element::Element(Vertex *v[3], const Basis *basis) : 
  Element::Element(v, nullptr, nullptr, nullptr, basis) {} 

void Element::printRecursive() {
  cout << *this << endl;
  if( !isLeaf()) {
    left()->printRecursive();
    right()->printRecursive();
  }
}
