#pragma once

#include <iostream>
#include <array>
#include <vector>
#include <cassert>

#include "basis.h"
#include "dofs.h"
#include "degree.h"
#include "basis.h"
#include "math.h"
#include "tritype.h"
#include "triclass.h"
#include "matrix.h"
#include "system.h"
#include "triangle.h"
#include "edge.h"
#include "node.h"
#include "elementmatrices.h"

class ElementTree;
class Refinable;
class Solvable;
class Approximator;
class Partition;
class Matchable;
class NearBest;

class Element : public Triangle, public Node {
  public:
    friend std::ostream & operator<<( std::ostream &os, Element &elt);

    void printRecursive();

    const Basis *basis() { return _basis; }
    TriType type() { return (_basis == nullptr) ? TriType() : _basis->type(); }
    const TriClass &triclass();
    void setRoot() { _root = this; _parent = nullptr; _gen = 0; _tc = TriClass(0); }

    virtual Element *left() const { return (Element *)Node::left(); }
    virtual Element *right() const { return (Element *)Node::right(); } 
    virtual Element *parent() const { return (Element *)Node::parent(); }
    virtual Element *root() const { return (Element *)Node::root(); }
    virtual Element *child(bool rightChild) const { return (Element *)Node::child(rightChild); }

    Edge bisectionEdge() { return Edge(_v[1], _v[2]); }

    /* Element matrix & vector routines */
    Matrix elementMatrix( int dof);
    ElementMatrix massMatrix( int dof, int dof2 = -1);
    ElementVector elementVector( int dof);
    Matrix transferMatrix( bool rightChild, int dof);

    Element(Vertex *v[3], Element *parent, Element *left, Element *right, const Basis *basis);
    Element(Vertex *v[3], Element *parent, Element *left, Element *right);
    Element(Vertex *v[3], const Basis *basis);

    virtual ~Element() = default;

    friend class ElementTree;
    friend class Refinable;
    friend class Solvable;
    friend class Approximator;
    friend class Partition;
    friend class Matchable;
    friend class NearBest;

  protected:
    const static int _childclass[4][2];

    /* the elements class, i.e. the equivalence class this element belongs to */
    TriClass _tc;

    /* a pointer to this elements basis */
    const Basis *_basis = nullptr;

    ElementMatrices *_eltmats;

    // for NearBest (blegh hacky)
    Element(Element *left, Element *right, long long index) : 
      Triangle(), Node(nullptr, left, right) { _index = index; }
};
