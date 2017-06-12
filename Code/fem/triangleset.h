#pragma once
#include <set>

#include "triangle.h"
#include "dofs.h"

class Element;

struct TriangleSetCompare {
  bool operator() (const Triangle *me, const Triangle *you) const {
    return me->index() < you->index();
  }
};

template<class T>
struct TrianglePairSetCompare {
  bool operator() (const std::pair<Triangle *, T> me, const std::pair<Triangle *, T> you) const {
    return me.first->index() < you.first->index();
  }
};

template <class T>
using ElementPairSet = std::set<std::pair<Element *, T>, TrianglePairSetCompare<T>>;

typedef std::set<Element *, TriangleSetCompare> ElementSet;
typedef ElementPairSet<Dofs> ElementDofsSet;
typedef ElementPairSet<int> ElementDimsSet;
typedef ElementPairSet<scalar> ElementScalarSet;
typedef std::set<std::pair<int,int>> IndexDofSet;
typedef std::set<int> IndexSet;

template <class T>
using ElementPairMap = std::map<Element *, T, TriangleSetCompare>;
typedef ElementPairMap<Dofs> ElementDofsMap;

class TriangleSet {
public:
  static ElementSet unionSets(ElementSet s1, ElementSet s2) {
    ElementSet res = s1;
    res.insert(s2.begin(), s2.end());
    return res;
  }

  static void unionSetsInto(ElementSet &s1, const ElementSet &&s2) {
    s1.insert(s2.begin(), s2.end());
  }

};
