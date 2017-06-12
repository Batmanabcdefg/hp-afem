#pragma once

class TriClass {
  const static int _childclass[4][2];
  int _class = -1;
public:
  TriClass() : _class(-1) {}
  TriClass(int tc) : _class(tc) {}

  bool valid() const { return _class != -1; }
  operator int () const { return _class; }

  TriClass left()  const { return TriClass(_childclass[_class][0]); }
  TriClass right() const { return TriClass(_childclass[_class][0]); }
};
