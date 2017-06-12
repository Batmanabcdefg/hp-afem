#pragma once

class ElementTree;

class Indexed {
  protected:
    long long _index = -1;

  public:
    const long long index() const { return _index; }
    Indexed( long long index) : _index(index) {}
    Indexed() : _index(-1) {}

    friend class ElementTree;
};
