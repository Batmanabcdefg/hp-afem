#pragma once

class Index {
  long long _index, _root;

public:
  Index(long long index, long long root) : _index(index), _root(root) {}

  Index left() const {
    return Index(2L*(_index + 1L) - 1L, _root);
  }

  Index right() const {
    return Index(2L*(_index + 1L), _root);
  }

  Index root() const {
    return Index(0, _root);
  }

  bool operator == (Index &other) {
    return (_root == other._root) && (_index == other._index);
  }

  bool operator < (Index &other) {
    return (_index < other._index) || (_index == other._index && _root < other._root);
  }

  std::string toString() const {
    return std::to_string(_root) << "," << std::to_string(_index);
  }

  friend std::ostream &operator<<(std::ostream& out, const Index &i) {
    os << i.toString();
    return os;
  }
};
