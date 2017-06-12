#pragma once
#include <cassert>

class ElementTree;

class Node {
  protected:
    Node *_left;
    Node *_right;
    Node *_parent;
    Node *_root;
    int _gen;
    bool _isLeaf;

  public:
    virtual Node *left() const { return _left; }
    virtual Node *right() const { return _right; } 
    virtual Node *parent() const { return _parent; }
    virtual Node *root() const { return _root; }

    virtual Node *child(bool rightChild) const { assert( !isLeaf()); return (rightChild ? _right : _left); }
    int gen() const { return _gen; }
    bool isRoot() const { return _parent == nullptr; }
    bool isLeaf(bool force = false) const { if(!force) return _isLeaf; return _left == nullptr; }

    Node(Node *parent, Node *left, Node *right) : 
      _left(left), _right(right), _parent(parent), _root(parent == nullptr ? this : parent->root()), _gen(parent == nullptr ? 0 : parent->gen() + 1), _isLeaf(left == nullptr) {}

    friend class ElementTree;
};
