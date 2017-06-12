#pragma once

#include <vector>
#include <iostream>
#include <string>
#include <array>

#include "hasdof.h"
#include "matrix.h"
#include "tritype.h"

class MatrixCombine;

class Basis : public HasDOF {
 public:
  Basis(std::string dirname, int tt, bool bin = true);

  virtual int dof() const override { return _dim; }
  const TriType &type() const { return _type; }
  const Matrix &eltmat(int i) const { return _eltmat[i]; }
  const Matrix &transfermat(int i) const { return _transfermat[i]; }
  const Matrix &massmat() const { return _massmat; }
  const Vector &eltvec() const { return _eltvec; }
  int p() const { return _dim; }

  friend class Bases;
  friend class MatrixCombine;

 private:
  std::vector<Matrix> _eltmat;
  std::vector<Matrix> _transfermat;
  Matrix _massmat;
  Vector _eltvec;
  int _dim;
  TriType _type;
};

class Bases : public HasDOF {
 public:
  // loads a set of bases from file
  Bases(std::string dirname, bool bin = true);

  const Basis &basis(int tt) { return _basis[tt]; }

  virtual int dof() const override { return _dim; }
  friend class MatrixCombine;

 protected:
  int _dim;
  std::vector<Basis> _basis;
};
