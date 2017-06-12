#include "basis.h"

#include <boost/numeric/ublas/symmetric.hpp>
#include <cstdlib>
#include <fstream>
#include <string>
#include <iostream>

#include "matrix.h"

using namespace std;

Basis::Basis(std::string dirname, int tt, bool bin) {
  string ext = to_string(tt) + ".mat";
  _type = TriType(tt);

  // load the element vector from file
  _eltvec = readVector(dirname + "eltvec_" + ext, bin);
  _dim = _eltvec.rows();

  // load the parts of the element matrix from file
  _eltmat.push_back(readMatrix(dirname + "eltmat0_" + ext, bin));
  _eltmat.push_back(readMatrix(dirname + "eltmat1_" + ext, bin));
  _eltmat.push_back(readMatrix(dirname + "eltmat2_" + ext, bin));
  assert(_eltmat[0].rows() == _dim);
  assert(_eltmat[1].rows() == _dim);
  assert(_eltmat[2].rows() == _dim);

  // load transfer matrices from file
  _transfermat.push_back(readMatrix(dirname + "transfermat0_" + ext, bin));
  _transfermat.push_back(readMatrix(dirname + "transfermat1_" + ext, bin));
  assert(_transfermat[0].rows() == _dim);
  assert(_transfermat[1].rows() == _dim);

  // load mass matrix from file
  _massmat = readMatrix(dirname + "massmat_" + ext, bin);
  assert(_massmat.rows() == _dim);
}

Bases::Bases(string dirname, bool bin) {
  // just load each basis into the vector
  for( auto i = 0; i < 8; i++) {
    _basis.push_back(Basis(dirname, i, bin));
    if( i == 0) _dim = _basis[i]._dim;
    else assert( _basis[i]._dim == _dim);
  }
}
