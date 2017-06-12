#pragma GCC diagnostic ignored "-Wformat"
#pragma once

#include <cstdio>
#include <string>
#include <cassert>
#include <iostream>
#include <fstream>
#include <ginac/ginac.h>
#include "ginacsupportforeigen.h"
#include "basis.h"

using namespace std;
using namespace GiNaC;

template< typename... Args >
string formatted( const string &format, Args... args ) {
  size_t size = snprintf(nullptr, 0, format.c_str(), args ...);
  unique_ptr<char[]> buf(new char[size]);
  snprintf(buf.get(), size, format.c_str(), args ...);
  return string(buf.get(), buf.get() + size-1);
}

namespace Precomputations {
class Saver {
public:
  static void save(string file, const Eigen::MatrixXd &mat) {
    cout << "Saving to " << file << endl;
    ofstream ofs(file, ofstream::trunc);
    cout << ofs.is_open() << endl;
    if( !ofs) {
      cout << "could not open file for writing" << endl;
      return;
    }
    auto olddigs = Digits;
    Digits = 17;
    size_t r = mat.rows();
    size_t c = mat.cols();
    ofs << r << " " << c << endl;
    for( size_t i = 0; i < r; i++) {
      for( size_t j = 0; j < c; j++) {
        ofs << formatted("%.70f", mat(i,j));
        if( j < c-1) ofs << " ";
      }
      ofs << endl;
    }
    ofs << endl;
    ofs.flush();
    ofs.close();
    Digits = olddigs;
  }

  static void saveSymbolic(string file, const matrix &mat) {
    cout << "Saving symbolic to " << file << endl;
    ofstream ofs(file, ofstream::trunc);
    cout << ofs.is_open() << endl;
    if( !ofs) {
      cout << "could not open file for writing" << endl;
      return;
    }
    size_t r = mat.rows();
    size_t c = mat.cols();
    ofs << r << " " << c << endl;
    for( size_t i = 0; i < r; i++) {
      for( size_t j = 0; j < c; j++) {
        ofs << mat(i,j);
        if( j < c-1) ofs << " ";
      }
      ofs << endl;
    }
    ofs << endl;
    ofs.flush();
    ofs.close();
  }

  static void save(string file, const matrix &mat) {
    cout << "Saving to " << file << endl;
    ofstream ofs(file, ofstream::trunc);
    cout << ofs.is_open() << endl;
    if( !ofs) {
      cout << "could not open file for writing" << endl;
      return;
    }
    auto olddigs = Digits;
    Digits = 70;
    size_t r = mat.rows();
    size_t c = mat.cols();
    ofs << r << " " << c << endl;
    for( size_t i = 0; i < r; i++) {
      for( size_t j = 0; j < c; j++) {
        ofs << mat(i,j).evalf();
        if( j < c-1) ofs << " ";
      }
      ofs << endl;
    }
    ofs << endl;
    ofs.flush();
    ofs.close();
    Digits = olddigs;
  }

  static void save(string file, Basis &b) {
    cout << "Saving basis to " << file << endl;
    ofstream ofs(file, ofstream::trunc);
    cout << ofs.is_open() << endl;
    if( !ofs) {
      cout << "could not open file for writing" << endl;
      return;
    }

    auto olddigs = Digits;
    Digits = 70;

    //print stuff here
    ofs << b.dim() << endl;
    for( unsigned int i = 0; i < b.dim(); i++) {
      symbol x = b.bf(i).x;
      symbol y = b.bf(i).y;
      ex bf = b.bf(i).bf();
      ex pol = bf.expand();
      int deg = max(pol.degree(x), pol.degree(y));
      ofs << deg;
      for( int j = 0; j < deg+1; j++) {
        for( int k = 0; k < deg+1; k++) {
          ofs << " " << pol.coeff(x, j).coeff(y, k).evalf();
        }
      }
      ofs << endl;
    }

    ofs.flush();
    ofs.close();
    Digits = olddigs;
  }
};
}
