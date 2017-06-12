#include <cstdlib>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <typeinfo>

#include "matrix.h"
#include "print.h"

using namespace std;

namespace Eigen {
  template<class Matrix>

    void write_binary(string filename, const Matrix& matrix){
      std::ofstream out(filename, ios::out | ios::binary | ios::trunc);
      typename Matrix::Index rows=matrix.rows(), cols=matrix.cols();
      out.write((char*) (&rows), sizeof(typename Matrix::Index));
      out.write((char*) (&cols), sizeof(typename Matrix::Index));
      out.write((char*) matrix.data(), rows*cols*sizeof(typename Matrix::Scalar) );
      out.close();
    }

  template<class Matrix>
    void read_binary(string filename, Matrix& matrix){
      cout << "reading binary from " << filename << endl;
      std::ifstream in(filename,ios::in | std::ios::binary);
      assert(in.good());
      typename Matrix::Index rows=0, cols=0;
      in.read((char*) (&rows),sizeof(typename Matrix::Index));
      in.read((char*) (&cols),sizeof(typename Matrix::Index));
      matrix.resize(rows, cols);
      in.read( (char *) matrix.data() , rows*cols*sizeof(typename Matrix::Scalar) );
      in.close();
    }
}

void writeMatrix(string filename, Matrix mat, bool binary) {
  if( binary) return write_binary(filename + scalarname, mat);

  cout << "Saving to " << filename << endl;
  ofstream ofs(filename, ofstream::trunc);
  cout << ofs.is_open() << endl;
  if( !ofs) {
    cout << "could not open file for writing" << endl;
    return;
  }
  int r = mat.rows();
  int c = mat.cols();
  ofs << r << " " << c << endl;
  for( int i = 0; i < r; i++) {
    for( int j = 0; j < c; j++) {
      ofs << Print::formatted("%.70f", mat(i,j));
      if( j < c-1) ofs << " ";
    }
    ofs << endl;
  }
  ofs << endl;
  ofs.flush();
  ofs.close();
}

void writeMatrix(string filename, Vector mat, bool binary) {
  if( binary) return write_binary(filename + scalarname, mat);

  cout << "Saving to " << filename << endl;
  ofstream ofs(filename, ofstream::trunc);
  cout << ofs.is_open() << endl;
  if( !ofs) {
    cout << "could not open file for writing" << endl;
    return;
  }
  int r = mat.rows();
  int c = mat.cols();
  ofs << r << " " << c << endl;
  for( int i = 0; i < r; i++) {
    for( int j = 0; j < c; j++) {
      ofs << Print::formatted("%.70f", mat(i,j));
      if( j < c-1) ofs << " ";
    }
    ofs << endl;
  }
  ofs << endl;
  ofs.flush();
  ofs.close();
}

Matrix readMatrix(string filename, bool binary) {
  string binfn = filename + ".bin"+ scalarname;
  cout << "gonna read " << binfn << endl;
  ifstream binfile(binfn);
  bool isbinfile = binfile.good();
  binfile.close();

  if( binary && isbinfile) {
    Matrix mat;
    Eigen::read_binary(binfn, mat);
    return mat;
  } else {
    ifstream file(filename);
    assert(file.good());

    int rows, cols;
    file >> rows >> cols;

    Matrix mat(rows, cols);
    for( auto i = 0; i < rows; i++) {
      scalar val;
      for( auto j = 0; j < cols; j++) {
        file >> val;
        mat(i,j) = val;
      }
    }

    file.close();
    if( binary) write_binary(binfn, mat);

    return mat;
  }
}

Vector readVector(string filename, bool binary) {
  string binfn = filename + ".bin"+ scalarname;
  ifstream binfile(binfn);
  bool isbinfile = binfile.good();
  binfile.close();

  if( binary && isbinfile) {
    Vector vec;
    Eigen::read_binary(binfn, vec);
    return vec;
  } else {
    ifstream file(filename);
    assert(file.good());

    int rows, cols;
    file >> rows >> cols;
    assert(cols == 1);

    Vector vec(rows);
    scalar val;
    for( auto i = 0; i < rows; i++) {
      file >> val;
      vec(i) = val;
    }

    file.close();
    if(binary) write_binary(binfn, vec);

    return vec;
  }
}
