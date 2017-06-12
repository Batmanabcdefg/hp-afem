#pragma once
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <vector>
#include <iostream>

#include "math.h"

typedef Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
typedef Eigen::SparseMatrix<scalar>                           SparseMatrix;
typedef Eigen::Matrix<scalar, Eigen::Dynamic, 1>              Vector;
typedef Eigen::SparseVector<scalar>                           SparseVector;
typedef Eigen::Triplet<scalar>                                Trip;
typedef std::vector<Trip>                                     TripVec;

/**
 *  These functions will read a matrix or vector from file.  If the binary flag
 *  is set to true, it will try to read a binary matrix file from filename + ".bin"
 *  and if it can't find that, will create it.
 */
Matrix        readMatrix(std::string filename, bool binary = true);
Vector        readVector(std::string filename, bool binary = true);
void          writeMatrix(std::string filename, Matrix mat, bool binary = true);
void          writeVector(std::string filename, Vector vec, bool binary = true);
