/*
* A class to solve linar system having matrix A and rhs vector b
* According to chap. 10 of Guide to Scientific Computing in Cpp
* Author JJ, Date Dec 2019
*/

#include<iostream>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<cassert>
#include<cstdlib>

#include "Vector.hpp"
#include "Matrix.hpp"
#include "LinearSystem.hpp"

// Default (and only) constructor
LinearSystem::LinearSystem(const Matrix& A, const Vector& b)
{
  int size = A.GetNumberOfCols();
  assert((size == A.GetNumberOfRows()) && (size == b.GetSize()));
  mSize = size;
  mpA = new Matrix(A);
  mpb = new Vector(b);
}

// Destructor - free matrix and vector memory
LinearSystem::~LinearSystem()
{
  delete mpA;
  delete mpb;
}

// Solve the system using Gauss-Jordan elimination with pivoting
Vector LinearSystem::SolveGJ() const
{
  Matrix rA(mSize, mSize);
  Vector rb(mSize);
  rA = *mpA;
  rb = *mpb;

  int i, j, k, pivot_row;
  double max, temp;
  Vector result(mSize);

  for (k = 0; k < mSize; k++) {
    // pivoting if necessary
    pivot_row = k;
    max = 0.0;
    for (i = k; i < mSize; i++) {
      pivot_row = (fabs(rA(i,k)) > max) ? i : pivot_row;
    }
    if (pivot_row != k) {
      // swap matrix rows
      for (j = k; j < mSize; j++) {
        temp = rA(pivot_row, j);
        rA(pivot_row, j) = rA(k, j);
        rA(k, j) = temp;
      }
      // swap vector entries
      temp = rb(pivot_row);
      rb(pivot_row) = rb(k);
      rb(k) = temp;
    }
    // zeroing column part below diagonal
    for (i = k + 1; i < mSize; i++) {
      temp = -rA(i, k) / rA(k, k);
      rA(i, k) = 0.0;
      for (j = k + 1; j < mSize; j++) {
        rA(i, j) += temp * rA(k, j);
      }
      rb(i) += temp * rb(k);
    }
    // zeroing column part above diagonal
    for (i = 0; i < k; i++) {
      temp = -rA(i, k) / rA(k, k);
      rA(i, k) = 0.0;
      for (j = k + 1; j < mSize; j++) {
        rA(i, j) += temp * rA(k, j);
      }
      rb(i) += temp * rb(k);
    }
  }
  
  // Now matrix is diagonalized and we can compute the result
  for (k = 0; k < mSize; k++) {
    result(k) = rb(k) / rA(k, k);
  }
  return result;
}

// Solving Linear System using Crammer's rule
Vector LinearSystem::SolveC() const
{
  Matrix rA(mSize, mSize);
  Vector rb(mSize);
  rb = *mpb;
  
  int i, k;
  double det, det_k;

  Vector result(mSize);
  det = mpA->CalculateDeterminant();

  for (k = 0; k < mSize; k++) {
    rA = *mpA;
    for (i = 0; i < mSize; i++) {
      rA(i, k) = rb(i);
    }
    det_k = rA.CalculateDeterminant();
    result(k) = det_k / det;
  }
  return result;
}
