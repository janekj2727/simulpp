/*
 * A class to deal with matrices of double values in different manners
 * overloading enables vector multiplying, summing, indexing etc.
 * According to chap. 10 of Guide to Scientific Computing in Cpp
 * Author JJ, Date Dec 2019
 * slightly modified Jan 2021
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <cassert>
#include <cstdlib>
#include <vector>

#include "Vector.hpp"
#include "Matrix.hpp"

// Constructor for a matrix of a given size
Matrix::Matrix(int rows, int cols) : mNRows(rows), mNCols(cols)
{
  int i;
  if ((mData = (double *)calloc(mNRows * mNCols, sizeof(double))) == NULL)
  {
    std::cerr << "Warning: Allocation of memory failed!!!\n";
  }

  for (i = 0; i < mNRows * mNCols; i++)
  {
    mData[i] = 0.0;
  }
}

// Constructor for a matrix of a given size with preallocated memory
Matrix::Matrix(int rows, int cols, double *memory) : mNRows(rows), mNCols(cols)
{
  int i;
  mData = memory;

  for (i = 0; i < mNRows * mNCols; i++)
  {
    mData[i] = 0.0;
  }
}

// Overridden copy constructor
Matrix::Matrix(const Matrix &M) : mNRows(M.mNRows), mNCols(M.mNCols)
{
  int i;
  if ((mData = (double *)calloc(mNRows * mNCols, sizeof(double))) == NULL)
  {
    std::cerr << "Warning: Allocation of memory failed!!!\n";
  }

  for (i = 0; i < mNRows * mNCols; i++)
  {
    mData[i] = M.mData[i];
  }
}

// Overridden default destructor
Matrix::~Matrix()
{
  if (mData != nullptr)
  {
    free((void *)mData);
    mData = nullptr;
  }
}

// Public methods to get size
int Matrix::GetNumberOfRows() const
{
  return mNRows;
}

int Matrix::GetNumberOfCols() const
{
  return mNCols;
}

// Overloading operators
Matrix &Matrix::operator=(const Matrix &M)
{
  assert((mNRows == M.mNRows) && (mNCols == M.mNCols));
  int i;
  for (i = 0; i < mNRows * mNCols; i++)
  {
    mData[i] = M.mData[i];
  }
  return *this;
}

Matrix &Matrix::operator+=(const Matrix &M)
{
  assert((mNRows == M.mNRows) && (mNCols == M.mNCols));
  int i;
  for (i = 0; i < mNRows * mNCols; i++)
  {
    mData[i] += M.mData[i];
  }
  return *this;
}

Matrix Matrix::operator+(const Matrix &M) const
{
  assert((mNRows == M.mNRows) && (mNCols == M.mNCols));
  Matrix result(mNRows, mNCols);
  int i;
  for (i = 0; i < mNRows * mNCols; i++)
  {
    result.mData[i] = mData[i] + M.mData[i];
  }
  return result;
}

Matrix Matrix::operator-(const Matrix &M) const
{
  assert((mNRows == M.mNRows) && (mNCols == M.mNCols));
  Matrix result(mNRows, mNCols);
  int i;
  for (i = 0; i < mNRows * mNCols; i++)
  {
    result.mData[i] = mData[i] - M.mData[i];
  }
  return result;
}

Matrix Matrix::operator+() const
{
  Matrix result(mNRows, mNCols);
  int i;
  for (i = 0; i < mNRows * mNCols; i++)
  {
    result.mData[i] = mData[i];
  }
  return result;
}

Matrix Matrix::operator-() const
{
  Matrix result(mNRows, mNCols);
  int i;
  for (i = 0; i < mNRows * mNCols; i++)
  {
    result.mData[i] = -mData[i];
  }
  return result;
}

Matrix Matrix::operator*(const double scalar) const
{
  Matrix result(mNRows, mNCols);
  int i;
  for (i = 0; i < mNRows * mNCols; i++)
  {
    result.mData[i] = mData[i] * scalar;
  }
  return result;
}

// Matrix-vector multiplying
Vector Matrix::operator*(const Vector &v) const
{
  assert(v.GetSize() == mNCols);
  Vector result(mNRows);
  int i, j;
  for (i = 0; i < mNRows; i++)
  {
    for (j = 0; j < mNCols; j++)
    {
      result[i] += v.Read(j) * mData[i * mNCols + j];
    }
  }
  return result;
}

Matrix Matrix::operator*(const Matrix &M) const
{
  assert(mNCols == M.mNRows);
  Matrix result(mNRows, M.mNCols);
  int i, j, k;
  for (i = 0; i < result.mNRows; i++)
  {
    for (j = 0; j < result.mNCols; j++)
    {
      for (k = 0; k < mNCols; k++)
      {
        result.mData[i * M.mNCols + j] += mData[i * mNCols + k] * M.mData[k * M.mNCols + j];
      }
    }
  }
  return result;
}

int Matrix::DeleteMData()
{
  mData = nullptr;
  return 0;
}

// Calculate determinant
double Matrix::CalculateDeterminant() const
{
  assert(mNRows == mNCols);
  double det = 0.0;
  int i, j, k, l;
  if (mNRows == 1)
  {
    return mData[0];
  }
  else
  {
    for (i = 0; i < mNCols; i++)
    {
      Matrix Minor(mNRows - 1, mNCols - 1);
      for (j = 1; j < mNRows; j++)
      {
        l = 0;
        for (k = 0; k < mNCols; k++)
        {
          if (k == i)
          {
            if (++k >= mNCols)
              break;
          }
          Minor(j, l + 1) = mData[j * mNCols + k];
          l++;
        }
      }
      det += pow(-1.0, (double)i) * mData[i] * Minor.CalculateDeterminant();
    }
  }
  return det;
}

// Overloading print << operator
std::ostream &operator<<(std::ostream &output, const Matrix &matrix)
{
  int i, j;
  std::ios_base::fmtflags f(output.flags()); // store initial state of output stream
  // modify stream to print numbers in scientific format with fixed width
  output.setf(std::ios::scientific);
  output.setf(std::ios::showpos);
  output.precision(8);
  output << "/ ";
  for (i = 0; i < matrix.GetNumberOfCols(); i++)
  {
    output << matrix.Read(0, i) << "    ";
  }
  output << "\\\n";
  for (j = 1; j < matrix.GetNumberOfRows() - 1; j++)
  {
    output << "| ";
    for (i = 0; i < matrix.GetNumberOfCols(); i++)
    {
      output << matrix.Read(j, i) << "    ";
    }
    output << "|\n";
  }

  if (matrix.GetNumberOfRows() < 2)
  {
    // nothing
  }
  else
  {
    output << "\\ ";
    for (i = 0; i < matrix.GetNumberOfCols(); i++)
    {
      output << matrix.Read(j, i) << "    ";
    }
    output << "/\n";
  }

  output.flags(f); // restore original state of output stream

  return output;
}

// Scalar product of two rows
double Matrix::RowRowDotProduct(int i, int j) const
{
  double result = 0;
  int k;

  for (k = 0; k < mNCols; k++)
  {
    result += mData[i * mNCols + k] * mData[j * mNCols + k];
  }

  return result;
}

// Sum of one column
double Matrix::SumColumn(int column) const
{
  assert((column < mNCols) && (column >= 0));
  double sum = 0.0;
  int i;
  for (i = 0; i < mNRows; i++)
  {
    sum += mData[i * mNCols + column];
  }
  return sum;
}

// Sum of the first firstN values of the column
double Matrix::SumColumn(int column, int firstN) const
{
  assert((column < mNCols) && (column >= 0));
  double sum = 0.0;
  int i;
  for (i = 0; i < firstN; i++)
  {
    sum += mData[i * mNCols + column];
  }
  return sum;
}

// Print matrix as a markdown table
void Matrix::PrintMDTable(std::ostream &stream, std::vector<std::string> col_head, std::vector<std::string> row_head) const
{
  bool has_row_heads = (row_head.size() > 0) ? true : false;
  int i, j;
  int width, add = 0;
  // int no_cols = mNCols + ((has_row_heads) ? 1 : 0);
  int max_rh_length = 3;
  std::ios_base::fmtflags ff;
  ff = stream.flags();
  std::vector<std::string>::const_iterator it;

  // try to calculate width of columns
  width = stream.precision() + 2;
  if (ff & std::ios_base::scientific)
  {
    add = 5;
  }
  else
  {
    for (i = 0; i < mNRows; i++)
    {
      for (j = 0; j < mNCols; j++)
      {
        if ((mData[i * mNCols + j] > 1e6) || (mData[i * mNCols + j] < 1e-6))
        {
          add = 4;
        }
      }
    }
  }
  width += add;
  if (has_row_heads)
  {
    for (it = row_head.begin(); it != row_head.end(); it++)
    {
      max_rh_length = std::max((int)it->size(), max_rh_length);
    }
    max_rh_length += 4; // to enable ** ** for bold
  }
  if (col_head.size() > 0)
  {
    for (it = col_head.begin(); it != col_head.end(); it++)
    {
      width = std::max((int)it->size(), width);
    }
  }

  // print column headers
  stream << "| ";
  if (has_row_heads)
  {
    stream << std::string(max_rh_length, ' ') << " | ";
  }
  for (i = 0; i < mNCols; i++)
  {
    if (i < (int)col_head.size())
      stream << col_head[i] << std::string(width - col_head[i].size(), ' ') << " | ";
    else
      stream << std::string(width, ' ') << " | ";
  }
  stream << std::endl;
  // print second table row (justify marks)
  stream << "| ";
  if (has_row_heads)
  {
    stream << ":" << std::string(max_rh_length - 1, '-') << " | ";
  }
  for (i = 0; i < mNCols; i++)
  {
    stream << std::string(width - 1, '-') << ": | ";
  }
  stream << std::endl;
  // print table itself
  for (j = 0; j < mNRows; j++)
  {
    stream << "| ";
    if (has_row_heads)
    {
      if (j < (int)row_head.size())
        stream << "**" << row_head[j] << "**" << std::string(max_rh_length - row_head[j].size() - 4, ' ') << " | ";
      else
        stream << std::string(max_rh_length, ' ') << " | ";
    }
    for (i = 0; i < mNCols; i++)
    {
      stream << std::setw(width) << std::right << mData[j * mNCols + i] << " | ";
    }
    stream << std::endl;
  }

  stream.setf(ff);
}

Vector CalculatePartialProduct(const Matrix &A, const Matrix &B, int row)
{
  assert(A.mNCols == B.mNRows);
  Vector result(B.mNCols);
  int i, j;

  for (i = 0; i < result.mSize; i++)
  {
    for (j = 0; j < A.mNCols; j++)
    {
      result.mData[i] += A.mData[row * A.mNCols + j] * B.mData[j * B.mNCols + i];
    }
  }
  return result;
}