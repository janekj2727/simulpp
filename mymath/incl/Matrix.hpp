#ifndef MATRIXHEADERFILE
#define MATRIXHEADERFILE

/**
 * @class Matrix
 * @ingroup math
 * @brief Implementation of matrix for basic linear algebra operations.
 * Custom implementation of matrices. The class provides multiple two constructors that allow to create a zero matrix of any size or a copy of the existing Matrix.
 * Apart from the basic algebra (addition, multiplication by a scalar, multiplication by a Vector or by a Matrix), determinant calculation, dot product of one line with another or sum of one column is implemented.
 * It was created to be compatible with other classes from @ref math group and is widely used throughout the `simul++` package.
 */
#include <vector>
#include <string>
#include <cassert>
#include "Vector.hpp"

class Matrix
{
private:
  double *mData;    ///< Data (2-dimensional array of doubles saved as one-dimensional) stored in the Matrix.
  const int mNRows; ///< Number of rows of the Matrix
  const int mNCols; ///< Number of columns of the Matrix
  /**
   * @brief Default constructor.
   * Made private in order not to be used because a Matrix of unknown size cannot be instantiated.
   */
  Matrix();

public:
  /**
   * @brief Simple constructor.
   * Creates a new instance of the Matrix class with defined size.
   * The new Matrix contains all zeros.
   * @param rows Number of rows of the new Matrix.
   * @param cols Number of columns of the new Matrix.
   */
  Matrix(int rows, int cols);
  /**
   * @brief Copy constructor.
   * Creates a new instance of the Matrix class that copies M's Size and Data.
   * @param M Template Matrix to be copied.
   */
  Matrix(const Matrix &M);
  /**
   * @brief Matrix constructor with preallocated memory
   * Creates a new instance of the Matrix class with defined size.
   * No memory is allocated for @ref mData, instead a preallocated @p memory is used.
   *
   * @note No checks of preallocated memory size is done, use with great caution!!!
   *
   * @param rows Number of rows of the new Matrix
   * @param cols Number of columns of the new Matrix
   * @param memory Pointer to memory allocated in advance
   */
  Matrix(int rows, int cols, double *memory);
  /**
   * @brief Default destructor.
   * Deallocates the mData array and deletes the Matrix instance.
   */
  ~Matrix();

  /**
   * @brief Get the number of rows of the Matrix.
   * Simple function that provides read-only access to the private member mRows.
   * Once instantiated the size of a Matrix cannot change.
   * @return The number of rows.
   */
  int GetNumberOfRows() const;
  /**
   * @brief Get the number of rows of the Matrix.
   * Simple function that provides read-only access to the private member mRows.
   * Once instantiated the size of a Matrix cannot change.
   * @return The number of rows.
   */
  int GetNumberOfCols() const;
  /**
   * @brief Get the element of the Matrix.
   * The @ref Read(int, int) const function provides zero-based indexing access to read the Matrix element at row \f$row\f$ and column \f$col\f$ (both numbered from zero).
   * It is a save alternative to other element access functions as it does not allow to alter data in the Vector.
   * @param row Row index (zero-based indexing).
   * @param col Column index (zero-based indexing).
   * @return The (value of the) element of the Matrix.
   */
  inline double Read(int row, int col) const
  {
    return mData[row * mNCols + col];
  }; // Read-only zero-based indexing
  /**
   * @brief Get or set the element of the Matrix.
   * The () operator provides zero-based indexing to both read & write the Matrix element at row \f$row\f$ and column \f$col\f$ (both numbered from zero).
   * Thus, it provides access to the private @ref mData.
   * @note Working with pointer to Matrix, one has to use `M->operator()(0,0)` to access or set the element (0,0) of Matrix M, which is admittedly quite impractical.
   * @param row Row index (zero-based indexing).
   * @param col Column index (zero-based indexing).
   * @return The (reference to the) element of the Matrix.
   */
  inline double &operator()(int row, int col)
  {
    assert(row >= 0);
    assert(row < mNRows);
    assert(col >= 0);
    assert(col < mNCols);
    return mData[row * mNCols + col];
  }; // C style zero-based indexing
  /**
   * @brief Assign new data to the Matrix by copying the data of another Matrix.
   * The = operator works as an assignment operator. Thus, the Matrix on the left-hand side is assigned the data of the Matrix resulting from the evaluation of the right-hand side.
   * Both Matrices must have the same size (number of rows and columns, assertion is made inside).
   * @note This is (obviously) not a constructor. Left-hand side Matrix must be instantiated before and have the correct size.
   * @param M The Matrix whose data are to be assigned to the left-hand side Matrix (return value). The source of data.
   * @return The (reference to the) Matrix to which new data are assigned.
   */
  Matrix &operator=(const Matrix &M);
  /**
   * @brief Add two matrices and assign the result to the original matrix.
   * The += operator works as usual. Thus, the Matrix on the left-hand side is assigned the sum with the Matrix resulting from the evaluation of the right-hand side.
   * Both Matrices must have the same size (number of rows and columns, assertion is made inside).
   * @param M The Matrix whose values are to be added to the left-hand side Matrix (return value).
   * @return The (reference to the) Matrix to which the sum is assigned.
   */
  Matrix &operator+=(const Matrix &M);
  /**
   * @brief Binary plus (\f$+\f$) operator (addition).
   * The binary \f$+\f$ operator yields the sum of the two Matrices (the original one to which the method belongs and the other which is specified by the argument of the method).
   * Addition is done as usual, i.e. elementwise.
   * @param M The Matrix to be added to the original Matrix (addend).
   * @return The sum of the two Matrices.
   * @par Example
   * ```cpp
   * Matrix A(2,2), B(2,2), C(2,2); // declare three Matrices of the same size
   * ...                            // do something to alter the data in A and B
   * C = A + B;                     // C now contains the sum of A and B (assignment operator=(const Vector&) also used)
   * ```
   */
  Matrix operator+(const Matrix &M) const;
  /**
   * @brief Binary minus (\f$-\f$) operator (subtraction).
   * The binary \f$-\f$ operator yields the difference of the two Matrices (the original one to which the method belongs and the other which is specified by the argument of the method).
   *
   * @param M The Matrix to be subtracted from the original Matrix (subtrahend).
   * @return The difference of the two Matrices.
   * @par Example
   * ```cpp
   * Matrix A(2,2), B(2,2), C(2,2); // declare three Matrices of the same size
   * ...                            // do something to alter the data in A and B
   * C = A - B;                     // C now contains the difference A - B (assignment operator=(const Vector&) also used)
   * ```
   */
  Matrix operator-(const Matrix &M) const;
  /**
   * @brief Unary plus (\f$+\f$) operator.
   * The unary \f$+\f$ operator is rarely used, but it is defined to make the implementation complete. It returns the copy of the Matrix.
   * For more info, see @ref operator-() const.
   * @return The copy of the Matrix.
   */
  Matrix operator+() const;
  /**
   * @brief Unary minus (\f$-\f$) operator.
   * The unary \f$-\f$ operator is used to invert the sign of the Matrix elements (data). A copy of the original Matrix, but with elements of the opposite sign, is returned.
   *
   * @return Minus the original Matrix.
   */
  Matrix operator-() const;
  /**
   * @brief Multiplication by scalar.
   * The * operator is used to multiply the Matrix by a scalar value.
   *
   * @note Multiplication by zero can be used to set all Matric elements to zero.
   * @param scalar A double-precision number by which the Matrix is multiplied.
   * @return The product of the Matrix-scalar multiplication.
   * @par Example
   * ```cpp
   * Matrix A(3, 3);      // declare a Matrix
   * ...                  // do something to alter the data in a
   * A = A * 10.0         // A is now 10-fold bigger (assignment operator=(const Vector&) also used)
   * ```
   */
  Matrix operator*(const double scalar) const;
  /**
   * @brief Matrix-Vector multiplication.
   * The * operator is used to multiply the Vector by the Matrix, i.e., apply the linear transformation given by the Matrix on the Vector.
   * The result of the multiplication \f$\mathbf{M} \cdot \mathbf{v}\f$ is a Vector \f$\mathbf{a}\f$ whose \f$i\f$-th element is
   * \f[
   *     \sum_{j=0}^{N_{\mathrm{cols}}-1} M_{i,j} \cdot v_j
   * \f]
   * where \f$M_{i,j}\f$ denotes the element of Matrix \f$\mathbf{M}\f$ at \f$i\f$-th row and \f$j\f$-th column and \f$v_j\f$ denotes the \f$j\f$-th element of the Vector \f$\mathbf{v}\f$.
   *
   * @param v The Vector to be multiplied by the Matrix.
   * @return The product (Vector) of the Matrix-Vector multiplication.
   * @par Example
   * ```cpp
   * Matrix A(3, 2);      // declare a Matrix
   * Vector v(2), a(3);   // declare 2 Vectors
   * A(0, 0) = 1;
   * A(0, 1) = 2;
   * A(1, 0) = 3;
   * ...                  // assign values 1, 2, 3,... to Matrix A (A=[[1,2],[3,4],[5,6]])
   * v(0) = -1;
   * v(1) = 1;            // assign values to Vector v (v=[-1,1])
   * a = A * v            // a is now a=[1,1,1] (assignment operator=(const Vector&) also used)
   * ```
   */
  Vector operator*(const Vector &v) const;
  /**
   * @brief Matrix-Matrix multiplication.
   * The * operator is used to multiply the Matrix by another Matrix, i.e., to obtain the composition of two linear maps represented by the Matrices.
   * The result of the multiplication \f$\mathbf{A} \cdot \mathbf{M}\f$ is a Matrix \f$\mathbf{B}\f$ (Matrix product) whose element is
   * \f[
   *     B_{i,j} = \sum_{k=0}^{N-1} A_{i,k} \cdot M_{k,j}
   * \f]
   * where \f$B_{i,j}\f$ denotes the element of Matrix \f$\mathbf{B}\f$ at \f$i\f$-th row and \f$j\f$-th column and \f$N\f$ is the number of columns
   * in Matrix \f$\mathbf{A}\f$ that must be equal to the number of rows in Matrix \f$\mathbf{M}\f$.
   *
   * @param M Matrix to be multiplied by the original Matrix (which the method formally belongs to).
   * @return The product of the Matrix-Matrix multiplication (Matrix-product).
   * @par Example
   * ```cpp
   * Matrix A(3, 2);      // declare Matrices
   * Matrix M(2, 2);
   * Matrix B(3, 2);
   * A(0, 0) = 1;
   * A(0, 1) = 2;
   * A(1, 0) = 3;
   * ...                  // assign values 1, 2, 3,... to Matrix A (A=[[1,2],[3,4],[5,6]])
   * M(0, 0) = 1;
   * M(0, 1) = 2;
   * M(1, 0) = 3;
   * M(1, 1) = 4;         // assign values 1, 2, 3, 4 to Matrix M (M=[[1,2],[3,4]])
   * B = A * M;           // a is now B=[[7,10],[15,22],[23,34]] (assignment operator=(const Matrix&) also used)
   * ```
   */
  Matrix operator*(const Matrix &M) const;
  /**
   * @brief Set @ref mData pointer to `nullptr`
   * Deallocation is performed elsewhere.
   *
   * @return 0 upon success
   */
  int DeleteMData();
  /**
   * @brief Calculate the determinant of the Matrix.
   * Calculate the determinant of the Matrix using the definition (see <a href="https://en.wikipedia.org/wiki/Determinant"> Wiki page </a>) (not suitable for large Matrices).
   * Useful for solving sets of linear equations by the Cramer rule.
   *
   * @return The determinant of the Matrix.
   * @par Example
   * ```cpp
   * Matrix A(2, 2);             // declare a Matrix
   * A(0, 0) = 1;
   * A(0, 1) = 2;
   * A(1, 0) = 3;
   * A(1, 1) = 4;                // assign values 1, 2, 3,... to Matrix A (A=[[1,2],[3,4]])
   * std::cout << A.CalculateDeterminant() << std::endl;
   *                             // prints '-2'
   * ```
   */
  double CalculateDeterminant() const;
  /**
   * @brief Calculate dot product (vector-vector multiply, scalar product) of two Matrix' rows
   * Take row \f$i\f$ of the Matrix as a Vector, row \f$j\f$ as another Vector and calculate the dot product of these two Vectors.
   * The result is a double value calculated by:
   * \f[
   *     \sum_{k=0}{N_{\mathrm{col}}-1} M_{i, k} \cdot M{j, k}
   * \f]
   * where \f$N_{\mathrm{col}}\f$ is the number of columns of Matrix \f$\mathbf{M}\f$.
   *
   * In `simul++` used e.g. to obtain the square of Atom velocity to calculate the kinetic energy.
   *
   * @param i Row-index 1.
   * @param j Row-index 2.
   * @return Dot product of the two Matrix rows.
   */
  double RowRowDotProduct(int i, int j) const;
  /**
   * @brief Sum up one column of the Matrix.
   * Sum up the elements of one column (given by index column).
   *
   * @param column Column-index.
   * @return The sum of the column.
   */
  double SumColumn(int column) const;
  /**
   * @brief Sum up first firstN elements of one column of the Matrix.
   * Sum up the first firstN elements of one column (given by index column).
   *
   * @param column Column-index.
   * @param firstN Only the first firstN terms are added.
   * @return The sum of the first firstN terms in the column.
   */
  double SumColumn(int column, int firstN) const;
  /**
   * @brief Print matrix as a markdown table to @p stream.
   * Headers for columns an rows can be specified by parameters @p col_head and @p row_head.
   * Row headers are **bold** and left-justified, other columns are right-justified.
   * Any modifications (precision, format...) of @p stream must be made prior calling this method and restored after it if desired.
   *
   * @param stream Output stream to which the matrix should be printed.
   * @param col_head Column headers.
   * @param row_head Row headers.
   */
  void PrintMDTable(std::ostream &stream, std::vector<std::string> col_head = std::vector<std::string>(), std::vector<std::string> row_head = std::vector<std::string>()) const;
  friend Matrix CalculateOuterProduct(const Vector &a, const Vector &b);
  friend Vector CalculatePartialProduct(const Matrix &A, const Matrix &B, int row);
};

/**
 * @relates Matrix
 * @brief Overloading the '<<' operator to pretty-print the Matrix @p matrix
 * Pretty-printing of the Matrix @p matrix by `std::cout <<` as:
 * ```
 * / m_0,0      m_0,1     ... \
 * | m_1,0      m_1,1     ... |
 * | ...        ...       ... |
 * \ m_Nrows,0  m_Nrows,1 ... /
 * ```
 *
 * @param output (Reference to) output stream to which the printing is redirected.
 * @param matrix Matrix to be printed.
 * @return The (reference to the) same output stream to enable concatenating (with `<<`).
 *
 * @par Example
 * ```cpp
 * Matrix A(2, 2);             // declare a Matrix
 * A(0, 0) = 1;
 * A(0, 1) = 2;
 * A(1, 0) = 3;
 * A(1, 1) = 4;                // assign values 1, 2, 3,... to Matrix A (A=[[1,2],[3,4]])
 * std::cout << A << std::endl;// print matrix to std::cout
 * ```
 * yields to
 * ```
 * / 1.0000  2.0000 \
 * \ 3.0000  4.0000 /
 * ```
 */
std::ostream &operator<<(std::ostream &output, const Matrix &matrix);
/**
 * @relates Matrix
 * @brief Calculate one row of Matrix-Matrix (dot) product
 *
 * @param A Left Matrix of the product
 * @param B Right Matrix of the product
 * @param row Index of desired row
 * @return Vector containing the @p row -th row of the product @p A . @p B.
 */
Vector CalculatePartialProduct(const Matrix &A, const Matrix &B, int row);

/**
 * @relates Vector
 * @brief Calculate outer product of two Vectors (result is a Matrix).
 *
 * The outer product of two Vectors is given by
 * \f{equation}{
 *     \mathbf{a}\otimes \mathbf{b} =
 *     \mathbf{a} \cdot \mathbf{b}^\mathrm{T} =
 *     \begin{bmatrix}
 *         a_0 b_0 & a_0 b_1 & \cdots \\
 *         a_1 b_0 & a_1 b_1 & \cdots \\
 *         \vdots  & \vdots & \ddots \\
 *     \end{bmatrix}
 * \f}
 *
 * @param a Left Vector of the outer product
 * @param b Right Vector of the outer product
 * @return The outer product of Vectors @p a and @p b.
 */
inline Matrix CalculateOuterProduct(const Vector &a, const Vector &b)
{
  Matrix result(a.mSize, b.mSize);
  int i, j;
  for (i = 0; i < result.mNRows; i++)
  {
    for (j = 0; j < result.mNCols; j++)
    {
      result.mData[i * result.mNCols + j] = a.mData[i] * b.mData[j];
    }
  }
  return result;
}

#endif