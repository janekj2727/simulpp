#ifndef VECTORHEADERFILE
#define VECTORHEADERFILE

#include <iostream>
#include <vector>
#include <cassert>
class Matrix;

/**
 * @defgroup math Mathematics
 * @brief Set of classes and methods for various math operations.
 * The set of classes and functions that provides basic linear algebra operations, equations and equation systems solving (both linear and non-linear), statistics of the data sets etc.
 * Originally based on the book Guide to Scientific Computing in C++ @cite GuideToCpp2012, but largely extended to cover the needs of the `simul++` package.
 *
 * @class Vector
 * @ingroup math
 * @brief Implementation of column vector for basic linear algebra operations.
 * Custom implementation of column vector. The class provides multiple constructors that allow to create a zero vector of any size, the copy of the existing Vector and even to read the vector data from a file.
 * Apart from the basic linear algebra operations, various statistics operations are provided including various estimators of the standard error of mean.
 * It was created to be compatible with other classes from @ref math group and is widely used throughout the `simul++` package.
 */

class Vector
{
private:
  double *mData; ///< Data (array of doubles) stored in the Vector.
  int mSize;     ///< Size (length) of the Vector. All Vectors are assumed to be column vectors.
  /**
   * @brief Default constructor.
   * Made private in order not to be used because a Vector of unknown length cannot be instantiated.
   */
  Vector();

public:
  /**
   * @brief Copy constructor.
   * Creates a new instance of the Vector class that copies otherVector's Size and Data.
   * @param otherVector Template Vector to be copied.
   */
  Vector(const Vector &otherVector);
  /**
   * @brief Simple constructor.
   * Creates a new instance of the Vector class with defined length.
   * The new Vector contains all zeros.
   * @param size Length (mSize) of the new Vector.
   */
  Vector(int size);
  /**
   * @brief Constructor that reads the Vector from file.
   * Creates a new instance of the Vector class with data stored in a file.
   * Each line of the file should contain one floating-point number.
   * The size is inferenced automatically from the number of lines in the file.
   * The file can contain comment lines starting by `#`, `!` or `%` which are neither processed nor counted.
   * @param fin C-style file stream opened to read.
   */
  Vector(FILE *fin); // initialize Vector by reading file..
  /**
   * @brief Constructor that reads the Vector from file with given name.
   * Creates a new instance of the Vector class with data stored in a given file.
   * Similar to @ref Vector(FILE*) constructor but includes the opening and closing of the file.
   * Thus, only the file name is sufficient.
   *
   * Each line of the file should contain one floating-point number.
   * The size is inferenced automatically from the number of lines in the file.
   * The file can contain comment lines starting by `#`, `!` or `%` which are neither processed nor counted.
   * @param name_of_file The name of the file that stores the data.
   */
  Vector(const std::string name_of_file); // initialize Vector from file... (by name)
  /**
   * @brief Construct a new Vector object from `std::vector<double>`
   * Creates a new instance of the Vector class with data stored in `std::vector<double>` object.
   * Conversion from standard vector to Vector class.
   *
   * @param orig_vec original vector to be copied.
   */
  Vector(const std::vector<double> orig_vec); // initialize Vector from std::vector<double>
  /**
   * @brief Default destructor.
   * Deallocates the mData array and deletes the Vector instance.
   */
  ~Vector();

  /**
   * @brief Get the size (length) of the Vector.
   * Simple function that provides read-only access to the private member mSize.
   * @return The length of the Vector.
   */
  int GetSize() const;
  /**
   * @brief Get or set the \f$i\f$-th element of the Vector.
   * The [] operator provides zero-based indexing to both read & write the Vector element.
   * Thus, it works like a shortcut to provide access to the private @ref mData.
   * Mimics the [] operator of a basic array in C++ (works as `mData[i]`).
   * @param i Element index (zero-based indexing).
   * @return The (reference to) \f$i\f$-th element of the Vector.
   */
  double &operator[](int i); // zero-based indexing
  /**
   * @brief Get the \f$i\f$-th element of the Vector.
   * The @ref Read(int) const function provides zero-based indexing access to read the Vector element.
   * It is a save alternative to other element access functions as it does not allow to alter data in the Vector.
   * @param i Element index (zero-based indexing).
   * @return The (value of) \f$i\f$-th element of the Vector.
   */
  inline double Read(int i) const
  {
    assert((i > -1) && (i < mSize));
    return mData[i];
  }; // read-only zero-based indexing
  /**
   * @brief Get or set the \f$i\f$-th element of the Vector.
   * The () operator provides zero-based indexing to both read & write the Vector element.
   * Thus, it works exactly the same as @ref operator[](int).
   *
   * @note Originally, in @cite GuideToCpp2012 it provided Matlab-style one-based indexing, but this behaviour was changed to prevent confusion.
   * Therefore, all indexing in `simul++` package is zero-based (C-style).
   * @param i Element index (zero-based indexing).
   * @return The (reference to) \f$i\f$-th element of the Vector.
   */
  inline double &operator()(int i)
  {
    assert((i > -1) && (i < mSize));
    return mData[i];
  }; // zero-based indexing
  /**
   * @brief Assign new data to the Vector by copying the data of another Vector.
   * The = operator works as an assignment operator. Thus, the Vector on the left-hand side is assigned the data of the Vector resulting from the evaluation of the right-hand side.
   * Both Vectors must have the same length (assertion is made inside).
   * @par Example
   * ```cpp
   * Vector a(2), b(2);       // declare two Vectors of the same size (2)
   * a(0) = 1.0, a(1) = 2.0;  // assign data to a (a is now [1.0, 2.0])
   * b = a;                   // b now contains the same elements (data) as a
   * ```
   *
   * @note This is (obviously) not a constructor. Left-hand side Vector must be instantiated before and have the correct length.
   * @param otherVector The Vector whose data are to be assigned to the left-hand side Vector (return value). The source of data.
   * @return The (reference to) Vector to which new data are assigned.
   */
  Vector &operator=(const Vector &otherVector);
  /**
   * @brief Unary plus (\f$+\f$) operator.
   * The unary \f$+\f$ operator is rarely used, but it is defined to make the implementation complete. It returns the copy of the Vector.
   * For more info, see @ref operator-().
   *
   * @return The (copy of) original Vector.
   */
  Vector operator+() const; // unary +
  /**
   * @brief Unary minus (\f$-\f$) operator.
   * The unary \f$-\f$ operator is used to invert the sign of the Vector elements (data). A copy of the original Vector, but with elements of the opposite sign, is returned.
   *
   * @return Minus the original Vector.
   */
  Vector operator-() const; // unary -
  /**
   * @brief Binary plus (\f$+\f$) operator (addition).
   * The binary \f$+\f$ operator yields the sum of the two Vectors (the original one to which the method belongs and the other which is specified by the argument of the method).
   *
   * @param otherVector The Vector to be added to the original Vector (addend).
   * @return The sum of the two Vectors.
   * @par Example
   * ```cpp
   * Vector a(3), b(3), c(3); // declare three Vectors of the same size
   * ...                      // do something to alter the data in a and b
   * c = a + b;               // c now contains the sum of a and b (assignment @ref operator=(const Vector&) also used)
   * ```
   */
  Vector operator+(const Vector &otherVector) const; // binary +
  /**
   * @brief Binary minus (\f$-\f$) operator (subtraction).
   * The binary \f$-\f$ operator yields the difference of the two Vectors (the original one to which the method belongs and the other which is specified by the argument of the method).
   *
   * @param otherVector The Vector to be subtracted from the original Vector (subtrahend).
   * @return The difference of the two Vectors.
   * @par Example
   * ```cpp
   * Vector a(3), b(3), c(3); // declare three Vectors of the same size
   * ...                      // do something to alter the data in a and b
   * c = a - b;               // c now contains the difference a - b (assignment @ref operator=(const Vector&) also used)
   * ```
   */
  Vector operator-(const Vector &otherVector) const; // binary -
  /**
   * @brief Multiplication by scalar.
   * The * operator is used to multiply the Vector by a scalar value.
   *
   * @note Multiplication by zero can be used to set all Vector values to zero.
   * @param scalar A double-precision number by which the Vector is multiplied.
   * @return The product of the Vector-scalar multiplication.
   * @par Example
   * ```cpp
   * Vector a(3);         // declare a Vector
   * ...                  // do something to alter the data in a
   * a = a * 10.0         // a is now 10-fold bigger
   * ```
   */
  Vector operator*(const double scalar) const;
  /**
   * @brief Multiplication by scalar in place.
   * The *= operator is used to multiply the Vector by a scalar value in place.
   *
   * @note Multiplication by zero can be used to set all Vector values to zero (see @ref Clear()).
   * @param scalar A double-precision number by which the Vector is multiplied.
   * @return The original vector multiplied by the scalar.
   * @par Example
   * ```cpp
   * Vector a(3);         // declare a Vector
   * ...                  // do something to alter the data in a
   * a *= 10.0            // a is now 10-fold bigger
   * ```
   */
  Vector &operator*=(const double scalar);
  /**
   * @brief Clear vector (assign 0.0 to all indeces)
   *
   * @return 0
   */
  int Clear();
  /**
   * @brief Calculate the \f$p\f$-norm of the Vector.
   * The method returns the \f$p\f$-norm of the Vector given by (for Vector \f$\mathbf{a}\f$):
   * \f{equation}{
   * \left|\left| \mathbf{a} \right|\right|_p = \left(\sum_{i=0}^{\mathtt{mSize}-1} a_i^p \right)^{(1/p)}
   * \f}
   * The standard Euclidean norm is recovered by using the default version (\f$p=2\f$).
   * Using \f$p=1\f$, the function behaves as @ref Sum().
   * The maximum norm can theoretically be reached by using some large value of \f$p\f$, but this cannot be recommended.
   *
   * @param p The norm specification (default: 2 (the Euclidean norm)).
   * @return The norm of the Vector.
   * @par Example
   * ```cpp
   * Vector a(2);            // declare a Vector of size 2 (mSize=2)
   * a(0) = 3.0, a(1) = 4.0; // a is now [3; 4]
   * a.CalculateNorm();      // returns 5.0
   * a.CalculateNorm(1);     // returns 7.0
   * ```
   */
  double CalculateNorm(int p = 2) const; // p-norm of a vector
  /**
   * @brief Print Vector to file.
   * Print the Vector to a file specified by the C-like stream, each Vector element on a new line.
   * Useful to save Vector for future use. The resulted file can be loaded by one of the two *read-file* constructors
   * @ref Vector(FILE*) or @ref Vector(const std::string).
   *
   * Similar to @ref PrintToFile(const std::string) const, but does not handle the opening and closing of the file.
   *
   * @param fout The pointer to the C-like output stream to which the Vector is to be saved.
   * @return
   * -  0 upon success
   * - -1 if output stream not defined (`fout == NULL`)
   */
  int PrintToFile(FILE *fout) const; // print vector to a file
  /**
   * @brief Print Vector to file.
   * Print the Vector to a file specified by the file name, each Vector element on a new line.
   * Useful to save Vector for future use. The resulted file can be loaded by one of the two *read-file* constructors
   * @ref Vector(FILE*) or @ref Vector(const std::string).
   *
   * Similar to @ref PrintToFile(FILE*) const, but opens the file first, then reads it and then close it again.
   * Moreover, the time stamp is added as a comment on the first line of the file.
   *
   * @param filename Name of the file where the Vector is to be saved.
   * @return
   * -  0 upon success
   * - -1 if the file cannot be opened
   */
  int PrintToFile(const std::string filename) const; // print vector to a file
  /**
   * @brief Calculate the mean value of the Vector elements.
   * Calculate the arithmetic mean (average) of the Vector \f$\mathbf{v}\f$ elements using:
   * \f{equation}{
   *     \frac{1}{N}\sum_{i=0}^{N-1} v_i
   * \f}
   * where \f$N=\f$ @ref mSize.
   *
   * @par See also
   * @ref MeanWithError(double&, bool) const, @ref MeanWithErrorEstimate(double&) const, @ref BlockMean(double&, int) const, @ref BlockMeanFlyvbjerg(double&, double&, int) const.
   *
   * @return Arithmetic mean.
   */
  double Mean() const; // calculate the mean value of a vector
  /**
   * @brief Calculate the mean value of the Vector elements and its error.
   * Calculate the arithmetic mean (average) of the Vector \f$\mathbf{v}\f$ elements using:
   * \f{equation}{
   *     \overline{v}=\frac{1}{N}\sum_{i=0}^{N-1} v_i
   * \f}
   * where \f$N=\f$ @ref mSize.
   *
   * The standard error of mean is calculated by (assuming uncorrelated data, uncorrelated Vector elements):
   * \f{equation}{
   *     \sigma_{\overline{v}}=\sqrt{\frac{\sum_{i=0}^{N-1}(v_i-\overline{v})^2}{N M}}=
   *     \sqrt{\frac{\sum_{i=0}^{N-1}v_i^2 - \left(\sum_{i=0}^{N-1}v_i\right)^2}{N M}}
   * \f}
   * where \f$M=N-1\f$ if Bessel corrected (default), \f$M=N\f$ otherwise.
   *
   * @par See also
   * @ref Mean() const, @ref MeanWithErrorEstimate(double&) const, @ref BlockMean(double&, int) const, @ref BlockMeanFlyvbjerg(double&, double&, int) const.
   *
   * @param sigma Standard error of the mean (intended as output, uncorrelated data).
   * @param Bessel Use the Bessel correction to calculate the error (default: true).
   * @return Arithmetic mean.
   */
  double MeanWithError(double &sigma, bool Bessel = true) const; // calculate mean and its standard error (with Bessel correction)
  /**
   * @brief Calculate the mean value of the Vector elements (correlated) and its error.
   * Calculate the arithmetic mean (average) of the Vector \f$\mathbf{v}\f$ elements using:
   * \f{equation}{
   *     \overline{v}=\frac{1}{N}\sum_{i=0}^{N-1} v_i
   * \f}
   * where \f$N=\f$ @ref mSize.
   *
   * The easiest way to estimate the standard error of mean for correlated data makes use of cumulative running avarage @cite kolafaMolModSim p 71:
   * \f{equation}{
   *     \overline{\overline{v}}_n=\frac{1}{n}\sum_{i=0}^{n-1} v_i
   * \f}
   * where \f$n\f$ `runs' from 0 to @ref mSize. To estimate the error of mean for correlated data, the second half of the sequence \f$\overline{\overline{v}}_n\f$
   * is used and the error calculated by taking the difference between the maximum and minimum of it and multiply it by 0.6:
   * \f{equation}{
   *     \sigma_{\overline{v}} \approx 0.6 \left(\max_{i \in [N/2, N-1]}\overline{\overline{v}}_i-\min_{i \in [N/2, N-1]}\overline{\overline{v}}_i\right)
   * \f}
   *
   * @par See also
   * @ref Mean() const, @ref MeanWithError(double&, bool) const, @ref BlockMean(double&, int) const, @ref BlockMeanFlyvbjerg(double&, double&, int) const.
   *
   * @param sigma Standard error of the mean (intended as output, correlated data).
   * @return Arithmetic mean.
   */
  double MeanWithErrorEstimate(double &sigma) const; // calculate mean of correlated data and estimate its stderr (roughly) using the running avarage method
  /**
   * @brief Calculate the mean value of the Vector elements and variance of block means.
   * Calculate the arithmetic mean (average) of the Vector \f$\mathbf{v}\f$ elements using:
   * \f{equation}{
   *     \overline{v}=\frac{1}{N}\sum_{i=0}^{N-1} v_i
   * \f}
   * where \f$N=\f$ @ref mSize.
   *
   * One of the methods how to calculate the error of mean (of correlated data) is the blocking method @cite FrenkelUnderstanding (Appendix D) @cite AllenCompSimLiq (p 172 (pdf 196)).
   * The Vector is divided into @p Nblocks subsequent blocks and mean values are calculated from each block:
   * \f{equation}{
   *     \overline{v}_j=\frac{1}{N_\mathrm{b}}\sum_{i=j\cdot N_b}^{(j+1)\cdot N_b-1} v_i
   * \f}
   * where \f$N_\mathrm{b}=N/\f$ @p Nblocks (rounded to integer) is the length of one block.
   * Having thus @p Nblocks values of mean estimate (\f$\overline{v}_j\f$), we can calculate their variance:
   * \f{equation}{
   *     \sigma^2_{\overline{v}_j} = \frac{1}{\mathtt{Nblocks}}\sum_{j=0}^{\mathtt{Nblocks}-1}(\overline{v}_j - \overline{v})^2
   * \f}
   * which is assigned to @p variance by this function.
   *
   * To complete the error analysis, the value of characteristic decay time (\f$t_v^\mathrm{c}\f$) @cite FrenkelUnderstanding (or statistic efficiency @cite AllenCompSimLiq)
   * can be estimated as a limit:
   * \f{equation}{
   *     t_v^\mathrm{c}=\lim_{N_\mathrm{b} \rightarrow \infty} \frac{N_\mathrm{b} \sigma^2_{\overline{v}_j}}{\sigma^2_v}
   * \f}
   * where \f$\sigma^2_v\f$ is the (uncorrected) variance of all data (see @ref Variance(bool) const).
   * The estimate of standard error of mean (of correlated data) is then:
   * \f{equation}{
   *     \sigma_{\overline{v}} \approx \sqrt{\sigma^2_v \frac{t_v^\mathrm{c}}{N}}
   * \f}
   * To perform this task, the script `interactive_statistics.mms` can be used.
   *
   * @par See also
   * @ref Mean() const, @ref MeanWithError(double&, bool) const, @ref MeanWithErrorEstimate(double&) const, @ref BlockMeanFlyvbjerg(double&, double&, int) const,
   * @ref Variance(bool) const.
   *
   * @note Methods @ref BlockMean(double&, int) const and @ref BlockMeanFlyvbjerg(double&, double&, int) const assign @b different things to @p variance !
   *
   * @param variance Variance of the means calculated from @p Nblocks values (\f$\sigma^2_{\overline{v}_j}\f$) (intended as output).
   * @param Nblocks Number of data blocks.
   * @return Arithmetic mean of data.
   */
  double BlockMean(double &variance, int Nblocks) const; // calculate mean and blocked averages variance via blocking method with Nblocks blocks
  /**
   * @brief Calculate the mean value of the Vector elements and estimates variance of the mean and its error using Flyvbjerg blocking method.
   * Calculate the arithmetic mean (average) of the Vector \f$\mathbf{v}\f$ elements using:
   * \f{equation}{
   *     \overline{v}=\frac{1}{N}\sum_{i=0}^{N-1} v_i
   * \f}
   * where \f$N=\f$ @ref mSize.
   *
   * One of the methods how to calculate the error of mean (of correlated data) is the Flyvbjerg blocking method @cite FrenkelUnderstanding (Appendix D) (see also @ref BlockMean(double&, int) const).
   * The Vector is divided into @p Nblocks subsequent blocks and mean values are calculated from each block:
   * \f{equation}{
   *     \overline{v}_j=\frac{1}{N_\mathrm{b}}\sum_{i=j\cdot N_b}^{(j+1)\cdot N_b-1} v_i
   * \f}
   * where \f$N_\mathrm{b}=N/\f$ @p Nblocks (rounded to integer).
   * Having thus @p Nblocks values of mean estimate (\f$\overline{v}_j\f$), we can calculate:
   * \f{equation}{
   *     \sigma^2_{\overline{v}_j} = \frac{1}{\mathtt{Nblocks}}\sum_{i=0}^{\mathtt{Nblocks}-1}(\overline{v}_j - \overline{v})^2.
   * \f}
   * If the data were uncorrelated (or the blocks sufficiently long), the ratio:
   * \f{equation}{
   *     \frac{\sigma^2_{\overline{v}_j}}{\mathtt{Nblocks}-1}
   * \f}
   * (this value is assigned to @p variance ) would be constant and equal to the variance of mean of original data.
   * Therefore, to estimate the error of mean, sufficiently long blocks should be used to ensure that consequent blocks are uncorrelated.
   * Nevertheless, the longer blocks, the less number of them and the error of the estimate larger.
   * The error of the estimated variance of mean can be estimated by:
   * \f{equation}{
   *     \sqrt{\frac{2 \sigma^4_{\overline{v}_j}}{(\mathtt{Nblocks}-1)^3}}
   * \f}
   * (this value is assigned to @p varianceError ).
   * The best estimate of \f$\sigma^2_{\overline{v}}\f$ is thus obtained by calculating its value together with its error for multiple block lengths
   * (multiple number of blocks) and to pick the @e plateau value
   * (values have already converged, but the error is still quite small).
   * To perform this task, the script `interactive_statistics.mms` can be used.
   *
   * @par See also
   * @ref Mean() const, @ref MeanWithError(double&, bool) const, @ref MeanWithErrorEstimate(double&) const, @ref BlockMeanFlyvbjerg(double&, double&, int) const,
   * @ref Variance(bool) const.
   *
   * @note Methods @ref BlockMean(double&, int) const and @ref BlockMeanFlyvbjerg(double&, double&, int) const assign @b different things to @p variance !
   *
   * @param variance Estimate of the variance of mean \f$\sigma^2_{\overline{v}}\f$.
   * @param varianceError Estimated error of the @p variance above.
   * @param Nblocks Number of data blocks.
   * @return Arithmetic mean of data.
   */
  double BlockMeanFlyvbjerg(double &variance, double &varianceError, int Nblocks) const; // calculate mean and estimate of variance of ensemble average
  /**
   * @brief Calculate the variance of the Vector elements.
   * Calculate the variance of data stored in the Vector by formula:
   * \f{equation}{
   *     \sigma_v^2=\frac{1}{N} \sum_{i=0}^{\mathtt{mSize-1}} \left(v_i - \overline{v}\right)^2
   * \f}
   * where v_i is the \f$i\f$-th element of the Vector (numbering from zero), \f$\overline{v}\f$ is the mean value of the Vector elements (see @ref Mean() const) and
   * \f$N=\f$ @ref mSize if @p Bessel = false (default) or \f$N=\f$ @ref mSize \f$-1\f$ if Bessel correction demanded ( @p Bessel = true).
   *
   * @par See also
   * @ref Mean() const, @ref Variance(double, bool) const, @ref SampleVariance(bool) const.
   *
   * @note Same as @ref SampleVariance(bool) const but with the different default (not using Bessel correction).
   *
   * @param Bessel Use Bessel correction (default false)
   * @return Variance of the Vector elements (\f$\sigma_v^2\f$).
   */
  double Variance(bool Bessel = false) const; // calculate the variance of the vector
  /**
   * @brief Calculate the variance of the Vector elements, use the given mean value.
   * Calculate the variance of data stored in the Vector by formula:
   * \f{equation}{
   *     \sigma_v^2=\frac{1}{N} \sum_{i=0}^{\mathtt{mSize-1}} \left(v_i - \overline{v}_\mathrm{given}\right)^2
   * \f}
   * where v_i is the \f$i\f$-th element of the Vector (numbering from zero), \f$\overline{v}_\mathrm{given}\f$ is the @b given @p mean value and
   * \f$N=\f$ @ref mSize if @p Bessel = false (default) or \f$N=\f$ @ref mSize \f$-1\f$ if Bessel correction demanded ( @p Bessel = true).
   *
   * @par See also
   * @ref Variance(bool) const, @ref VarianceWithError(double&, int) const, @ref SampleVariance(bool) const.
   *
   * @param mean Mean value used in calculation (\f$\overline{v}_\mathrm{given}\f$).
   * @param Bessel Use Bessel correction (default false)
   * @return Variance of the Vector elements (\f$\sigma_v^2\f$).
   */
  double Variance(double mean, bool Bessel = true) const; // calculate the variance when mean is known
  /**
   * @brief Calculate the sample variance of the Vector elements.
   * Calculate the variance of data stored in the Vector by formula:
   * \f{equation}{
   *     \sigma_v^2=\frac{1}{N} \sum_{i=0}^{\mathtt{mSize-1}} \left(v_i - \overline{v}\right)^2
   * \f}
   * where v_i is the \f$i\f$-th element of the Vector (numbering from zero), \f$\overline{v}\f$ is the mean value of the Vector elements (see @ref Mean() const) and
   * \f$N=\f$ @ref mSize \f$-1\f$ if @p Bessel = true (default) or \f$N=\f$ @ref mSize if @p Bessel = true.
   *
   * @par See also
   * @ref Mean() const, @ref Variance(double, bool) const, @ref Variance(bool) const.
   *
   * @note Same as @ref Variance(bool) const but with the different default (using Bessel correction).
   *
   * @param Bessel Use Bessel correction (default false)
   * @return Sample variance of the Vector elements (\f$\sigma_v^2\f$).
   */
  double SampleVariance(bool Bessel = true) const; // calculate the sample variance (Bessel corrected)
  /**
   * @brief Calculate the variance of the Vector elements and estimate its error by simple blocking method with given number of blocks.
   * Calculate the mean value of the Vector elements:
   * \f{equation}{
   *     \overline{v}=\frac{1}{N}\sum_{i=0}^{N-1} v_i
   * \f}
   * where \f$N=\f$ @ref mSize.
   *
   * To calculate the variance of the data stored in the Vector and the error of this variance, simple blocking method is used (assuming uncorrelated blocks).
   * The Vector is divided into @p Nblocks subsequent blocks and variances are calculated from each block:
   * \f{equation}{
   *     \sigma^2_j=\frac{1}{N_\mathrm{b}}\sum_{i=j\cdot N_\mathrm{b}}^{(j+1)\cdot N_\mathrm{b}-1} \left(v_i - \overline{v}\right)^2
   * \f}
   * where \f$N_\mathrm{b}=N/\f$ @p Nblocks (rounded to integer).
   * Note that the overall mean value is used (not the mean value of the block, see @ref Variance(double, bool) const.
   * Having thus @p Nblocks estimates of variance, we can calculate the mean variance
   * \f{equation}{
   *     \sigma_v^2 = \frac{1}{\mathtt{Nblocks}}\sum_{j=0}^{\mathtt{Nblocks}-1}\sigma^2_j
   * \f}
   * and estimate its error by:
   * \f{equation}{
   *     \sigma_{\sigma^2_v} = \sqrt{\frac{\sum_{j=0}^{\mathtt{Nblocks}-1}(\sigma^2_j-\sigma^2_v)}{\mathtt{Nblocks}(\mathtt{Nblocks}-1)}}
   * \f}
   *
   * @par See also
   * @ref Variance(bool) const, @ref MeanWithError(double&, bool) const, @ref SampleVariance(bool) const.
   *
   * @param sigma Estimated standard error of the data variance (\f$\sigma_{\sigma^2_v}\f$, intended as output).
   * @param Nblocks Number of blocks for blocking method.
   * @return Variance of the Vector elements (\f$\sigma_v^2\f$).
   */
  double VarianceWithError(double &sigma, int Nblocks) const; // calculate variance and estimate the stderr of variance using Nblocks blocks
  /**
   * @brief Calculate excess kurtosis of the data stored in the Vector and estimate its error by blocking method.
   * To calculate the excess kurtosis and estimate its error, simple blocking method is used (assuming uncorrelated subsequent blocks).
   * For each block the kurtosis is calculated by one of the two following formulae.
   * The uncorrected (for default @p Bessel = false) excess kurtosis is calculated by:
   * \f{equation}{
   *     \mathrm{Kurt}(v)_j=\frac{\frac{1}{\mathtt{N_\mathrm{b}}} \sum_{i=j\cdot N_\mathrm{b}}^{(j+1)\cdot N_\mathrm{b}-1}
   *     (v_i - \overline{v})^4}{\sigma_v,j^4} - 3 = k_{v,j}
   * \f}
   * where \f$N_\mathrm{b}=N/\f$ @p Nblocks (rounded to integer), \f$v_i\f$ is the \f$i\f$-th element of the Vector (numbering from zero), \f$\overline{v}\f$ is the mean value of the data (global mean value, see @ref Mean() const)
   * and \f$\sigma^4_v = \left(\sigma_v,j^2\right)^2\f$ is the squared variance of the block (Vector elements at indices from \f$j\cdot N_\mathrm{b}\f$ to \f$(j+1)\cdot N_\mathrm{b}-1\f$) (uncorrected, see @ref Variance(bool) const).
   * If Bessel correction demanded ( @p Bessel = true), the corrected value is calculated by (see @cite Harding2014 and <a href="https://en.wikipedia.org/wiki/Kurtosis"> Wikipedia article </a>):
   * \f{equation}{
   *     \mathrm{Kurt}(v)_j = \frac{n-1}{(n-2)(n-3)}\left[(n+1)k_{v,j} + 6\right] = \frac{(n+1)n(n-1)}{(n-2)(n-3)}
   *     \frac{\sum_{i=j\cdot N_\mathrm{b}}^{(j+1)\cdot N_\mathrm{b}-1} (v_i - \overline{v})^4}
   *     {\left(\sum_{i=j\cdot N_\mathrm{b}}^{(j+1)\cdot N_\mathrm{b}-1} (v_i - \overline{v})^2\right)^2} - 3\frac{(n-1)^2}{(n-2)(n-3)}
   * \f}
   * Note that the overall mean value is used (not the mean value of the block, see @ref Variance(double, bool) const.
   *
   * Having thus @p Nblocks estimates of kurtosis, we can calculate the mean kurtosis
   * \f{equation}{
   *     \mathrm{Kurt}(v) = \frac{1}{\mathtt{Nblocks}}\sum_{j=0}^{\mathtt{Nblocks}-1}\mathrm{Kurt}(v)_j
   * \f}
   * and estimate its error by:
   * \f{equation}{
   *     \sigma_{\mathrm{Kurt}(v)} = \sqrt{\frac{\sum_{j=0}^{\mathtt{Nblocks}-1}(\mathrm{Kurt}(v)_j-\mathrm{Kurt}(v))}{\mathtt{Nblocks}(\mathtt{Nblocks}-1)}}
   * \f}
   *
   * @par See also
   * @ref ExcessKurtosis(bool) const, @ref Variance(bool) const, @ref VarianceWithError(double&, int) const, @ref Mean() const.
   *
   * @note Excess kurtosis is used here (3 subtracted), thus the value should be zero for normally distributed data.
   * Bessel correction is always used for estimating the error from @p Nblocks values of kutosis.
   * The parameter @p Bessel determines only is the Bessel correction is used within calculation of each individual block kurtoses.
   *
   * @param sigma Estimated error of excess kurtosis (\f$\sigma_{\mathrm{Kurt}(v)}\f$, intended as output).
   * @param Nblocks Number of blocks for blocking method.
   * @param Bessel Use Bessel correction for block kurtoses (default false – do not use correction).
   * @return Estimated excess kurtosis of the data stored in the Vector (\f$\mathrm{Kurt}(v)\f$).
   */
  double ExcessKurtosisWithError(double &sigma, int Nblocks, bool Bessel = false) const; // calculate excess kurtosis and estimate the error by blocking method
  /**
   * @brief Calculate excess kurtosis of the data stored in the Vector.
   * The uncorrected (for default @p Bessel = false) excess kurtosis is calculated by:
   * \f{equation}{
   *     \mathrm{Kurt}(v)=\frac{\frac{1}{\mathtt{mSize}} \sum_{i=0}^{\mathtt{mSize}-1} (v_i - \overline{v})^4}{\sigma_v^4} - 3 = k_v
   * \f}
   * where \f$v_i\f$ is the \f$i\f$-th element of the Vector (numbering from zero), \f$\overline{v}\f$ is the mean value of the data (see @ref Mean() const)
   * and \f$\sigma^4_v = \left(\sigma_v^2\right)^2\f$ is the squared variance of the data (uncorrected, see @ref Variance(bool) const).
   *
   * If Bessel correction demanded ( @p Bessel = true), the corrected value is calculated by (see @cite Harding2014 and <a href="https://en.wikipedia.org/wiki/Kurtosis"> Wikipedia article </a>):
   * \f{equation}{
   *     \mathrm{Kurt}(v) = \frac{n-1}{(n-2)(n-3)}\left[(n+1)k_v + 6\right] = \frac{(n+1)n(n-1)}{(n-2)(n-3)} \frac{\sum_{i=0}^{\mathtt{mSize}-1} (v_i - \overline{v})^4}{\left(\sum_{i=0}^{\mathtt{mSize}-1} (v_i - \overline{v})^2\right)^2} - 3\frac{(n-1)^2}{(n-2)(n-3)}
   * \f}
   *
   * @par See also
   * @ref ExcessKurtosisWithError(double&, int, bool) const, @ref Variance(bool) const, @ref Mean() const.
   *
   * @note Excess kurtosis is used here (3 subtracted), thus the value should be zero for normally distributed data.
   *
   * @param Bessel Use Bessel correction (default false – do not use correction).
   * @return Estimated excess kurtosis of the data stored in the Vector (\f$\mathrm{Kurt}(v)\f$).
   */
  double ExcessKurtosis(bool Bessel = false) const; // calculate excess kurtosis
  /**
   * @brief Calculate the sum of Vector entries.
   * Sum up the elements of the Vector.
   *
   * @par See also
   * @ref Sum(int) const.
   *
   * @return Sum of Vector entries.
   */
  double Sum() const; // calculate sum of vector entries
  /**
   * @brief Calculate the sum of the first @p firstN Vector entries.
   * Sum up the first @p firstN elements of the Vector.
   *
   * @par See also
   * @ref Sum() const.
   *
   * @param firstN Number of entries to be summed up.
   * @return Sum of the first @p firstN Vector entries.
   */
  double Sum(int firstN) const; // calculate sum of the first 'firstN' vector entries
  /**
   * @brief Calculate @p firstN (auto)correlation coefficients.
   * Autocorrelation coefficients are defined by:
   * \f{equation}{
   *     c_j = \frac{\mathrm{Cov}(v_i, v_{i+j})}{\mathrm{Var}(v)} = \frac{\mathrm{Cov}(v_i, v_{i+j})}{\sigma_v^2}
   *     \label{eq:Vector-autocorr_coefs}
   * \f}
   * where \f$\sigma_v^2\f$ is without the Bessel correction and:
   * \f{equation}{
   *     \mathrm{Cov}(v_i, v_{i+j}) = \frac{1}{N_j} \sum_{i=0}^{N_j-1} (v_i - \overline{v})(v_{i+j} - \overline{v})
   *     \label{eq:Vector-covariance}
   * \f}
   * where \f$\overline{v}\f$ is the mean value of Vector elements and
   * \f$N\f$ is the number of pairs taken into account (practically limited by the length of the Vector,
   * because \f$N\f$ can be only taken such that \f$N-1+j\f$ is still a valid index).
   * Note that \f$\mathrm{Var}(v)=\mathrm{Cov}(v_i, v_i)\f$.
   * The upper limit of \f$j\f$ is given by the parameter @p firstN (\f$j=\mathtt{firstN}-1\f$), which is set to @ref mSize if @p firstN = 0 is given.
   *
   * The calculation is carried out as follows:
   * 1. Calculate average value of the Vector elements (see @ref Mean() const)
   * 2. Calculate the sums in \f$\eqref{eq:Vector-covariance}\f$ for each \f$j\f$ counting the number of pairs contributing to it (\f$N_j\f$)
   * 3. Calculate the variance (\f$\sigma_v^2\f$) from the sum obtained in the previous step with \f$j=0\f$ dividing it by the number of pairs (which should be \f$N_0=\f$ @ref mSize)
   * 4. Calculate the (auto)correlation coefficients by dividing the sum (from step 2) by the number of pairs (\f$N_j\f$) and by the variance obtained in step 3.
   *
   * @param firstN Desired number of autocorrelation coefficients to be computed (default 0 means @p firstN = @ref mSize)
   * @return Vector of (auto)correlation coefficients (estimates).
   */
  Vector CorrelationCoefficients(int firstN = 0) const; // calculate the correlation coefficients
  /**
   * @brief Calculate the cumulative sum (running total) of the Vector elements.
   * Returns a Vector with the same size as the original ( @ref mSize) whose \f$j\f$-th element is:
   * \f{equation}{
   *     \sum_{i=0}^j v_i
   * \f}
   * where \f$v_i\f$ is the \f$i\f$-th element of the original Vector (which the method is called on).
   *
   * @par See also
   * @ref Sum() const, @ref Sum(int) const.
   *
   * @return Vector filled with cumulative sum.
   */
  Vector CumulativeSum() const; // calculate cumulative sum
  /**
   * @brief Calculate minimum value among Vector entries.
   *
   * @return Minimum value in the Vector.
   */
  double Min() const;
  /**
   * @brief Calculate maximum value among Vector entries.
   *
   * @return Maximum value in the Vector.
   */
  double Max() const;
  friend double CalculateScalarProduct(const Vector &a, const Vector &b); // scalar product
  friend double CalculateAngle(const Vector &a, const Vector &b);         // angle between
  friend int length(const Vector &vector);
  friend Vector CalculateVectorProduct(const Vector &a, const Vector &b); // cross-product
  friend Matrix CalculateOuterProduct(const Vector &a, const Vector &b);
  friend Vector CalculatePartialProduct(const Matrix &A, const Matrix &B, int row);
};

// prototype of friend and other functions
/**
 * @brief Get the length of the Vector @p vector.
 * Provides read-only acces to the private Vector member @ref Vector::mSize.
 *
 * @par See also
 * @ref Vector::GetSize() const.
 *
 * @param vector The instance of Vector whose length should be returned.
 * @return The lenght of the Vector @p vector.
 */
int length(const Vector &vector);
/**
 * @brief Calculate scalar product (dot product, inner product) of two Vectors
 * The scalar product of two Vectors @p a and @p b is given by:
 * \f{equation}{
 *     \mathbf{a}\cdot \mathbf{b} = \sum_{i=0}^{\mathtt{mSize}-1} a_i\cdot b_i
 * \f}
 * Assertion is made that the two Vectors are of the same length.
 * The result is a scalar value.
 *
 * @par See also
 * @ref CalculateAngle(const Vector&, const Vector&).
 *
 * @param a @e Bra-vector of the scalar product.
 * @param b @e Ket-vector of the scalar product.
 * @return Scalar product (dot, inner product) of the Vectors @p a and @p b.
 */
double CalculateScalarProduct(const Vector &a, const Vector &b);
/**
 * @brief Calculate angle between two Vectors.
 * The angle between Vectors @p a and @p b is calculated by:
 * \f{equation}{
 *     \alpha = \arccos \left(\frac{\mathbf{a}\cdot\mathbf{b}}{\left|\left|\mathbf{a}\right|\right|\cdot \left|\left|\mathbf{b}\right|\right|}\right)
 * \f}
 * where \f$\mathbf{a}\cdot\mathbf{b}\f$ is the scalar product of the two Vectors (see @ref CalculateScalarProduct(const Vector&, const Vector&))
 * and \f$\left|\left|\mathbf{a}\right|\right|\f$ is the 2-norm of Vector @p a (see @ref Vector::CalculateNorm(int) const).
 *
 * @par See also
 * @ref CalculateScalarProduct(const Vector&, const Vector&), @ref Vector::CalculateNorm(int) const.
 *
 * @param a First Vector.
 * @param b Second Vector.
 * @return The angle between the two Vectors (\f$\alpha\f$).
 */
double CalculateAngle(const Vector &a, const Vector &b);
/**
 * @relates Vector
 * @brief Overloading the '<<' operator to pretty-print the Vector @p vector
 * Pretty-printing of the Vector @p vector by `std::cout <<`
 * ```
 * [ a_0  a_1  a_2  ... ]
 * ```
 *
 * @param output (Reference to) output stream to which the printing is redirected.
 * @param vector Vector to be printed.
 * @return The (reference to the) same output stream to enable concatenating (with `<<`).
 */
std::ostream &operator<<(std::ostream &output, const Vector &vector);
/**
 * @relates Vector
 * @brief Calculate vector product (cross product) of two Vectors
 * The scalar product of two Vectors @p a and @p b is given by:
 * \f{equation}{
 *     \mathbf{a}\times \mathbf{b} = \left|
 *     \begin{array}
 *         \mathbf{i} & \mathbf{j} & \mathbf{k} \\
 *         a_0 & a_1 & a_2 \\
 *         b_0 & b_1 & b_2
 *     \end{array}
 *     \right| =
 *     \begin{bmatrix}
 *         a_1 b_2 - a_2 b_1 \\
 *         a_2 b_0 - a_0 b_2 \\
 *         a_0 b_1 - a_1 b_0
 *     \end{bmatrix}
 * \f}
 * Assertion is made that the two Vectors have length 3.
 * The result is a Vector of length 3.
 *
 * @par See also
 * @ref CalculateScalarProduct(const Vector&, const Vector&).
 *
 * @param a First Vector to be multiplied.
 * @param b Second Vector to be multiplied.
 * @return Vector product (cross product) of the Vectors @p a and @p b.
 */
Vector CalculateVectorProduct(const Vector &a, const Vector &b);

#endif
