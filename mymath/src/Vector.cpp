/*
 * A class to deal with vectors of double values in different manners
 * overloading enables vector summing, zero- and one-based indexing etc.
 * According to chap. 10 of Guide to Scientific Computing in Cpp
 * Author JJ, Date Dec 2019
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <ctime>
#include <sys/time.h>

#include "Vector.hpp"

// Overridden copy constructor
Vector::Vector(const Vector &otherVector)
{
  mSize = otherVector.GetSize();
  if ((mData = (double *)calloc(mSize, sizeof(double))) == NULL)
  {
    std::cerr << "Warning: Allocation of memory failed!!!\n";
  }
  int i;
  for (i = 0; i < mSize; i++)
  {
    mData[i] = otherVector.mData[i];
  }
}

// Constructor for a vector of given size
Vector::Vector(int size)
{
  mSize = size;
  if ((mData = (double *)calloc(mSize, sizeof(double))) == NULL)
  {
    std::cerr << "Warning: Allocation of memory failed!!!\n";
  }
  int i;
  for (i = 0; i < mSize; i++)
  {
    mData[i] = 0.0;
  }
}

// Read vector from file
Vector::Vector(FILE *fin)
{
  double auxfl;
  char line[50];
  int i;
  std::vector<double> auxvec;

  if (fin != NULL)
  {
    while (fscanf(fin, "%s\n", line) > 0)
    {
      if ((line[0] == '#') || (line[0] == '!') || (line[0] == '%'))
      {
        continue; // comment line
      }
      sscanf(line, "%lf", &auxfl);
      auxvec.push_back(auxfl);
      if (feof(fin))
      {
        break;
      }
    }
    mSize = (int)auxvec.size();
    if ((mData = (double *)calloc(mSize, sizeof(double))) == NULL)
    {
      std::cerr << "Warning: Allocation of memory failed!!!\n";
    }
    for (i = 0; i < mSize; i++)
    {
      mData[i] = auxvec.at(i);
    }
  }
  else
  {
    std::cerr << "ERROR: Pointer to file not valid!!!\n";
  }
}

// Read vector from file
Vector::Vector(const std::string name_of_file)
{
  double auxfl;
  std::vector<double> auxvec;
  std::string line;
  std::ifstream file;
  int i;

  file.open(name_of_file);

  if (file)
  {
    while (getline(file, line))
    {
      if ((line[0] == '#') || (line[0] == '!') || (line[0] == '%'))
      {
        continue; // comment line
      }
      sscanf(line.c_str(), "%lf", &auxfl);
      auxvec.push_back(auxfl);
    }
    mSize = (int)auxvec.size();
    if ((mData = (double *)calloc(mSize, sizeof(double))) == NULL)
    {
      std::cerr << "Warning: Allocation of memory failed!!!\n";
    }
    for (i = 0; i < mSize; i++)
    {
      mData[i] = auxvec.at(i);
    }
  }
  else
  {
    std::cerr << "ERROR: File name: " << name_of_file << " not valid (cannot open file)!!!\n";
  }
}

// New Vector from std::vector<double>
Vector::Vector(const std::vector<double> orig_vec)
{
  int i;
  mSize = orig_vec.size();
  if ((mData = (double *)calloc(mSize, sizeof(double))) == NULL)
  {
    std::cerr << "Warning: Allocation of memory failed!!!\n";
  }
  for (i = 0; i < mSize; i++)
  {
    mData[i] = orig_vec.at(i);
  }
}

// Overridden destructor
Vector::~Vector()
{
  free((void *)mData);
  mData = NULL;
}

// Public method to get size
int Vector::GetSize() const
{
  return mSize;
}

// C style zero-based indexing
double &Vector::operator[](int i)
{
  assert((i > -1) && (i < mSize));
  return mData[i];
}

// Overloading = operator
Vector &Vector::operator=(const Vector &otherVector)
{
  assert(mSize == otherVector.mSize);
  int i;
  for (i = 0; i < mSize; i++)
  {
    mData[i] = otherVector.mData[i];
  }
  return *this;
}

// Overloading + operator unary
Vector Vector::operator+() const
{
  Vector result(mSize);
  int i;
  for (i = 0; i < mSize; i++)
  {
    result.mData[i] = mData[i];
  }
  return result;
}

// Overloading - operator unary
Vector Vector::operator-() const
{
  Vector result(mSize);
  int i;
  for (i = 0; i < mSize; i++)
  {
    result.mData[i] = -mData[i];
  }
  return result;
}

// Overloading + operator binary
Vector Vector::operator+(const Vector &otherVector) const
{
  assert(mSize == otherVector.mSize);
  Vector result(mSize);
  int i;
  for (i = 0; i < mSize; i++)
  {
    result.mData[i] = mData[i] + otherVector.mData[i];
  }
  return result;
}

// Overloading - operator binary
Vector Vector::operator-(const Vector &otherVector) const
{
  assert(mSize == otherVector.mSize);
  Vector result(mSize);
  int i;
  for (i = 0; i < mSize; i++)
  {
    result.mData[i] = mData[i] - otherVector.mData[i];
  }
  return result;
}

// Overloading scalar multiplication
Vector Vector::operator*(const double scalar) const
{
  Vector result(mSize);
  int i;
  for (i = 0; i < mSize; i++)
  {
    result.mData[i] = mData[i] * scalar;
  }
  return result;
}

// Overloading *= operator
Vector &Vector::operator*=(const double scalar)
{
  int i;
  for (i = 0; i < mSize; i++)
  {
    mData[i] *= scalar;
  }
  return *this;
}

// Clear vector
int Vector::Clear()
{
  memset(mData, 0, mSize*sizeof(double)); // expect all zero bytes to represent 0.0
  assert(mData[0] == 0.0); // test if the above assumption is correct
  return 0;
}

// Calculate p-norm (default Eucleidian norm)
double Vector::CalculateNorm(int p) const
{
  double norm = 0;
  int i;
  for (i = 0; i < mSize; i++)
  {
    norm += pow(fabs(mData[i]), p);
  }
  return pow(norm, 1.0 / ((double)p));
}

// Print vector to file
int Vector::PrintToFile(FILE *fout) const
{
  int i;

  if (fout == NULL)
  {
    return -1;
  }
  for (i = 0; i < mSize; i++)
  {
    fprintf(fout, "%12f\n", mData[i]);
  }
  return 0;
}

// Mean value
double Vector::Mean() const
{
  int i;
  double sum = 0.0;

  for (i = 0; i < mSize; i++)
  {
    sum += mData[i];
  }
  return sum / mSize;
}

// print vector to the file
int Vector::PrintToFile(const std::string filename) const
{
  std::ofstream file;

  file.open(filename);
  if (!file)
  {
    std::cerr << "ERROR: cannot open file " << filename << "\n";
    return -1;
  }

  int i;
  timeval curr_time;
  gettimeofday(&curr_time, NULL);
  file << "# Vector printed by Vector::PrintToFile " << ctime((time_t *)&(curr_time.tv_sec));
  for (i = 0; i < mSize; i++)
  {
    file << std::setprecision(14) << mData[i] << "\n";
  }
  return 0;
}

// Mean value with standard error (Bessel corrected by default)
double Vector::MeanWithError(double &sigma, bool Bessel) const
{
  int i;
  double sum = 0.0, sumsq = 0.0;
  int N;

  for (i = 0; i < mSize; i++)
  {
    sum += mData[i];
    sumsq += mData[i] * mData[i];
  }

  if (Bessel)
  {
    N = mSize - 1;
  }
  else
  {
    N = mSize;
  }

  // sum /= mSize;
  // sumsq /=mSize;
  if (N < 1)
  {
    sigma = -1.0; // error, not enough data
  }
  else
  {
    sigma = sqrt((sumsq - sum * sum / mSize) / N / mSize);
  }

  return sum / mSize;
}

// calculate mean of correlated data and estimate its stderr (roughly) using the running avarage method
double Vector::MeanWithErrorEstimate(double &sigma) const
{
  double mean = Mean();
  Vector cumsum(mSize);
  cumsum = CumulativeSum();
  double max = mean, min = mean;
  int i;

  for (i = mSize / 2; i < mSize; i++)
  {
    max = std::max(max, cumsum(i) / (i + 1));
    min = std::min(min, cumsum(i) / (i + 1));
  }

  sigma = 0.6 * (max - min);

  return mean;
}

// mean by blocking method (for more info see Frenkel&Smit Appendix D or Allen&Tildesley p192(pdf 207))
double Vector::BlockMean(double &variance, int Nblocks) const
{
  double mean = Mean();
  Vector means(Nblocks);
  // double means_variance;
  int blocksize = (int)round((double)mSize / Nblocks);
  int i, j;

  for (i = 0; i < Nblocks - 1; i++)
  {
    for (j = i * blocksize; j < (i + 1) * blocksize; j++)
    {
      means(i) += mData[j];
    }
    means(i) /= blocksize;
  }
  for (j = i * blocksize; j < mSize; j++)
  {
    means(i) += mData[j];
  }
  means(i) /= j - i * blocksize;

  variance = means.Variance();
  // sigma = sqrt(means_variance/Nblocks);
  // std::cout << std::setw(14) << std::scientific << blocksize << "\t" << blocksize * variance / Variance() << "\n";

  return mean;
}

// calculate mean and estimate of variance of ensemble avarage (Flyvbjerg method described in F&S App. D)
double Vector::BlockMeanFlyvbjerg(double &variance, double &varianceError, int Nblocks) const
{
  double mean = Mean();
  Vector means(Nblocks);
  // double means_variance;
  int blocksize = (int)round((double)mSize / Nblocks);
  int i, j;

  for (i = 0; i < Nblocks - 1; i++)
  {
    for (j = i * blocksize; j < (i + 1) * blocksize; j++)
    {
      means(i) += mData[j];
    }
    means(i) /= blocksize;
  }
  for (j = i * blocksize; j < mSize; j++)
  {
    means(i) += mData[j];
  }
  means(i) /= j - i * blocksize;

  variance = means.Variance() / (Nblocks - 1.0);
  varianceError = variance * sqrt(2.0 / (Nblocks - 1.0));

  return mean;
}

// calculate the variance of the vector
double Vector::Variance(bool Bessel) const
{
  int i, N;
  double mean, sum = 0.0;

  mean = Mean();

  for (i = 0; i < mSize; i++)
  {
    sum += (mData[i] - mean) * (mData[i] - mean);
  }

  if (Bessel)
  {
    N = mSize - 1;
  }
  else
  {
    N = mSize;
  }

  return sum / N;
}

// calculate the variance of the vector enforcing the mean
double Vector::Variance(double mean, bool Bessel) const
{
  int i, N;
  double sum = 0.0;

  for (i = 0; i < mSize; i++)
  {
    sum += (mData[i] - mean) * (mData[i] - mean);
  }

  if (Bessel)
  {
    N = mSize - 1;
  }
  else
  {
    N = mSize;
  }

  return sum / N;
}

// calculate the sample variance (Bessel corrected)
double Vector::SampleVariance(bool Bessel) const
{
  int i, N;
  double mean, sum = 0.0;

  mean = Mean();

  for (i = 0; i < mSize; i++)
  {
    sum += (mData[i] - mean) * (mData[i] - mean);
  }

  if (Bessel)
  {
    N = mSize - 1;
  }
  else
  {
    N = mSize;
  }

  return sum / N;
}

// calculate the correlation coefficients
Vector Vector::CorrelationCoefficients(int firstN) const
{
  int i, j;
  double mean;
  int N;

  if (firstN <= 0)
  {
    N = mSize;
  }
  else
  {
    N = firstN;
  }

  Vector res(N);
  int k[N];

  mean = Mean();

  k[0] = 0; // only for make not to print warning (may be used uninitialized)
  for (i = 1; i < N; i++)
  {
    k[i] = 0;
  }

  for (i = 0; i < mSize; i++)
  {
    for (j = i; j < std::min(mSize, i + N); j++)
    {
      res(j - i) += (mData[i] - mean) * (mData[j] - mean);
      k[j - i]++;
    }
  }
  res(0) /= ((k[0] == 0) ? 1.0 : k[0]);
  for (i = 1; i < N; i++)
  {
    res(i) /= ((k[i] == 0) ? 1.0 : k[i]);
    // std::cout << "res before: " << res(i) << ", res after:";
    res(i) /= res(0);
    // std::cout << res(i) << "\n";
  }

  res(0) = 1.0;

  return res;
}

// Excess kurtosis
double Vector::ExcessKurtosis(bool Bessel) const
{
  long double fourthMoment = 0.0L;
  double mean = Mean();
  int i;
  double secMoment = Variance() * mSize;

  for (i = 0; i < mSize; i++)
  {
    fourthMoment += (long double)powl((long double)(mData[i] - mean), 4.0L);
  }

  if (Bessel)
  { // sumetc -4 -2 -0 | tabproc "C*A/B/B-3" "C*(C+1)*(C-1)/((C-2)*(C-3))*A/B/B-3*(C-1)*(C-1)/((C-2)*(C-3))
    // "Sample (excess) kurtosis" from Harding et al. 2014
    return (double)mSize * (mSize + 1.0) * (mSize - 1.0) / ((mSize - 2.0) * (mSize - 3.0)) * fourthMoment / secMoment / secMoment - 3.0 * (mSize - 1.0) * (mSize - 1.0) / ((mSize - 2.0) * (mSize - 3.0));
  }
  else
  {
    // "Kurtosis" simple
    return (double)mSize * fourthMoment / secMoment / secMoment - 3.0;
  }
}

// Excess kurtosis with error estimate from blocking method
double Vector::ExcessKurtosisWithError(double &sigma, int Nblocks, bool Bessel) const
{
  long double fourthMoment = 0.0L;
  double mean = Mean();
  Vector kurtoses(Nblocks);
  // double means_variance;
  int blocksize = (int)round((double)mSize / Nblocks);
  long double secMoment;
  int i, j;
  double kurtosis;

  for (i = 0; i < Nblocks - 1; i++)
  {
    fourthMoment = 0.0L;
    secMoment = 0.0L;
    for (j = i * blocksize; j < (i + 1) * blocksize; j++)
    {
      secMoment += powl((long double)(mData[j] - mean), 2.0L);
      fourthMoment += powl((long double)(mData[j] - mean), 4.0L);
    }
    if (Bessel)
    {
      kurtoses(i) = (double)blocksize * (blocksize + 1.0) * (blocksize - 1.0) / ((blocksize - 2.0) * (blocksize - 3.0)) * fourthMoment / secMoment / secMoment - 3.0 * (blocksize - 1.0) * (blocksize - 1.0) / ((blocksize - 2.0) * (blocksize - 3.0));
    }
    else
    {
      kurtoses(i) = (double)blocksize * fourthMoment / secMoment / secMoment - 3.0;
    }
  }
  fourthMoment = 0.0L;
  secMoment = 0.0L;
  for (j = i * blocksize; j < mSize; j++)
  {
    secMoment += powl((long double)(mData[j] - mean), 2.0L);
    fourthMoment += powl((long double)(mData[j] - mean), 4.0L);
  }
  blocksize = mSize - i * blocksize;
  if (Bessel)
  {
    kurtoses(i) = (double)blocksize * (blocksize + 1.0) * (blocksize - 1.0) / ((blocksize - 2.0) * (blocksize - 3.0)) * fourthMoment / secMoment / secMoment - 3.0 * (blocksize - 1.0) * (blocksize - 1.0) / ((blocksize - 2.0) * (blocksize - 3.0));
  }
  else
  {
    kurtoses(i) = (double)blocksize * fourthMoment / secMoment / secMoment - 3.0;
  }

  kurtosis = kurtoses.MeanWithError(sigma);
  // std::cerr << "kurtoses: " << kurtoses << "\n";

  return kurtosis;
}

// Overloading print << operator
std::ostream &operator<<(std::ostream &output, const Vector &vector)
{
  int i;
  output << "[ ";
  for (i = 0; i < vector.GetSize() - 1; i++)
  {
    output << vector.Read(i) << "    ";
  }
  output << vector.Read(vector.GetSize() - 1) << " ]";
  return output;
}

// Friend function - scalar product
double CalculateScalarProduct(const Vector &a, const Vector &b)
{
  assert(a.mSize == b.mSize);
  int i;
  double result = 0.0;
  for (i = 0; i < a.mSize; i++)
  {
    result += a.mData[i] * b.mData[i];
  }
  return result;
}

// Friend function - Matlab style length of vector
int length(const Vector &vector)
{
  return vector.mSize;
}

// Sum of vector entries
double Vector::Sum() const
{
  int i;
  double sum = 0.0;
  for (i = 0; i < mSize; i++)
  {
    sum += mData[i];
  }
  return sum;
}

// Sum of the first 'firstN' vector entries
double Vector::Sum(int firstN) const
{
  assert((firstN > 0) && (firstN <= mSize));
  int i;
  double sum = 0.0;
  for (i = 0; i < firstN; i++)
  {
    sum += mData[i];
  }
  return sum;
}

// Cumulative sum of the entries
Vector Vector::CumulativeSum() const
{
  Vector result(mSize);
  int i;
  double sum = 0.0;

  for (i = 0; i < mSize; i++)
  {
    sum += mData[i];
    result(i) = sum;
  }

  return result;
}

// min and max
double Vector::Min() const
{
  double min = mData[0];
  int i;
  for (i = 1; i < mSize; i++)
  {
    min = fmin(min, mData[i]);
  }
  return min;
}

double Vector::Max() const
{
  double max = mData[0];
  int i;
  for (i = 1; i < mSize; i++)
  {
    max = fmax(max, mData[i]);
  }
  return max;
}

// Friend function - angle
double CalculateAngle(const Vector &a, const Vector &b)
{
  assert(a.mSize == b.mSize);
  double alpha = 0.0;

  alpha = acos(CalculateScalarProduct(a, b) / (a.CalculateNorm() * b.CalculateNorm()));
  return alpha;
}

// Variance with error from block analysis
double Vector::VarianceWithError(double &sigma, int Nblocks) const
{
  double mean = Mean();
  double variance;
  Vector variances(Nblocks); // variances of each block
  int blocksize = (int)ceil((double)mSize / Nblocks);
  int i, j;

  for (i = 0; i < Nblocks - 1; i++)
  {
    for (j = i * blocksize; j < (i + 1) * blocksize; j++)
    {
      variances(i) += pow(mData[j] - mean, 2.0);
    }
    variances(i) /= blocksize;
  }
  for (j = i * blocksize; j < mSize; j++)
  {
    variances(i) += pow(mData[j] - mean, 2.0);
  }
  variances(i) /= j - i * blocksize;
  // std::cerr << variances << "\n";
  variances.MeanWithError(sigma, true);
  variance = Variance(mean);
  return variance;
}

/**
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
Vector CalculateVectorProduct(const Vector &a, const Vector &b)
{
  assert(a.mSize == b.mSize);
  assert(a.mSize == 3);
  Vector result(3);

  result[0] = a.Read(1) * b.Read(2) - a.Read(2) * b.Read(1);
  result[1] = a.Read(2) * b.Read(0) - a.Read(0) * b.Read(2);
  result[2] = a.Read(0) * b.Read(1) - a.Read(1) * b.Read(0);

  return result;
}