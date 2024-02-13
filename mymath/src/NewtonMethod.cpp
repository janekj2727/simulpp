/*
* A class to solve nonlinear system by Newton method
* Derivatives to Jacobian are replaced by differences to be general.
* Author JJ, Date Nov 2020
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
#include "NewtonMethod.hpp"

// #define DEBUG

// Default constructor
NewtonMethod::NewtonMethod(Vector (*pF)(Vector&), Vector& x, const int nMaxSteps, const double tol)
{
  int size = x.GetSize();
  mSize = size;
  mNMaxSteps = nMaxSteps;
  mTolerance = tol;
  mpx = new Vector(x);
  mpF = pF;
}

/*
// Constructor specifying size
NewtonMethod::NewtonMethod(int size, Vector (*pF)(Vector&), Vector& x, const int nMaxSteps, const double tol)
{
  mSize = size;
  mNMaxSteps = nMaxSteps;
  mTolerance = tol;
  mpx = new Vector(x);
  mpF = pF;
}
*/

// Destructor - free matrix and vector memory
NewtonMethod::~NewtonMethod()
{
  delete mpx;
  mpF = NULL;
}

// Solve the system of equations using Newton method
Vector NewtonMethod::Solve() const
{
  Vector rx(*mpx);

  int i, j, k;
  double delta = 1e-5;
  Vector rF(mSize);
  Vector rb(mSize);
  Vector rh(mSize);
  Vector rplus(mSize);
  Vector rminus(mSize);
  Vector rdiff(mSize);
  Matrix rJ(mSize, mSize);
  LinearSystem *LS;

  for (i = 0; i < mNMaxSteps; i++)
  {
      rF = -mpF(rx);
    //   #ifdef DEBUG
    //       std::cout << "|F| = " << rF.CalculateNorm() << "\n";
    //       std::cout << "-F(x):\n" << rF << "\n\n";
    //   #endif
      if (rF.CalculateNorm() < mTolerance)
      {
          break;
      }
      
      for (j = 0; j < mSize; j++)
      {
          rh(j) = delta;
          rplus = rx + rh;
          rminus = rx - rh;
          rdiff = (mpF(rplus) - mpF(rminus))*(1/(2*delta));
          for (k = 0; k < mSize; k++)
          {
              rJ(k, j) = rdiff(k);
          }
          rh(j) = 0.0;
      }
    //   #ifdef DEBUG
    //       std::cout << "\nJacobian:\n" << rJ << "\n";
    //   #endif
      LS = new LinearSystem(rJ, rF);
      rb = LS->SolveGJ();
      delete LS;
      
    //   #ifdef DEBUG
    //       std::cout << "|b| = " << rb.CalculateNorm() << "\n";
    //   #endif
      if (rb.CalculateNorm() < mTolerance)
      {
          break;
      }
      rx = rx + rb;
    //   #ifdef DEBUG
    //       std::cout << "rb: " << rb << "\n";
    //       std::cout << "rx: " << rx << "\n";
    //       std::cout << "niter: " << i << "\n";
    //   #endif
  }
  
  return rx;
}
