#ifndef NEWTONMETHODHEADERFILE
#define NEWTONMETHODHEADERFILE

class NewtonMethod
{
  private:
    int mSize; // number of equations
    int mNMaxSteps; // maximum number of steps
    double mTolerance;
    Vector *mpx; // solution vector
    Vector (*mpF)(Vector&); // system of equtions to solve (pointer to function) 
    NewtonMethod(const NewtonMethod& NM) {}; // Keeping copy constructor private disables its use...
  public:
    NewtonMethod(Vector (*pF)(Vector&), Vector& x, const int nMaxSteps = 20, const double tol = 1e-18);
//    NewtonMethod(int size, Vector (*pF)(Vector&), Vector& x, const int nMaxSteps = 20, const double tol = 1e-18);
    ~NewtonMethod();
    Vector Solve() const;
};

#endif