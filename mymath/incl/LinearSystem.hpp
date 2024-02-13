#ifndef LINEARSYSTEMHEADERFILE
#define LINEARSYSTEMHEADERFILE

class LinearSystem
{
  private:
    int mSize; // linear system size
    Matrix *mpA; // matrix
    Vector *mpb; // rhs vector
    LinearSystem(const LinearSystem& LS) {}; // Keeping copy constructor private disables its use...
  public:
    LinearSystem(const Matrix& A, const Vector& b);
    ~LinearSystem();
    Vector SolveGJ() const;
    Vector SolveC() const;
};

#endif