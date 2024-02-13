#ifndef ABSTRACTGEARINTEGRATORHEADER
#define ABSTRACTGEARINTEGRATORHEADER

#include "AbstractIntegrator.hpp"
#ifdef PARALLEL
extern int thread_count;
#include <omp.h>
#endif

class AbstractGearIntegrator : public AbstractIntegrator
{
protected:
    Matrix *A;           // predictor
    Vector *r;           // corrector
    const int size;      // size of atom->R
    double sumr;         // "sum" or corrector coefficients for SHAKE
    double lambda_Shake; // scaling factor for SHAKE
    double omega_Shake;  // SHAKE overrelaxation
    double eps_Shake;    // maximum error for SHAKE
    int size_Shake;      // SHAKE "type" (size not including TRVP)
    int maxiter_Shake;   // maximum number of iterations for SHAKE

protected:
    virtual int Integrate(int Nsteps) = 0;
    int Initialize(bool startWithShake = true) override;
    int PrepareConfig() override;
    int Correction(Molecule *mol);
    double ShakeHook(Molecule *mol);
    virtual int CalculateRHS(Atom *at) = 0;
    virtual double EnergyOfExtraDOF() const override { return 0.0; };

public:
    AbstractGearIntegrator(double step, SimulatedSystem *simsystem, int Rsize, char *integName, Vector trvpB, bool startWithShake = true);
    ~AbstractGearIntegrator();
};

#endif