#ifndef ABSTRACTINTEGRATORHEADER
#define ABSTRACTINTEGRATORHEADER

#ifdef PARALLEL
extern int thread_count;
#include <omp.h>
#endif

class SimulatedSystem;

// #define STARTWITHSHAKE  // uncomment if SHAKE is to be performed during initialization (breaks exact restart, but improves stability)

class AbstractIntegrator
{
private:
    AbstractIntegrator(): predVel(1) {};

protected:
    std::vector<int> positionsAt; // indices of positions in R (to be moved when whole molecule is moved)
    // double epsShake; // maximum relative error of constr. bonds
    const int predVel; // position where the predicted velocity is stored (default 1 – predicted velocity = corrected velocity (for Gear))

public:
    double h;                // integration step
    SimulatedSystem *system; // pointer to simulated system of molecules

    AbstractIntegrator(double step, SimulatedSystem *simsystem, int predV)
    : predVel(predV)
    {
        h = step;
        system = simsystem;
        positionsAt.clear();
    };
    virtual ~AbstractIntegrator(){};
    virtual int Integrate(int Nsteps) = 0;
    virtual int Initialize(bool startWithShake = true) = 0;
    virtual int PrepareConfig() = 0;
    virtual double EnergyOfExtraDOF() const = 0;

    int SetSystem(SimulatedSystem *newsystem)
    {
        system = newsystem;
        return 0;
    }
    virtual int RestoreExtendedDOFs()
    {
        return 0;
    }
};

#endif