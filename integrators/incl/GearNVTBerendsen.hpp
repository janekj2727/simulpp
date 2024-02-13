#ifndef GEARNVTBERENDSENHEADER
#define GEARNVTBERENDSENHEADER

#include "AbstractIntegrator.hpp"
#include "AbstractGearIntegrator.hpp"
class Matrix;
class Vector;

class GearNVTBerendsen : public AbstractGearIntegrator
{
private:
    int Integrate(int Nsteps) override;
    inline int CalculateRHS(Atom *at) override;
    // thermostat variables and functions
    double Tf;       // final (desired) temperature
    double tau_T;    // thermostat constant
    double lambda_T; // velocity scaling due to thermostat
    double Tinst;    // instanteneous temperature
    inline double GetTempScaling();

public:
    GearNVTBerendsen(double step, SimulatedSystem *simsystem, int Rsize, char *integName, double Tfinal, double tauT, Vector trvpB, bool startWithShake = true);
    ~GearNVTBerendsen();
};

#endif