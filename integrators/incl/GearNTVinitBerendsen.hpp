#ifndef GEARNTVINITBERENDSENHEADER
#define GEARNTVINITBERENDSENHEADER

#include "AbstractIntegrator.hpp"
#include "AbstractGearIntegrator.hpp"
class Matrix;
class Vector;

class GearNTVinitBerendsen : public AbstractGearIntegrator
{
private:
    double scaling_factor; // scaling of box in one step...
    int scaling_steps;     // number of steps in which scaling will be performed
    int Integrate(int Nsteps) override;
    inline int CalculateRHS(Atom *at) override;
    // thermostat variables and functions
    double Tf;       // final (desired) temperature
    double tau_T;    // thermostat constant
    double lambda_T; // velocity scaling due to thermostat
    double Tinst;    // instanteneous temperature
    inline double GetTempScaling();
    // barostat variables and functions
    double Pf;            // final (desired) pressure
    double tau_P;         // barostat constant
    bool molecular_based; // scaling is molecular/atomic based (NTVinitBerendsen implemented with molecular scaling -> true)
    double lambda_P;      // predicted value of box size
    inline double GetScalingFactor(int order);
    inline double GetShakeScaling() {return 1.0;};

public:
    GearNTVinitBerendsen(double step, SimulatedSystem *simsystem, int Rsize, char *integName, double Tfinal, double tauT, double rhofinal, double taurho, Vector trvpB, bool startWithShake = true);
    ~GearNTVinitBerendsen();
};

#endif