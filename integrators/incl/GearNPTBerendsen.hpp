#ifndef GEARNPTBERENDSENHEADER
#define GEARNPTBERENDSENHEADER

#include "AbstractIntegrator.hpp"
#include "AbstractGearIntegrator.hpp"
class Matrix;
class Vector;
class Atom;

class GearNPTBerendsen : public AbstractGearIntegrator
{
private:
    int Integrate(int Nsteps) override;
    inline int CalculateRHS(Atom *at) override;
    // thermostat variables and functions
    double Tf; // final (desired) temperature
    double tau_T; // thermostat constant
    double lambda_T; // velocity scaling due to thermostat
    double Tinst; // instanteneous temperature
    inline double GetTempScaling();
    // barostat variables and functions
    double Pf; // final (desired) pressure
    double tau_P; // barostat constant
    bool molecular_based; // scaling is molecular/atomic based (NPTBerendsen implemented with molecular scaling -> true)
    double lambda_P; // predicted value of box size
    inline double GetScalingFactor(int order);
    inline double GetShakeScaling() {return 1.0;};

public:
    GearNPTBerendsen(double step, SimulatedSystem *simsystem, int Rsize, char *integName, double Tfinal, double tauT, double Pfinal, double tauP, Vector trvpB, bool startWithShake = true);
    ~GearNPTBerendsen();    
};

#endif